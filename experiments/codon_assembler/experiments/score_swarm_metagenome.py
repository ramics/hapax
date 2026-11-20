"""Score every assembler's contigs on one multi-swarm metagenome run.

Contig-only scoring, so hapax and the incumbents are judged by the same rules.
Every contig of at least MIN_LEN nt is scored; shorter ones are dropped for
all tools alike.

Metrics, and why each exists:

* **Window entity assignment.** Each contig is cut into 250-nt windows (step
  125) and each window is mapped to the full truth set (every haplotype of
  every swarm plus every background genome). A window belongs to the entity
  (swarm, or background genome) of its best hit; windows whose best hits in two
  entities score within 1% of each other are ambiguous and ignored, which
  exempts genuinely shared sequence. A contig's entity is its majority entity.
* **Cross-entity chimaera**: a contig where a second entity owns at least 10%
  of the assigned windows and at least 2 windows. Swarms are separate hosts,
  so any such join is an error. The two HIV env hosts, the two HBV hosts and
  HXB2-versus-subtype-C are the sister pairs that make this non-trivial.
* **Contigs per swarm / completeness**: contigs whose entity is the swarm, and
  the longest single contig's aligned span on its best haplotype divided by
  that haplotype's length. hapax's design target is 1 contig at 1.0.
* **Haplotypes resolved**: a haplotype h counts if some contig's best full
  alignment among the swarm's haplotypes is h, covers at least 95% of h, and
  has divergence below min(1%, half the distance from h to its nearest other
  haplotype). The half-distance rule means a consensus that merely sits near a
  cluster of close haplotypes does not get credit for all of them. This is
  what haplotype assemblers (Haploflow) try to do and hapax does not.
* **Variant recall / precision**: truth is the abundance-weighted minor-allele
  frequency at each column of the swarm's MSA; true sites have MAF >= 5%.
  A tool's calls are sites with second-allele frequency >= 5% at depth >= 20,
  read from (a) hapax's own variants.tsv, and (b) for every tool, a pileup of
  all reads mapped back to that tool's contigs (minimap2 -x sr, primary
  alignments). (b) is the standard assemble-then-map pipeline. Calls are
  lifted to MSA columns through an alignment of the contig to the swarm's
  majority consensus. A call is a true positive if its column has any real
  variation (MAF >= 1%). Recall counts true sites not covered by any contig as
  missed.

    python experiments/score_swarm_metagenome.py results/swarm_metagenome/work/seed1
"""

from __future__ import annotations

import collections
import itertools
import json
import math
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENV = Path(os.environ.get("HAPAX_TOOLS", "/home/user/mamba/envs/asm/bin"))
MAFFT = Path(os.environ.get("HAPAX_MAFFT", "/home/user/mamba/envs/aux/bin/mafft"))
MIN_LEN = 500
WIN, STEP = 250, 125
MAF_TRUE, MAF_ANY, CALL_FREQ, CALL_DEPTH = 0.05, 0.01, 0.05, 20
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def read_fasta(path) -> list[tuple[str, str]]:
    out, name, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out.append((name, "".join(seq).upper()))
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name:
        out.append((name, "".join(seq).upper()))
    return out


def paf(cmd: list[str]) -> list[list[str]]:
    r = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return [l.split("\t") for l in r.stdout.splitlines() if l]


def tag(row, key):
    for t in row[12:]:
        if t.startswith(key + ":"):
            return t.split(":", 2)[2]
    return None


# ---------------------------------------------------------------- truth

def load_truth(work: Path):
    truth = json.loads((work / "truth.json").read_text())
    entity = {}
    for sw, haps in truth["swarms"].items():
        for h in haps:
            entity[h] = sw
    for acc in truth["background"]:
        entity[acc] = acc
    allfa = work / "truth_all.fa"
    if not allfa.exists():
        with open(allfa, "w") as f:
            for n in entity:
                for nm, s in read_fasta(work / "sources" / f"{n}.fa"):
                    f.write(f">{nm}\n{s}\n")
    return truth, entity, allfa


def swarm_msa(work: Path, swarm: str, haps: dict):
    """MSA of a swarm's haplotypes (cached), majority consensus, and per-column truth MAF."""
    cache = work.parent / "msa"
    cache.mkdir(exist_ok=True)
    fa = cache / f"{swarm}.fa"
    aln = cache / f"{swarm}.aln.fa"
    if not aln.exists():
        with open(fa, "w") as f:
            for h in haps:
                (n, s), = read_fasta(work / "sources" / f"{h}.fa")
                f.write(f">{n}\n{s}\n")
        r = subprocess.run([str(MAFFT), "--quiet", "--auto", "--thread", "4", str(fa)], capture_output=True, text=True, check=True)
        aln.write_text(r.stdout)
    rows = {n: s for n, s in read_fasta(aln)}
    L = len(next(iter(rows.values())))
    # leading/trailing gaps are amplicon extent, not deletions: treat as missing
    ext = {}
    for n, s in rows.items():
        a = len(s) - len(s.lstrip("-"))
        b = len(s.rstrip("-"))
        ext[n] = (a, b)
    cons, cons_col, maf = [], [], []
    for c in range(L):
        w = collections.Counter()
        tot = 0.0
        for n, s in rows.items():
            if not (ext[n][0] <= c < ext[n][1]):
                continue
            f = haps[n]["frac"]
            tot += f
            w[s[c]] += f
        if tot == 0:
            continue
        base, top = max(w.items(), key=lambda x: x[1])
        if base == "-":
            continue
        cons.append(base)
        cons_col.append(c)
        bases = {k: v for k, v in w.items() if k != "-"}
        bt = sum(bases.values())
        srt = sorted(bases.values(), reverse=True)
        maf.append(srt[1] / bt if len(srt) > 1 else 0.0)
    return "".join(cons), maf


def nearest_neighbour_div(work: Path, swarm: str, haps: dict) -> dict:
    rows = {n: s for n, s in read_fasta(work.parent / "msa" / f"{swarm}.aln.fa")}
    names = list(rows)
    nn = {n: 1.0 for n in names}
    for a, b in itertools.combinations(names, 2):
        sa, sb = rows[a], rows[b]
        d = n = 0
        for x, y in zip(sa, sb):
            if x == "-" or y == "-":
                continue
            n += 1
            d += x != y
        dv = d / max(n, 1)
        nn[a] = min(nn[a], dv)
        nn[b] = min(nn[b], dv)
    return nn


# ---------------------------------------------------------------- per-contig assignment

def assign_windows(contigs, allfa: Path, entity: dict, tmp: Path, threads=4):
    wfa = tmp / "windows.fa"
    with open(wfa, "w") as f:
        for n, s in contigs:
            starts = list(range(0, max(len(s) - WIN, 0) + 1, STEP))
            if starts[-1] + WIN < len(s):
                starts.append(len(s) - WIN)
            for st in starts:
                f.write(f">{n}@{st}\n{s[st:st + WIN]}\n")
    rows = paf([str(ENV / "minimap2"), "-c", "-x", "asm20", "-N", "50", "--secondary=yes", "-p", "0.5", "-t", str(threads), str(allfa), str(wfa)])
    best = collections.defaultdict(dict)  # window -> entity -> best identity*cov
    for r in rows:
        q, qlen, qs, qe, t = r[0], int(r[1]), int(r[2]), int(r[3]), r[5]
        alen = int(r[10])
        nm = int(tag(r, "NM") or 0)
        score = (alen - nm) / qlen  # matches over full window: penalises partial hits
        e = entity[t]
        if score > best[q].get(e, 0):
            best[q][e] = score
    per = collections.defaultdict(collections.Counter)
    nwin = collections.Counter()
    for n, s in contigs:
        starts = list(range(0, max(len(s) - WIN, 0) + 1, STEP))
        if starts[-1] + WIN < len(s):
            starts.append(len(s) - WIN)
        nwin[n] = len(starts)
        for st in starts:
            hits = best.get(f"{n}@{st}")
            if not hits:
                per[n]["_unmapped"] += 1
                continue
            srt = sorted(hits.items(), key=lambda x: -x[1])
            if srt[0][1] < 0.80:
                per[n]["_unmapped"] += 1
            elif len(srt) > 1 and srt[0][1] - srt[1][1] < 0.01:
                per[n]["_ambiguous"] += 1
            else:
                per[n][srt[0][0]] += 1
    out = {}
    for n, s in contigs:
        c = per[n]
        ents = [(e, k) for e, k in c.most_common() if not e.startswith("_")]
        assigned = sum(k for _, k in ents)
        major = ents[0][0] if ents else None
        chim = len(ents) > 1 and ents[1][1] >= 2 and ents[1][1] >= 0.10 * assigned
        out[n] = {"len": len(s), "entity": major, "windows": nwin[n], "assigned": assigned,
                  "unmapped": c["_unmapped"], "ambiguous": c["_ambiguous"],
                  "chimaera": chim, "entities": dict(ents)}
    return out


def full_alignments(contigs_fa: Path, ref_fa: Path, threads=4):
    rows = paf([str(ENV / "minimap2"), "-c", "--cs", "-x", "asm20", "-N", "100", "--secondary=yes", "-p", "0.3", "-t", str(threads), str(ref_fa), str(contigs_fa)])
    return rows


# ---------------------------------------------------------------- variants

def cs_qmap(row) -> dict:
    """Map contig (query, forward coordinates) positions to target positions from a cs tag."""
    qlen, qs, qe, strand, ts = int(row[1]), int(row[2]), int(row[3]), row[4], int(row[7])
    cs = tag(row, "cs")
    m = {}
    q = qs if strand == "+" else qlen - qe  # position in the (possibly reverse-complemented) query
    t = ts
    i = 0
    while i < len(cs):
        op = cs[i]
        j = i + 1
        while j < len(cs) and cs[j] not in ":*+-~=":
            j += 1
        val = cs[i + 1:j]
        if op == ":":
            for _ in range(int(val)):
                m[q] = t
                q += 1
                t += 1
        elif op == "*":
            m[q] = t
            q += 1
            t += 1
        elif op == "+":
            q += len(val)
        elif op == "-":
            t += len(val)
        i = j
    if strand == "-":
        m = {qlen - 1 - k: v for k, v in m.items()}
    return m


def pileup_calls(bam: Path, contigs: set[str], tmp: Path) -> dict:
    """Per-contig {pos0: (depth, maf)} from samtools mpileup, same rule for every tool."""
    if not contigs:
        return {}
    bed = tmp / "regions.txt"
    bed.write_text("\n".join(contigs) + "\n")
    out = collections.defaultdict(dict)
    for ctg in contigs:
        r = subprocess.run([str(ENV / "samtools"), "mpileup", "-B", "-Q", "13", "-q", "0", "-d", "100000", "-r", ctg, str(bam)],
                           capture_output=True, text=True)
        for line in r.stdout.splitlines():
            f = line.split("\t")
            pos, ref, bases = int(f[1]) - 1, f[2].upper(), f[4]
            cnt = collections.Counter()
            k = 0
            while k < len(bases):
                ch = bases[k]
                if ch == "^":
                    k += 2
                    continue
                if ch in "+-":
                    j = k + 1
                    while j < len(bases) and bases[j].isdigit():
                        j += 1
                    k = j + int(bases[k + 1:j])
                    continue
                if ch in ".,":
                    cnt[ref] += 1
                elif ch.upper() in "ACGT":
                    cnt[ch.upper()] += 1
                k += 1
            d = sum(cnt.values())
            if d >= CALL_DEPTH:
                srt = sorted(cnt.values(), reverse=True)
                out[ctg][pos] = (d, srt[1] / d if len(srt) > 1 else 0.0)
    return out


def hapax_calls(var_tsv: Path) -> dict:
    out = collections.defaultdict(dict)
    if not var_tsv.exists():
        return out
    with open(var_tsv) as f:
        next(f)
        for line in f:
            x = line.rstrip("\n").split("\t")
            ctg, pos = x[0], int(x[1]) - 1
            counts = sorted([int(v) for v in x[6:10]], reverse=True)
            d = sum(counts)
            if d >= CALL_DEPTH:
                out[ctg][pos] = (d, counts[1] / d)
    return out


def variant_scores(calls: dict, swarm_ctgs: dict, cons_maps: dict, maf: list) -> dict:
    """calls: contig -> pos -> (depth, maf). cons_maps: contig -> {cpos: cons_pos}."""
    true_sites = {i for i, m in enumerate(maf) if m >= MAF_TRUE}
    covered = set()
    detected = {}
    tp = fp = 0
    for ctg in swarm_ctgs:
        cm = cons_maps.get(ctg, {})
        covered.update(cm.values())
        for p, (d, f) in calls.get(ctg, {}).items():
            if f < CALL_FREQ:
                continue
            t = cm.get(p)
            if t is None:
                continue  # call outside the aligned part of the contig: not scorable
            if maf[t] >= MAF_ANY:
                tp += 1
                if t in true_sites:
                    detected[t] = max(detected.get(t, 0), f)
            else:
                fp += 1
    pairs = [(maf[t], f) for t, f in detected.items()]
    corr = None
    if len(pairs) > 2:
        xs, ys = zip(*pairs)
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
        sy = math.sqrt(sum((y - my) ** 2 for y in ys))
        corr = sum((x - mx) * (y - my) for x, y in pairs) / (sx * sy) if sx and sy else None
    return {"true_sites": len(true_sites), "true_covered": len(true_sites & covered), "detected": len(detected),
            "recall": len(detected) / len(true_sites) if true_sites else None,
            "precision": tp / (tp + fp) if tp + fp else None, "calls_tp": tp, "calls_fp": fp, "freq_r": corr}


# ---------------------------------------------------------------- main

def score_tool(name: str, ctg_path: Path, work: Path, truth, entity, allfa, msas, nns, var_tsv=None) -> dict:
    tmp = work / "score" / name
    tmp.mkdir(parents=True, exist_ok=True)
    contigs = [(f"c{i}", s) for i, (n, s) in enumerate(read_fasta(ctg_path)) if len(s) >= MIN_LEN] if ctg_path.exists() else []
    orig = {f"c{i}": n for i, (n, s) in enumerate(read_fasta(ctg_path)) if len(s) >= MIN_LEN} if ctg_path.exists() else {}
    if not contigs:
        return {"n_contigs": 0}
    cfa = tmp / "contigs.fa"
    with open(cfa, "w") as f:
        for n, s in contigs:
            f.write(f">{n}\n{s}\n")
    asg = assign_windows(contigs, allfa, entity, tmp)
    full = full_alignments(cfa, allfa)
    res = {"n_contigs": len(contigs), "total_len": sum(len(s) for _, s in contigs),
           "chimaeras": sum(1 for a in asg.values() if a["chimaera"]),
           "chimaera_detail": {orig[n]: a["entities"] for n, a in asg.items() if a["chimaera"]},
           "unassigned_contigs": sum(1 for a in asg.values() if a["entity"] is None)}
    # background recovery, for context: fraction of each genome's length in its own contigs (capped)
    res["background"] = {}
    for acc, b in truth["background"].items():
        L = sum(a["len"] for a in asg.values() if a["entity"] == acc)
        res["background"][acc] = {"contigs": sum(1 for a in asg.values() if a["entity"] == acc),
                                  "len_frac": min(L / b["len"], 1.0)}
    # mapping reads back for the pileup calls
    bam = tmp / "reads.bam"
    p1, p2 = work / "reads_1.fq", work / "reads_2.fq"
    pile = {}
    if p1.exists():
        mm = subprocess.Popen([str(ENV / "minimap2"), "-a", "-x", "sr", "-t", "4", "--secondary=no", str(cfa), str(p1), str(p2)],
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        subprocess.run([str(ENV / "samtools"), "sort", "-@", "2", "-o", str(bam), "-"], stdin=mm.stdout, check=True, stderr=subprocess.DEVNULL)
        mm.wait()
        subprocess.run([str(ENV / "samtools"), "index", str(bam)], check=True)
        swarm_contigs = {n for n, a in asg.items() if a["entity"] in truth["swarms"]}
        pile = pileup_calls(bam, swarm_contigs, tmp)
        bam.unlink()
        Path(str(bam) + ".bai").unlink()
    hcalls = {}
    if var_tsv is not None:
        raw = hapax_calls(var_tsv)
        hcalls = {n: raw.get(orig[n], {}) for n, _ in contigs}
    res["swarms"] = {}
    for sw, haps in truth["swarms"].items():
        mine = {n: a for n, a in asg.items() if a["entity"] == sw}
        cons, maf = msas[sw]
        # contig -> best haplotype alignment
        hap_names = set(haps)
        best_hap = {}
        for r in full:
            q, t = r[0], r[5]
            if q not in mine or t not in hap_names:
                continue
            alen, nm = int(r[10]), int(tag(r, "NM") or 0)
            matches = alen - nm
            if q not in best_hap or matches > best_hap[q][1]:
                best_hap[q] = (t, matches, nm / alen, (int(r[8]) - int(r[7])) / int(r[6]))
        resolved = set()
        for q, (t, m, div, cov) in best_hap.items():
            if cov >= 0.95 and div < min(0.01, nns[sw][t] / 2):
                resolved.add(t)
        compl = max((cov for (_, _, _, cov) in best_hap.values()), default=0.0)
        # contig -> consensus coordinate maps
        cons_fa = tmp / f"{sw}.cons.fa"
        cons_fa.write_text(f">cons\n{cons}\n")
        sub = tmp / f"{sw}.ctg.fa"
        with open(sub, "w") as f:
            for n, s in contigs:
                if n in mine:
                    f.write(f">{n}\n{s}\n")
        cmaps = {}
        if mine:
            for r in paf([str(ENV / "minimap2"), "-c", "--cs", "-x", "asm20", "-t", "4", str(cons_fa), str(sub)]):
                if tag(r, "tp") != "P":
                    continue
                cmaps.setdefault(r[0], {}).update(cs_qmap(r))
        abund = sorted(haps.items(), key=lambda x: -x[1]["frac"])
        top = [h for h, _ in abund[: max(1, len(abund) // 4)]]
        entry = {"haplotypes": len(haps), "contigs": len(mine), "completeness": compl,
                 "chimaeric_contigs": sum(1 for a in mine.values() if a["chimaera"]),
                 "resolved": len(resolved), "resolved_top_quartile": len(resolved & set(top)),
                 "top_quartile": len(top),
                 "variants_remap": variant_scores(pile, mine, cmaps, maf)}
        if var_tsv is not None:
            entry["variants_own"] = variant_scores(hcalls, mine, cmaps, maf)
        res["swarms"][sw] = entry
    return res


def main() -> None:
    work = Path(sys.argv[1])
    truth, entity, allfa = load_truth(work)
    runs = json.loads((work / "runs.json").read_text())
    msas = {sw: swarm_msa(work, sw, haps) for sw, haps in truth["swarms"].items()}
    nns = {sw: nearest_neighbour_div(work, sw, haps) for sw, haps in truth["swarms"].items()}
    scores = {}
    for name, r in runs.items():
        var = Path(r["contigs"]).parent / "variants.tsv" if name.startswith("hapax") else None
        scores[name] = score_tool(name, Path(r["contigs"]), work, truth, entity, allfa, msas, nns, var)
        scores[name]["wall"] = r["wall"]
        scores[name]["status"] = r["status"]
        scores[name]["max_rss_gb"] = r.get("max_rss_gb")
        if "stats" in r:
            scores[name]["hmm_fraction"] = r["stats"].get("hmm_fraction")
        print(name, json.dumps({k: v for k, v in scores[name].items() if k not in ("swarms", "background", "chimaera_detail")}))
    (work / "scores.json").write_text(json.dumps(scores, indent=1))


if __name__ == "__main__":
    main()
