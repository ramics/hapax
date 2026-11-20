#!/usr/bin/env python3
"""Classify every base where a hapax contig disagrees with a comparator contig.

Why this exists. On the Berlin wastewater astrovirus (results/ww_berlin/REPORT.md) hapax's
contigs were 97.0% identical to the authors' SPAdes contig while MEGAHIT's were 99.2%, a gap
wider than the 1.6% polymorphism rate measured from the reads. A per-base identity number cannot
say whether that is consensus error, a different haplotype, or misaligned reads, so this script
asks the reads directly at every disagreeing site:

  * reads -> hapax contig (minimap2 -ax sr, independent of hapax's own aligner): is the base
    hapax emitted the majority of the reads that map there? If not, hapax's consensus call is
    wrong *by its own reads*, whatever the comparator says.
  * reads -> comparator contig: is the comparator's base the reads' majority, and is the site
    polymorphic (minor allele >= 20%)? A disagreement at a polymorphic site where hapax took the
    minor allele is a haplotype choice, not an error.

Each disagreement is classified as substitution / insertion / deletion, with its distance to
the nearer hapax contig end, and the read depth there.

With `--labels` (read name -> lineage, e.g. each read's best bait contig), every substitution is
also asked which lineage's reads carry hapax's base and which carry the comparator's. That is the
test that separates "consensus voted wrongly" from "the component held several lineages and the
vote followed whichever was locally commoner".

Deterministic; `--seed` is accepted for the experiment convention and unused.
"""
import argparse
import collections
import json
import os
import re
import subprocess
import sys
import time


def read_fa(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def pileup_counts(bam, ref_fa):
    """Per (contig, 0-based pos): Counter of A/C/G/T/'-'(deletion)/'+'(any insertion after)."""
    out = collections.defaultdict(collections.Counter)
    p = subprocess.run(["samtools", "mpileup", "-f", ref_fa, "-d", "100000", "-Q", "0", "-B", "-A", bam],
                       capture_output=True, text=True, check=True)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        ctg, pos, ref, s = f[0], int(f[1]) - 1, f[2].upper(), f[4]
        c = out[(ctg, pos)]
        i = 0
        while i < len(s):
            ch = s[i]
            if ch == "^":
                i += 2
                continue
            if ch == "$":
                i += 1
                continue
            if ch in "+-":
                m = re.match(r"[+-](\d+)", s[i:])
                n = int(m.group(1))
                if ch == "+":
                    c["+"] += 1
                i += len(m.group(0)) + n
                continue
            if ch in ".,":
                c[ref] += 1
            elif ch == "*":
                c["-"] += 1
            elif ch.upper() in "ACGT":
                c[ch.upper()] += 1
            i += 1
    return out


def map_reads(ref, reads, bam, threads):
    cmd = (f"minimap2 -ax sr -t {threads} {ref} {reads} 2>/dev/null | "
           f"samtools view -b -F 0x904 - | samtools sort -o {bam} - && samtools index {bam}")
    subprocess.run(cmd, shell=True, check=True)


def majority(c):
    bases = {b: n for b, n in c.items() if b in "ACGT-"}
    if not bases:
        return None, 0, 0
    tot = sum(bases.values())
    b = max(bases, key=bases.get)
    return b, bases[b], tot


def alleles_by_lineage(bam, ctg, pos, lab):
    """Base at contig position `pos` (0-based) in each read of `bam`, grouped by read label."""
    out = collections.defaultdict(collections.Counter)
    p = subprocess.run(["samtools", "view", bam, f"{ctg}:{pos + 1}-{pos + 1}"], capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        ri, qi, base = int(f[3]) - 1, 0, None
        for n, op in re.findall(r"(\d+)([MIDNSHX=])", f[5]):
            n = int(n)
            if op in "M=X":
                if ri <= pos < ri + n:
                    base = f[9][qi + pos - ri]
                    break
                ri += n
                qi += n
            elif op in "IS":
                qi += n
            elif op in "DN":
                if ri <= pos < ri + n:
                    base = "-"
                    break
                ri += n
        if base:
            out[lab.get(f[0], "?")][base] += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True, help="hapax contigs.fa")
    ap.add_argument("--target", required=True, help="comparator contig FASTA")
    ap.add_argument("--reads", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--labels", help="TSV read_name<TAB>lineage; enables per-lineage allele attribution")
    ap.add_argument("--target-lineage", help="label of the comparator's own lineage in --labels")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    hap = read_fa(a.contigs)
    tgt = read_fa(a.target)

    # 1. disagreements, from minimap2 --cs blocks (hapax contig = query, target = ref)
    paf = subprocess.run(["minimap2", "-c", "--cs=long", "-x", "asm10", a.target, a.contigs],
                         capture_output=True, text=True, check=True).stdout
    sites = []
    aligned_cols = 0
    for line in paf.splitlines():
        f = line.split("\t")
        q, qlen, qs, qe, strand, t, ts = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7])
        cs = [x for x in f[12:] if x.startswith("cs:Z:")][0][5:]
        # walk in target coords; query coords run forward (+) or backward (-)
        tpos = ts
        qpos = qs if strand == "+" else qe  # for '-', query position counts down
        for op, body in re.findall(r"([=:*+\-~])([A-Za-z0-9]+(?:[a-z]*)?)", cs):
            if op == "=":
                n = len(body)
                aligned_cols += n
                tpos += n
                qpos += n if strand == "+" else -n
            elif op == "*":
                aligned_cols += 1
                qp = qpos if strand == "+" else qpos - 1
                sites.append(dict(kind="sub", ctg=q, qpos=qp, qlen=qlen, tctg=t, tpos=tpos,
                                  tbase=body[0].upper(), qbase=body[1].upper(), strand=strand, n=1))
                tpos += 1
                qpos += 1 if strand == "+" else -1
            elif op == "+":  # bases in hapax not in target
                n = len(body)
                qp = qpos if strand == "+" else qpos - n
                sites.append(dict(kind="ins", ctg=q, qpos=qp, qlen=qlen, tctg=t, tpos=tpos, n=n,
                                  strand=strand, seq=body.upper()))
                qpos += n if strand == "+" else -n
            elif op == "-":  # bases in target missing from hapax
                n = len(body)
                qp = qpos if strand == "+" else qpos - 1
                sites.append(dict(kind="del", ctg=q, qpos=qp, qlen=qlen, tctg=t, tpos=tpos, n=n,
                                  strand=strand, seq=body.upper()))
                tpos += n

    # 2. reads -> hapax contigs, and reads -> target, by an aligner that is not hapax's
    hbam, tbam = f"{a.out}/reads_vs_hapax.bam", f"{a.out}/reads_vs_target.bam"
    map_reads(a.contigs, a.reads, hbam, a.threads)
    map_reads(a.target, a.reads, tbam, a.threads)
    hp = pileup_counts(hbam, a.contigs)
    tp = pileup_counts(tbam, a.target)

    # 3. whole-contig check: how often does hapax's own base lose the read vote?
    whole = collections.Counter()
    nonmaj_positions = []
    for name, s in hap.items():
        for i, b in enumerate(s):
            c = hp.get((name, i))
            if not c:
                whole["no_reads"] += 1
                continue
            m, n, tot = majority(c)
            whole["covered"] += 1
            if m is not None and m != b and c.get(b, 0) < n:
                whole["hapax_not_majority"] += 1
                nonmaj_positions.append((name, i, b, dict(c)))
            ins = c.get("+", 0)
            if tot and 2 * ins > tot:
                whole["majority_insertion_missing"] += 1

    # 4. classify each disagreement
    rows = []
    for st in sites:
        qp = st["qpos"]
        c = hp.get((st["ctg"], qp), collections.Counter())
        m, n, tot = majority(c)
        tc = tp.get((st["tctg"], st["tpos"]), collections.Counter())
        tm, tn, ttot = majority(tc)
        minor = 0.0
        if ttot:
            vals = sorted((v for k, v in tc.items() if k in "ACGT-"), reverse=True)
            minor = (vals[1] / ttot) if len(vals) > 1 else 0.0
        hb = hap[st["ctg"]][qp] if 0 <= qp < len(hap[st["ctg"]]) else "?"
        end_dist = min(qp, st["qlen"] - 1 - qp)
        if st["kind"] == "sub":
            # hapax base in target orientation is st['qbase']; compare in hapax orientation
            verdict = ("hapax_not_majority" if m and m != hb and c.get(hb, 0) < n else
                       "target_minor_allele" if tm and tm != st["tbase"] else
                       "polymorphic_hapax_took_other" if minor >= 0.2 else
                       "hapax_majority_target_majority_differ")
        else:
            verdict = "indel"
        rows.append(dict(st, hdepth=tot, hmaj=m, hbase=hb, hcounts=dict(c), tdepth=ttot, tmaj=tm,
                         tminor=round(minor, 3), tcounts=dict(tc), end_dist=end_dist, verdict=verdict))

    lineage = None
    if a.labels:
        lab = dict(l.split()[:2] for l in open(a.labels) if l.strip())
        lineage = collections.Counter()
        for r in rows:
            if r["kind"] != "sub":
                continue
            by = alleles_by_lineage(hbam, r["ctg"], r["qpos"], lab)
            own = by.get(a.target_lineage, collections.Counter())
            own_maj = own.most_common(1)[0][0] if own else None
            others = collections.Counter()
            for k, v in by.items():
                if k != a.target_lineage:
                    others.update(v)
            r["own_lineage"] = dict(own)
            r["other_lineages"] = {k: dict(v) for k, v in by.items() if k != a.target_lineage}
            if own_maj is None:
                lineage["no_target_lineage_reads"] += 1
            elif own_maj == r["hbase"]:
                lineage["hapax_base_is_target_lineage_majority"] += 1
            elif others.get(r["hbase"], 0) > own.get(r["hbase"], 0):
                lineage["hapax_base_from_other_lineage"] += 1
            else:
                lineage["hapax_base_unexplained"] += 1

    summ = dict(
        lineage_attribution=lineage,
        aligned_cols=aligned_cols,
        disagreements=collections.Counter(r["kind"] for r in rows),
        disagreement_bases=sum(r["n"] for r in rows),
        verdicts=collections.Counter(r["verdict"] for r in rows),
        sub_within_100_of_end=sum(1 for r in rows if r["kind"] == "sub" and r["end_dist"] < 100),
        indel_lengths=collections.Counter(r["n"] for r in rows if r["kind"] != "sub"),
        indel_len_mod3=collections.Counter(r["n"] % 3 for r in rows if r["kind"] != "sub"),
        whole_contig=whole,
        seconds=round(time.time() - t0, 1),
    )
    with open(f"{a.out}/sites.tsv", "w") as fh:
        keys = ["kind", "n", "ctg", "qpos", "qlen", "end_dist", "strand", "tpos", "tbase", "qbase", "hbase",
                "hmaj", "hdepth", "hcounts", "tmaj", "tdepth", "tminor", "tcounts", "verdict",
                "own_lineage", "other_lineages"]
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
    with open(f"{a.out}/nonmajority.tsv", "w") as fh:
        for name, i, b, c in nonmaj_positions:
            fh.write(f"{name}\t{i}\t{b}\t{c}\n")
    json.dump(dict(summary=summ, args=vars(a), git_sha=subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()),
        open(f"{a.out}/manifest.json", "w"), indent=1, default=dict)
    print(json.dumps(summ, indent=1, default=dict))


if __name__ == "__main__":
    sys.exit(main())
