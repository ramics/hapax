#!/usr/bin/env python3
"""What assembly path selection costs, measured on a WHOLE library against ground truth.

THE QUESTION. `results/pair_ladder/` left this project with a specific claim: below ~3.8%
divergence every assembler tested fuses lineages, and PenguiN's raw output holds clean
multi-kb pieces of distinct lineages among mosaics - so the reads carry the information and
*path selection* is what is missing. That claim has only ever been tested on simulated
two-lineage mixtures. This script tests it on a real library where the answer is known.

WHY THE WHOLE LIBRARY AND NOT THE RESIDUE (Imo, 2026-11-22). The obvious experiment is to
re-assemble the residue at different path settings and see whether it assembles better. It
is the wrong experiment: residue membership is decided by the IDENTIFICATION step that ran
on a previous assembly, so a residue-only arm bakes in that assembly's path choices. Worse,
the failure mode we are hunting is invisible to it: under laxer settings residue reads can be
absorbed into an abundant known contig rather than forming contigs of their own, which a
residue-only assembly cannot see happen. So every arm here assembles the WHOLE library and
anything about the residue is an OUTPUT of the arm, never an input to it.

WHY ABSORPTION IS THE FAILURE MODE AND NOT AN AFTERTHOUGHT. A setting that shrinks the
residue by fusing novel sequence into a known contig reads as an improvement on every
summary number this project computes. This repo has been caught by that shape twice:
subtracting unidentified contigs inflated removal to 85.9% against a true 77.0%
(`results/subtract_known/`), and every assembler tested fuses lineages below ~3.8%
divergence (`results/pair_ladder/`). So residue size never ships alone here. Each arm is
scored first against the ground truth below, and only then on how much residue it leaves.

THE GROUND TRUTH, WHICH NOBODY CHOSE. Library SRR23421544 (treated wastewater, Russia,
collected 2022-05-10, 4,290,429 single-end reads) is the library the post-cutoff holdout
comes from, and the depositors assembled 25 virus genomes out of these same waters and
deposited them on 2023-07-31. `results/holdout_candidates/depth.tsv` records, by mapping the
reads back, that all 25 are present in THIS run at mean depths from 0.4x to 700x, and
`results/holdout_candidates/shortlist.tsv` records each one's divergence from public RefSeq.
That gives three things no simulation can:

  - a recovery panel at 25 real depths (0.4x to 700x) and real distances from public
    reference, 0.08% to 30.4% divergent;
  - an ABSENT arm, which is the control and it is free: 53 more genomes in the same
    committed panel have exactly 0 reads mapped in this run. Any contig that aligns to one
    of those is a false join, measured in the same pass by the same code.

WHAT THE PANEL IS, MEASURED RATHER THAN ASSUMED, AFTER GETTING IT WRONG ONCE. The
divergences in `shortlist.tsv` are each genome's distance to its nearest *public* relative,
not to each other, so the panel's internal structure has to be measured. An all-against-all
pass over the 78 genomes finds **at least six near-identical pairs**, several co-present in
this library at high depth: OQ722321/OQ722330 99.98% identical over 6,063 nt at 499x and
700x, OQ722318/OQ722329 99.97%, OQ722324/OQ722335 99.91%, OQ722323/OQ722327 99.84%,
OQ722320/OQ722333 99.95%, OQ722317/OQ722328 99.52%, OQ722319/OQ722331 99.44%. That is
inside the band where `results/pair_ladder/` says every assembler fuses lineages, in a real
library, which the first version of this file reported as the opposite because
`--secondary=no` had hidden every alignment but each genome's self-hit (see `paf`).

CONSEQUENCE FOR WHAT CAN BE CLAIMED. At 99.9% identity a 200 nt alignment at >=95% identity
cannot say WHICH member of a pair a contig came from, so per-accession breadth is reported
with the pair structure beside it and the question "did both lineages come out separately"
is NOT answered here - it needs variant-level attribution, not alignment breadth. A contig
spanning two DIFFERENT groups is a join between unrelated molecules, which is the chimaera
number that ships beside every contiguity figure.

DIRECTION OF EVERY ARM, STATED BEFORE THE RUN. `--merge-level`'s similarity floor and
`--bubble-level 0` keep divergent paths apart: expect MORE and SHORTER contigs, less fusion,
possibly less breadth. A lower similarity floor merges more: expect higher N50, more
mosaics, and a residue that shrinks for the wrong reason - that is the absorption arm and it
is in the ladder deliberately. `--prune-level`/`--prune-depth`/`--low-local-ratio` govern
whether thin branches survive: keeping them should assemble more low-depth material and
more junk, pruning harder should grow the residue. `--no-mercy`/`--no-local` remove the two
mechanisms that rescue low-coverage paths, so they are the "less path rescue" end.

NOT IN THE LADDER, AND WHY. `--min-count` and the k-list are the two levers already tested
in this project and both moved nothing, so they are held FIXED at the committed baseline's
values for every arm. A ladder that varied them would not be measuring path selection.

THE CONTROL THAT MAKES THE LADDER READABLE. `baseline_repeat` is the baseline setting run
again with a different thread count. MEGAHIT is deterministic given its inputs, but that is
an assumption and this ladder's whole content is differences between arms - so the
same-setting difference is measured rather than assumed. If an arm differs from baseline by
less than baseline differs from its own repeat, the arm moved nothing.

SEEDS. One arm is one deterministic pass over the whole library, so there is no spread to
report on the assembly or recovery numbers and five of them would be five identical runs.
`--seed` is accepted and recorded because the convention requires it and because the
`cascade` stage, which samples read blocks, does use it. Where this script reports a
number with no spread, the number has no spread by construction, and that is stated rather
than hidden behind a single-seed directory.
"""
from __future__ import annotations

import argparse
import json
import re
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: The ladder. Every arm inherits BASE and adds its own flags, so the two dead levers
#: (--min-count, --k-list) are at MEGAHIT's defaults in every arm and never vary.
#:
#: WHAT THE COMMITTED BASELINE ACTUALLY WAS, AND THE GAP THIS CLOSES. The assembly behind
#: `results/holdout_cascade/` has no recorded command line: its manifest says
#: `assembled_at_run_time: false` and points at an `assembly.json` that does not exist, so
#: only the contig names survive (`megahit_k141_*`, i.e. a k sweep ending at 141, with the
#: shortest contigs ~310 nt). That is consistent with MEGAHIT's defaults plus a 200 nt
#: floor and NOT distinguishable from `--presets meta-sensitive`, which also ends at 141.
#: So the baseline here is MEGAHIT's defaults with `--min-contig-len 200`, this script
#: records the full argv and version, and `meta_sensitive` is carried as its own arm rather
#: than assumed - it changes `--min-count` to 1, which is one of the dead levers, so it
#: cannot sit in BASE without varying the thing the ladder holds fixed.
BASE = ["--min-contig-len", "200"]
ARMS: dict[str, dict] = {
    "baseline": {
        "flags": [],
        "direction": "the committed setting; everything else is read against it",
    },
    "baseline_repeat": {
        "flags": [],
        "threads_override": 2,
        "direction": "same setting, different thread count: the same-instrument "
                     "difference, which bounds what an arm has to beat to mean anything",
    },
    "strict_bubble": {
        "flags": ["--merge-level", "20,0.98"],
        "direction": "merge only near-identical bubbles. Expect more, shorter contigs and "
                     "less lineage fusion; breadth may fall",
    },
    "no_bubble": {
        "flags": ["--bubble-level", "0"],
        "direction": "no bubble merging at all. The extreme of keeping paths apart",
    },
    "lax_bubble": {
        "flags": ["--merge-level", "20,0.90"],
        "direction": "merge bubbles down to 90% similarity. THE ABSORPTION ARM: expect "
                     "longer contigs, more mosaics, and a residue that shrinks wrongly",
    },
    "keep_thin": {
        "flags": ["--prune-level", "0", "--prune-depth", "1"],
        "direction": "keep low-depth unitigs. Expect more low-depth contigs and more junk",
    },
    "keep_local": {
        "flags": ["--low-local-ratio", "0.05", "--disconnect-ratio", "0.02"],
        "direction": "keep low-coverage local branches that the defaults remove",
    },
    "prune_hard": {
        "flags": ["--prune-level", "3"],
        "direction": "prune harder. Expect a larger residue and little change deep down",
    },
    "no_rescue": {
        "flags": ["--no-mercy", "--no-local"],
        "direction": "remove both low-coverage rescue mechanisms. The least-path-rescue end",
    },
    "meta_sensitive": {
        "flags": ["--presets", "meta-sensitive"],
        "direction": "the preset the residue assembly used, never compared to the default "
                     "until now. NOT a path arm: it also sets --min-count 1, a lever this "
                     "project already measured as dead, so it is reported separately and "
                     "no path conclusion rests on it",
    },
}


def log(*a) -> None:
    print(*a, flush=True)


def sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(script: str) -> dict:
    """Script hash beside the git SHA, per CLAUDE.md: a -dirty SHA names HEAD, not the code."""
    try:
        sha = subprocess.run(["git", "-C", str(ROOT), "describe", "--always", "--dirty"],
                             capture_output=True, text=True).stdout.strip()
    except Exception:
        sha = "unknown"
    return {"git_sha": sha, "script_sha256": sha256(Path(script)),
            "python": sys.version.split()[0]}


def peak_rss_bytes() -> int:
    r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return int(r) * 1024


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def contig_stats(fa: Path) -> dict:
    lens = []
    for line in open(fa):
        if line.startswith(">"):
            m = re.search(r"len=(\d+)", line)
            lens.append(int(m.group(1)) if m else 0)
        elif not lens:
            continue
    if not lens or lens[0] == 0:  # header without len=: measure the sequence
        lens = [len(s) for s in read_fasta(fa).values()]
    lens.sort(reverse=True)
    total = sum(lens)
    acc, n50 = 0, 0
    for L in lens:
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    return {"contigs": len(lens), "bases": total, "n50": n50, "max": lens[0] if lens else 0,
            "ge500": sum(1 for L in lens if L >= 500),
            "ge1000": sum(1 for L in lens if L >= 1000),
            "ge3000": sum(1 for L in lens if L >= 3000)}


# ----------------------------------------------------------------- assemble
def assemble(args, arm: str) -> dict:
    spec = ARMS[arm]
    out = args.work / f"asm_{arm}"
    if (out / "final.contigs.fa").exists() and not args.force:
        log(f"[{arm}] assembly present, keeping it")
    else:
        if out.exists():
            shutil.rmtree(out)
        threads = spec.get("threads_override", args.threads)
        reads_args = (["-r", str(args.reads)] if args.r2 is None else
                      ["-1", str(args.reads), "-2", str(args.r2)])
        cmd = [str(args.megahit), *reads_args, "-o", str(out),
               "-t", str(threads), "--memory", str(args.memory), *BASE, *spec["flags"]]
        log(f"[{arm}] {' '.join(cmd)}")
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        secs = time.perf_counter() - t0
        (args.work / f"megahit_{arm}.log").write_text(proc.stdout + proc.stderr)
        if proc.returncode != 0:
            raise SystemExit(f"[{arm}] MEGAHIT exited {proc.returncode}; see "
                             f"{args.work / f'megahit_{arm}.log'}")
        log(f"[{arm}] assembled in {secs / 60:.1f} min")
        (out / "seconds.txt").write_text(str(secs))
    secs = float((out / "seconds.txt").read_text()) if (out / "seconds.txt").exists() else None
    st = contig_stats(out / "final.contigs.fa")
    return {"arm": arm, "flags": BASE + spec["flags"], "direction": spec["direction"],
            "threads": spec.get("threads_override", args.threads),
            "seconds": secs, "peak_rss_bytes": peak_rss_bytes(), **st,
            "contigs_fa": str(out / "final.contigs.fa")}


# ----------------------------------------------------------------- truth panel
def build_panel(args) -> tuple[Path, dict]:
    """The present and absent arms of the committed candidate panel for THIS library.

    Present = mapped_reads > 0 in `depth.tsv` for this run, which was measured by mapping
    the reads back rather than asserted. Absent = exactly 0 mapped reads in the same table,
    same panel, same code - so the false-join control costs nothing and cannot be accused of
    being drawn from a different distribution of sequence.
    """
    rows = [l.rstrip("\n").split("\t") for l in open(args.depth)]
    head, rows = rows[0], rows[1:]
    i = {k: head.index(k) for k in ("run", "genome", "length", "mapped_reads", "mean_depth")}
    present, absent = {}, {}
    for r in rows:
        if r[i["run"]] != args.run:
            continue
        (present if float(r[i["mapped_reads"]]) > 0 else absent)[r[i["genome"]]] = {
            "length": int(r[i["length"]]), "mapped_reads": float(r[i["mapped_reads"]]),
            "mean_depth": float(r[i["mean_depth"]])}
    seqs = read_fasta(args.candidates)
    div = {}
    if args.shortlist.exists():
        sl = [l.rstrip("\n").split("\t") for l in open(args.shortlist)]
        h = sl[0]
        for r in sl[1:]:
            d = dict(zip(h, r))
            if d.get("divergence_pct"):
                try:
                    div[d["acc"]] = float(d["divergence_pct"])
                except ValueError:
                    pass
    missing = [a for a in list(present) + list(absent) if a not in seqs]
    if missing:
        raise SystemExit(f"panel accessions absent from {args.candidates}: {missing[:5]}")
    panel = args.work / "panel.fa"
    with open(panel, "w") as w:
        for acc in list(present) + list(absent):
            w.write(f">{acc}\n{seqs[acc]}\n")
    meta = {"run": args.run, "present": present, "absent": absent,
            "divergence_pct": {a: div.get(a) for a in list(present) + list(absent)}}
    return panel, meta


def paf(args, query: Path, target: Path, preset: str) -> list[dict]:
    """Align with SECONDARY ALIGNMENTS ON, which is not the obvious setting and matters.

    THIS IS A DEFECT THIS FILE HAD AND THE NUMBERS IT PRODUCED WERE WRONG IN A WAY THAT
    LOOKED LIKE A FINDING. With minimap2's `--secondary=no`, each contig reports only its
    single best panel genome. The panel contains near-identical pairs - OQ722321/OQ722330
    are 99.98% identical over 6,063 nt and both sit in this library at 499x and 700x - so
    one member of each pair was credited and the other read 0.00 breadth. That is an
    artefact of the flag, and read as data it says "a 499x genome failed to assemble",
    which is exactly the kind of sentence this project has had to retract before. With
    secondaries on, a contig is credited against every panel genome it matches.

    The same flag broke the all-against-all pass over the panel itself: with self-hits as
    the primary alignment, every other hit was secondary and dropped, so the panel appeared
    to contain no two related genomes when in fact it contains at least six near-identical
    pairs.
    """
    cmd = [str(args.minimap2), "-t", str(args.threads), "-cx", preset,
           "--secondary=yes", "-N", "50", "-p", "0.01", str(target), str(query)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"minimap2 failed: {proc.stderr[-2000:]}")
    out = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        tags = dict(t.split(":", 2)[::2] for t in f[12:] if ":" in t)
        out.append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]),
                    "strand": f[4], "t": f[5], "tlen": int(f[6]), "ts": int(f[7]),
                    "te": int(f[8]), "matches": int(f[9]), "blocklen": int(f[10]),
                    "mapq": int(f[11]),
                    "identity": int(f[9]) / max(1, int(f[10]))})
    return out


def merge(spans: list[tuple[int, int]]) -> int:
    spans = sorted(spans)
    tot, cur_s, cur_e = 0, None, None
    for s, e in spans:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot


def lineage_groups(args, panel: Path) -> dict[str, str]:
    """Group panel genomes that are near-identical to each other.

    A contig covering two genomes from the SAME group is not a mosaic - those genomes are
    0.07-1.9% apart and a contig cannot be said to have joined two things that are the same
    thing to within sequencing error. A contig covering two DIFFERENT groups is a join
    between molecules that really are distinct, and that is the mosaic number. Groups are
    computed from the panel itself (single-linkage at >=95% identity over >=50% of the
    shorter genome), so the definition is in the data rather than in a hand-made list.
    """
    seqs = read_fasta(panel)
    hits = paf(args, panel, panel, "asm20")
    parent = {a: a for a in seqs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for h in hits:
        if h["q"] == h["t"]:
            continue
        short = min(len(seqs[h["q"]]), len(seqs[h["t"]]))
        if h["identity"] >= 0.95 and h["blocklen"] >= 0.5 * short:
            a, b = find(h["q"]), find(h["t"])
            if a != b:
                parent[a] = b
    return {a: find(a) for a in seqs}


def score(args, arm: str, contigs: Path, panel: Path, meta: dict,
          groups: dict[str, str]) -> dict:
    """Recovery, purity and mosaics for one arm, all from one contig-against-panel pass.

    PURITY IS THE ABSORPTION NUMBER. For the contig that carries most of a genome, the
    fraction of the CONTIG explained by that genome. A contig that is 96% the genome is a
    recovery; a 20 kb contig with the genome in 4 kb of it is the genome absorbed into
    something else, and both read as "recovered" on breadth alone.
    """
    hits = [h for h in paf(args, contigs, panel, "asm10")
            if h["identity"] >= args.min_identity and h["blocklen"] >= args.min_aln]
    by_genome: dict[str, list[dict]] = {}
    by_contig: dict[str, list[dict]] = {}
    for h in hits:
        by_genome.setdefault(h["t"], []).append(h)
        by_contig.setdefault(h["q"], []).append(h)

    genomes = {}
    for acc, info in {**meta["present"], **meta["absent"]}.items():
        hs = by_genome.get(acc, [])
        covered = merge([(h["ts"], h["te"]) for h in hs])
        best = max(hs, key=lambda h: h["te"] - h["ts"], default=None)
        row = {"length": info["length"], "mean_depth": info["mean_depth"],
               "present": acc in meta["present"],
               "divergence_pct": meta["divergence_pct"].get(acc),
               "covered_nt": covered, "covered_frac": covered / info["length"],
               "contigs_touching": len(hs),
               "contigs_for_90pct": None, "best_contig": None,
               "best_contig_len": None, "best_contig_covers_genome": None,
               "best_contig_purity": None, "best_contig_purity_accession_only": None,
               "group": groups.get(acc), "group_size": None}
        if best is not None:
            own = merge([(h["qs"], h["qe"]) for h in by_contig[best["q"]] if h["t"] == acc])
            # PURITY IS AGAINST THE GROUP, NOT THE ACCESSION, and the difference is not
            # cosmetic. OQ722323 is a 5,048 nt PARTIAL deposit 99.84% identical to the
            # 9,467 nt OQ722327; a contig carrying the full molecule explains only 0.53 of
            # itself with OQ722323's record and 0.92 with the group's, so per-accession
            # purity reports a correct recovery of the longer relative as absorption into
            # foreign sequence. The accession-level figure is kept beside it, because where
            # a genome has no near-identical relative the two are the same number.
            grp = {a for a, g in groups.items() if g == groups[acc]}
            own_grp = merge([(h["qs"], h["qe"]) for h in by_contig[best["q"]]
                             if h["t"] in grp])
            row.update({"best_contig": best["q"], "best_contig_len": best["qlen"],
                        "best_contig_covers_genome": (best["te"] - best["ts"]) / info["length"],
                        "best_contig_purity": own_grp / best["qlen"],
                        "best_contig_purity_accession_only": own / best["qlen"],
                        "group": groups[acc], "group_size": len(grp)})
            # How many contigs it takes to cover 90% of the genome, greedily largest
            # first. None means 90% is not reached by any number of them.
            order = sorted(hs, key=lambda h: h["te"] - h["ts"], reverse=True)
            row["contigs_for_90pct"] = None
            for need in range(1, len(order) + 1):
                if merge([(x["ts"], x["te"]) for x in order[:need]]) >= 0.9 * info["length"]:
                    row["contigs_for_90pct"] = need
                    break
        genomes[acc] = row

    mosaics = []
    for q, hs in by_contig.items():
        gs = {groups[h["t"]] for h in hs if (h["qe"] - h["qs"]) >= args.min_mosaic}
        if len(gs) > 1:
            mosaics.append({"contig": q, "qlen": hs[0]["qlen"], "groups": sorted(gs),
                            "spans": sorted((h["t"], h["qs"], h["qe"]) for h in hs)})
    false_joins = [acc for acc, r in genomes.items()
                   if not r["present"] and r["covered_nt"] >= args.min_aln]
    present_rows = [r for r in genomes.values() if r["present"]]
    return {"arm": arm,
            "genomes": genomes,
            "panel_present": len(meta["present"]), "panel_absent": len(meta["absent"]),
            "recovered_90pct": sum(1 for r in present_rows if r["covered_frac"] >= 0.9),
            "recovered_50pct": sum(1 for r in present_rows if r["covered_frac"] >= 0.5),
            "mean_covered_frac": sum(r["covered_frac"] for r in present_rows) / len(present_rows),
            "median_purity": sorted(r["best_contig_purity"] for r in present_rows
                                    if r["best_contig_purity"] is not None)[
                                 max(0, sum(1 for r in present_rows
                                            if r["best_contig_purity"] is not None) // 2)]
            if any(r["best_contig_purity"] is not None for r in present_rows) else None,
            "inter_group_mosaics": len(mosaics), "mosaic_detail": mosaics[:50],
            "false_joins_absent_arm": false_joins}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="recorded; one arm is one deterministic whole-library pass, so it "
                         "changes nothing in the assemble or score stages")
    ap.add_argument("--arms", default=",".join(ARMS),
                    help="comma-separated arm names from the ladder in this file")
    ap.add_argument("--stages", nargs="+", default=["assemble", "score"],
                    choices=["assemble", "score"])
    ap.add_argument("--reads", type=Path, required=True,
                    help="the library. Single-end on its own, or mate 1 with --r2")
    ap.add_argument("--r2", type=Path, default=None,
                    help="mate 2, ALREADY TRIMMED of its non-template 5' prefix where the "
                         "library has one. This script does not trim: Berlin's mate 2 "
                         "carries ~14 nt of non-template sequence that every experiment in "
                         "this repo removes, and a trimmer here would be a second one")
    ap.add_argument("--no-truth-panel", action="store_true",
                    help="skip the ground-truth scoring: for a library with no committed "
                         "panel, where the arm's only readouts are assembly shape and, "
                         "downstream, the residue it leaves")
    ap.add_argument("--run", default="SRR23421544", help="the run id as depth.tsv spells it")
    ap.add_argument("--candidates", type=Path,
                    default=ROOT / "results/holdout_candidates/candidates.fa")
    ap.add_argument("--depth", type=Path,
                    default=ROOT / "results/holdout_candidates/depth.tsv")
    ap.add_argument("--shortlist", type=Path,
                    default=ROOT / "results/holdout_candidates/shortlist.tsv")
    ap.add_argument("--megahit", type=Path, required=True)
    ap.add_argument("--minimap2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--memory", type=float, default=0.8, help="MEGAHIT --memory fraction")
    ap.add_argument("--min-identity", type=float, default=0.95,
                    help="a contig is said to carry a panel genome at this identity or "
                         "better; the genomes are the same molecules as the reads, so a "
                         "true recovery is ~1.0 and this is a floor, not a divergence band")
    ap.add_argument("--min-aln", type=int, default=200, help="shortest alignment counted")
    ap.add_argument("--min-mosaic", type=int, default=500,
                    help="span per lineage group before a contig counts as a mosaic")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    for a in arms:
        if a not in ARMS:
            raise SystemExit(f"unknown arm {a!r}; the ladder is {list(ARMS)}")

    t0 = time.perf_counter()
    asm = {}
    if "assemble" in args.stages:
        for a in arms:
            asm[a] = assemble(args, a)
            (args.out / f"assembly_{a}.json").write_text(json.dumps(asm[a], indent=1))
    else:
        for a in arms:
            asm[a] = json.loads((args.out / f"assembly_{a}.json").read_text())

    if "score" in args.stages and args.no_truth_panel:
        log("score stage skipped: --no-truth-panel")
    elif "score" in args.stages:
        panel, meta = build_panel(args)
        groups = lineage_groups(args, panel)
        (args.out / "panel.json").write_text(json.dumps(
            {"meta": meta, "groups": groups,
             "n_groups": len(set(groups.values()))}, indent=1))
        for a in arms:
            s = score(args, a, Path(asm[a]["contigs_fa"]), panel, meta, groups)
            (args.out / f"recovery_{a}.json").write_text(json.dumps(s, indent=1))
            log(f"[{a}] contigs {asm[a]['contigs']} n50 {asm[a]['n50']} | "
                f"recovered>=90% {s['recovered_90pct']}/{s['panel_present']} "
                f"median purity {s['median_purity']} mosaics {s['inter_group_mosaics']} "
                f"false joins {len(s['false_joins_absent_arm'])}")

    man = {"script": "experiments/path_settings.py", **provenance(__file__),
           "argv": sys.argv[1:], "arms": arms, "ladder": ARMS,
           "parameters": {k: str(v) for k, v in vars(args).items()},
           "seconds_total": time.perf_counter() - t0,
           "peak_rss_bytes": peak_rss_bytes(),
           "megahit_version": subprocess.run([str(args.megahit), "--version"],
                                             capture_output=True, text=True).stdout.strip(),
           "minimap2_version": subprocess.run([str(args.minimap2), "--version"],
                                              capture_output=True, text=True).stdout.strip(),
           "reads_sha256": sha256(args.reads) if args.reads.exists() else None}
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps(man, indent=1))
    log(f"manifest: {args.out / f'manifest_s{args.seed}.json'}")


if __name__ == "__main__":
    main()
