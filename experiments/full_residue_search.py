"""Search the WHOLE unplaced residue for named relatives, not a block of it.

Why this script exists. Every named-relative figure in this repo so far - the subject
table, the kind split, the leaderboard - rests on the 6,000 pairs the per-read report
samples, out of 1.37 million unplaced pairs. That is fine for a percentage and useless
for a leaderboard: a subject at the top of it is carried by two or three reads, and two
or three reads is a number with no shape. This runs the search over every unplaced pair
in the library so those rows have real counts behind them.

How it avoids a full cascade. A full cascade over 8.7 M pairs is tens of hours in this
container (measured: a warm 100,000-pair cascade does not finish inside 35 minutes), and
almost all of that work is spent re-deciding reads that the shared label table has
already decided. `library_label_classes.npz` holds one class per MATE for all 8,731,224
pairs in file order, and `unplaced` - both mates matching nothing - reproduces the
residue census's 1,372,021 pairs (15.71%) to the pair. So the pile is taken from there.

The bias that comes with that, with its direction. The label table was built from
nucleotide evidence only: no read-level Kraken2 and no translated search. Its `unplaced`
pile is therefore MORE GENEROUS than the cascade's residue - it contains reads the
cascade would have named - so a named-relative rate computed here is, if anything, a
small overestimate, and a "points at nothing" rate a small underestimate. The direction
is fixed and stated; it is not a two-sided uncertainty.

The control is run in the same pass, which is the only way it is worth anything. A
sample of the pairs the labeller DID place travels through the identical search at the
identical operating point, so "5.8% of the residue has a relative" can be read against
what the same instrument says about material that certainly has one. Its sample size is
its own and is reported rather than scaled.

Operating point. DIAMOND blastx at --sensitive, E-value reporting cutoff 100, best hit
per mate by bitscore, significance taken afterwards at 1e-5. --sensitive and not
--very-sensitive because the full pile at --very-sensitive does not fit the time this
container has; that is an operating-point change from the per-read report and every
number here is therefore reported as its own instrument, never merged with that one.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_facts
from subtract_known import provenance  # noqa: E402

#: Taken after the search, never during it: the chance band reads identically on residue
#: and already-identified reads, so a permissive reporting cutoff plus a strict
#: significance threshold is the only honest order.
SIGNIFICANT_EVALUE = 1e-5


def say(msg: str) -> None:
    """Progress on stderr. A run of this length is unreadable without it, and the
    per-index seconds are what price the next one."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def pair_class(class1: np.ndarray, class2: np.ndarray, codes: list[str]) -> np.ndarray:
    """Collapse the two mates to one class per pair, by the priority the census uses."""
    idx = {c: i for i, c in enumerate(codes)}
    out = np.full(class1.shape, idx["none"], dtype=np.int8)
    for name in ("unident", "other_ctg", "identified", "rrna"):
        if name not in idx:
            continue
        hit = (class1 == idx[name]) | (class2 == idx[name])
        out[hit] = idx[name]
    return out


def write_piles(r1: Path, r2: Path, wants: dict, trim_r2: int, paths: dict) -> dict:
    """Stream both FASTQs ONCE and write every pile in the same pass.

    Decompressing 8.7 M pairs twice is the dominant fixed cost of this script and it
    buys nothing: the piles are disjoint masks over one read order, so one pass can
    fill all of them. The pilot paid it per pile, which is why this is a function over
    a dict of masks rather than one called repeatedly.
    """
    tags = list(wants)
    n = {t: 0 for t in tags}
    handles = {t: open(paths[t], "w") for t in tags}
    try:
        with gzip.open(r1, "rt") as h1, gzip.open(r2, "rt") as h2:
            i = 0
            while True:
                a = h1.readline()
                if not a:
                    break
                s1 = h1.readline().rstrip("\n")
                h1.readline(); h1.readline()
                h2.readline()
                s2 = h2.readline().rstrip("\n")[trim_r2:]
                h2.readline(); h2.readline()
                for t in tags:
                    w = wants[t]
                    if i < w.size and w[i]:
                        name = a[1:].split()[0]
                        handles[t].write(f">{name}/1\n{s1}\n>{name}/2\n{s2}\n")
                        n[t] += 1
                i += 1
    finally:
        for h in handles.values():
            h.close()
    return n


def diamond_best(dbin: Path, db: Path, fa: Path, tsv: Path, threads: int,
                 mode: str, label: str) -> tuple[dict, float]:
    t0 = time.perf_counter()
    subprocess.run([str(dbin), "blastx", "--db", str(db), "--query", str(fa),
                    "--out", str(tsv), "--threads", str(threads), "--outfmt", "6",
                    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
                    "--max-target-seqs", "1", "--evalue", "100", mode, "--quiet"],
                   check=True)
    best: dict[str, dict] = {}
    with open(tsv) as h:
        for line in h:
            q, s, pid, length, ev, bits = line.rstrip("\n").split("\t")
            name = q.rsplit("/", 1)[0]
            bits = float(bits)
            if bits > best.get(name, {}).get("bitscore", -1.0):
                best[name] = {"subject": s, "pident": float(pid), "aln_len": int(length),
                              "evalue": float(ev), "bitscore": bits,
                              "mate": q.rsplit("/", 1)[1], "index": label}
    tsv.unlink(missing_ok=True)
    return best, time.perf_counter() - t0


def merge_best(a: dict, b: dict) -> dict:
    """Best hit per mate across the two indexes, ties to the viral index.

    Each hit records which index produced it at the point it was parsed. The pilot
    instead kept a side table keyed on id() of the hit dicts and then copied the dict
    before reading it, so the copy never matched either identity and a mate hit by
    only one index raised KeyError. A field on the row cannot come apart that way.
    """
    out = {}
    for name in set(a) | set(b):
        ha, hb = a.get(name), b.get(name)
        if hb is None:
            out[name] = dict(ha)
        elif ha is None:
            out[name] = dict(hb)
        else:
            out[name] = dict(ha if ha["bitscore"] >= hb["bitscore"] else hb)
    return out


def subjects_from(best: dict) -> dict:
    agg: dict[str, dict] = {}
    for name, h in best.items():
        if h["evalue"] >= SIGNIFICANT_EVALUE:
            continue
        a = agg.setdefault(h["subject"], {"subject": h["subject"], "index": h["index"],
                                          "pid": [], "aln": [], "ev": [], "reads": 0})
        a["reads"] += 1
        a["pid"].append(h["pident"])
        a["aln"].append(h["aln_len"])
        a["ev"].append(h["evalue"])
    rows = []
    for a in agg.values():
        rows.append({"subject": a["subject"], "index": a["index"],
                     "residue_reads": a["reads"],
                     "median_identity_percent": round(statistics.median(a["pid"]), 1),
                     "identity_range": [min(a["pid"]), max(a["pid"])],
                     "median_aligned_aa": int(statistics.median(a["aln"])),
                     "median_evalue": statistics.median(a["ev"]),
                     "worst_evalue": max(a["ev"])})
    rows.sort(key=lambda r: -r["residue_reads"])
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--classes", type=Path,
                    default=Path("/mnt/project-files/shared/resasm/library_label_classes.npz"))
    ap.add_argument("--r1", type=Path, required=True)
    ap.add_argument("--r2", type=Path, required=True)
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None)
    ap.add_argument("--dmnd", type=Path, nargs=2, required=True,
                    help="viral then sprot DIAMOND databases, in that order")
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--diamond-mode", default="--sensitive")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--control-pairs", type=int, default=200000,
                    help="pairs of already-placed material through the same search")
    ap.add_argument("--limit-residue", type=int, default=0,
                    help="0 means every unplaced pair; a small value is a pilot")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    # A trim length is a fact about a sequencing run, not a default. It is resolved
    # from experiments/libraries.json so a library this repo has no entry for is
    # refused by name rather than silently given Berlin's 14 nt, which on a library
    # with no non-template prefix deletes real bases from every mate, in the
    # direction of less recognition.
    library_facts.resolve_args(args, ["trim_r2"])

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    t_start = time.perf_counter()

    d = np.load(args.classes)
    codes = [str(c) for c in d["class_codes"]]
    pc = pair_class(d["class1"], d["class2"], codes)
    none_i = codes.index("none")
    unplaced = pc == none_i
    placed = ~unplaced
    n_lib = int(pc.size)

    want_res = unplaced.copy()
    if args.limit_residue:
        # A pilot draws at random across the whole library, never a prefix: the label
        # table is what showed this file's order to be spatially structured.
        idx = np.flatnonzero(unplaced)
        keep = rng.choice(idx, size=min(args.limit_residue, idx.size), replace=False)
        want_res = np.zeros_like(unplaced)
        want_res[keep] = True
    idx_p = np.flatnonzero(placed)
    keep_p = rng.choice(idx_p, size=min(args.control_pairs, idx_p.size), replace=False)
    want_ctl = np.zeros_like(placed)
    want_ctl[keep_p] = True

    wants = {"residue": want_res, "control": want_ctl}
    fas = {t: args.work / f"{t}_s{args.seed}.fa" for t in wants}
    t0 = time.perf_counter()
    npair = write_piles(args.r1, args.r2, wants, args.trim_r2, fas)
    extract_s = time.perf_counter() - t0
    say(f"extracted {npair} pairs in {extract_s:.0f}s")

    piles = {}
    for tag in wants:
        fa, n = fas[tag], npair[tag]
        secs = {}
        per_index = {}
        for db, label in zip(args.dmnd, ("viral", "sprot")):
            b, sec = diamond_best(args.diamond_bin, db, fa,
                                  args.work / f"{tag}_{label}_s{args.seed}.tsv",
                                  args.threads, args.diamond_mode, label)
            per_index[label] = b
            secs[label] = round(sec, 1)
            say(f"{tag}/{label}: {len(b)} mates hit, {sec:.0f}s over {2 * n} mates")
        best = merge_best(per_index["viral"], per_index["sprot"])
        sig = {k: v for k, v in best.items() if v["evalue"] < SIGNIFICANT_EVALUE}
        piles[tag] = {
            "pairs": n, "mates": 2 * n,
            "searched_mates_with_any_hit": len(best),
            "mates_with_significant_hit": len(sig),
            "subjects": subjects_from(best) if tag == "residue" else None,
            "subject_count": len({v["subject"] for v in sig.values()}),
            # Counts, not a name list. A residue subject that the already-placed pile
            # also hits is ALREADY SEEN, and the page says so on the row - so what is
            # needed here is how often, which a set of names cannot answer. It is a
            # sample of the placed pile, so it bounds from below: a subject absent here
            # may still be placed material that this sample missed.
            "control_subject_reads": (collections.Counter(
                v["subject"] for v in sig.values()) if tag == "control" else None),
            "extract_seconds_all_piles": round(extract_s, 1),
            "diamond_seconds": secs,
        }
        fa.unlink(missing_ok=True)

    out = {
        "unit": "one row per subject; counts are MATES, since a hit is per mate",
        "library_pairs": n_lib,
        "unplaced_pairs": int(unplaced.sum()),
        "unplaced_fraction": round(float(unplaced.mean()), 5),
        "operating_point": {"diamond_mode": args.diamond_mode,
                            "reporting_evalue": 100,
                            "significance_evalue": SIGNIFICANT_EVALUE,
                            "max_target_seqs": 1,
                            "trim_r2": args.trim_r2},
        "piles": piles,
        "bias": ("the pile comes from a nucleotide-only label table, so it is MORE "
                 "generous than the cascade's residue and the named-relative rate here "
                 "is if anything a small overestimate. One direction, not two."),
        "not_comparable_with": ("the per-read report's subject table, which ran "
                                "--very-sensitive on 6,000 pairs. Different instrument; "
                                "do not merge the two."),
        "wall_seconds": round(time.perf_counter() - t_start, 1),
        "seed": args.seed,
        **provenance(__file__),
    }
    dest = args.out / f"full_residue_search_s{args.seed}.json"
    dest.write_text(json.dumps(out, indent=1))
    r = piles["residue"]
    c = piles["control"]
    print(f"residue {r['pairs']:,} pairs -> {r['mates_with_significant_hit']:,} mates "
          f"with a significant relative across {r['subject_count']:,} subjects")
    print(f"control {c['pairs']:,} pairs -> {c['mates_with_significant_hit']:,} mates")
    print(f"{out['wall_seconds']}s -> {dest}")


if __name__ == "__main__":
    main()
