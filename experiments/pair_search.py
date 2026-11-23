#!/usr/bin/env python3
"""The non-overlapping mate pair as the unit a homology search scores.

WHAT IS NEW HERE, AND WHAT IS ALREADY DONE
------------------------------------------
Three things this repo has already measured are NOT this:

- **Merging** overlapping mates into one longer read (`results/profile_reach/`). That
  needs the mates to overlap, it covers 35.1% of Berlin's pairs, and it works: at 20%
  divergence recognition goes 10.5% -> 38.4%. The pair there becomes a longer SINGLE
  query.
- **The pair as evidence** (`results/crossplant/`, `results/residue_structure/`): both
  mates unwitnessed against one mate unwitnessed. The search still scores one read at a
  time; the pair enters afterwards, as a conjunction over per-read verdicts.
- **Mate width for novelty** (`results/mate_width/`): what a wider insert buys a
  per-read FLAG, with the rho_fp model for the conjunction's merit.

What no arm covers is the one most search tools advertise: a **paired mode**, where two
mates that do not overlap are handed to the searcher as ONE unit and the subject's
significance is computed over both. That is what this script measures, on the one
sequence whose true distance to the index is known and not chosen by us.

THE FOUR MODES, AND WHY THESE FOUR
----------------------------------
All four are computed from the same reads and the same two indexes, so the comparison is
mode against mode on one instrument.

1. `best_mate` - the production unit. The pair counts as identified when EITHER mate's
   own best HSP clears the threshold. This is the baseline every other mode has to beat,
   and it is `results/holdout_run/`'s per-read readout lifted to the pair.
2. `concat` - the two mates concatenated into one query, which is what a tool means by
   "paired input" when its scoring is k-mer based (Kraken2 joins mates with a spacer).
   For an alignment search the prediction is that this LOSES: the E-value of an HSP is
   proportional to the query length, the HSP itself is unchanged because no alignment
   spans the junction, so doubling the query doubles E and the pair clears less than the
   better mate did alone. The junction census below is what turns that from a prediction
   into a measurement: every concatenated HSP is classified as lying inside mate 1,
   inside mate 2, or crossing the junction. If nothing crosses, the searcher is scoring
   the mates independently inside one query while charging both of them for the length.
3. `concat_n` - the same with ten N's between the mates, the standard mitigation, which
   stops a spurious alignment running across the junction and (for blastx) stops a frame
   running through it. Reported separately because it is a different query length.
4. `fisher` - the honest joint statistic. Per subject, take each mate's best HSP against
   THAT subject, convert the two E-values to p-values, combine with Fisher's method, and
   charge the number of subjects the pair tried as a multiple-testing correction. This is
   the pair-level analogue of BLAST's own sum statistics over multiple HSPs of one query,
   and it is the only mode here that can clear something neither mate cleared alone. It
   is implemented in this file rather than taken from a tool's `--paired` switch
   precisely so that nothing can silently score the mates independently and call itself
   paired.
5. `agree` - both mates' best hit names the same subject, with each mate held to a
   RELAXED threshold. The consistency rule a mapper's paired rescue applies. Its
   threshold is set from the null, not chosen.

WHY THE READS ARE SIMULATED, AND WHAT THAT COSTS
------------------------------------------------
The holdout's containing library is single-end in practice (filed PAIRED, no second mate
present at ENA for any of the study's three runs), so the pair has to be constructed. The
genomes are real, the distance to the index is real and was not chosen by us, the insert
distribution and per-mate error rates are Berlin's measured ones, and the read length is
the holdout library's own 145 nt so that the `best_mate` rate can be read against
`results/holdout_run/`'s 5.25-7.00%. What is simulated is the pairing and the errors.
That is the weakest part of this design and it is stated here rather than in a caveat at
the bottom: a simulated pair cannot carry a chimaeric library molecule, so the `agree`
rule's false-positive rate on real material is NOT measured by the arms below. The
shuffled arm bounds chance alignment; it does not bound library artefact.

NON-OVERLAPPING BY CONSTRUCTION. An insert shorter than twice the read length makes the
mates overlap, and then the right instrument is the merge, which is already measured. So
inserts are drawn from Berlin's measured distribution TRUNCATED to the non-overlapping
band. Keeping that band costs real molecules - `results/mate_width/` measures insert >=300
at 15.5% of Berlin's pairs - and the cost is part of the answer, not an aside. Two bands
are run: `near`, just past the overlap limit, and `wide`, where mate separation is large
enough that `results/mate_width/` measured rho falling.

EVERY MODE IS READ AT TWO OPERATING POINTS
------------------------------------------
A nominal E-value is not a false-positive rate, and a joint statistic with a
multiple-testing correction invented in this file is exactly the kind of number that can
look good by being read at a threshold nobody calibrated. So each mode is reported at
E < 1e-5 (this repo's operating point, for comparability) AND at its own
matched-null point: the strictest threshold at which the mode's shuffled-arm rate is zero
over all seeds pooled. A mode that gains sensitivity only by moving its threshold has
gained nothing, and the matched-null column is where that shows.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import random
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from provenance import provenance
from search_hits import DEFAULT_COLUMNS, outfmt_args, read_hits
from shuffle_control import shuffle_arm
from berlin_io import revcomp
from spike_holdout import simulate_pairs

ROOT = Path(__file__).resolve().parents[1]

# Berlin's measured insert distribution (PR #10) and its measured per-mate error rates
# (PR #11). Named here rather than inline so the truncation below is visibly a truncation
# OF a measured distribution and not a new distribution.
INSERT_QUANTILES = [(0.50, 194), (0.75, 246), (0.95, 363), (0.99, 485), (1.0, 600)]
ERR1, ERR2 = 0.00089, 0.00117

# The arms. Distances are the finished-genome figures from
# results/holdout_candidates/shortlist.tsv, measured against viral RefSeq as it stands
# today, which makes every one of them a FLOOR (a later deposit can only have brought the
# nearest relative closer). All four genomes were first deposited 2023-07-31, after the
# 2022-05-10 cutoff, so all four are removed by the same date holdout the indexes carry.
# `kind` is the thing to read first. `real_mate1` arms take a REAL library read as mate 1,
# at its real position, and construct only its partner, so the arm inherits the library's
# own coverage and the single-mate rate reproduces the published one. `simulated` arms
# draw both mates uniformly off the genome; they are kept because the GAP between the two
# kinds is itself the measurement that uniform sampling is an easier query than the
# library gives (47% against 6.1% per mate at nucleotide level, first pass here).
ARMS = {
    "candidate": dict(acc="OQ722337.1", nt_div=30.28, prot_id=66.452, kind="real_mate1",
                      role="the blind: the genome results/holdout_run/ measured at read "
                           "level, 30.3% nucleotide divergent, with its real reads"),
    "leak": dict(acc="OQ722330.1", nt_div=0.08, prot_id=99.813, kind="real_mate1",
                 role="known material: real reads of a sequence 0.08% divergent from a "
                      "pre-cutoff relative. The comparison arm every rate is read "
                      "against"),
    "candidate_sim": dict(acc="OQ722337.1", nt_div=30.28, prot_id=66.452,
                          kind="simulated",
                          role="the same genome sampled uniformly, which is what a "
                               "simulated-pair experiment would have measured"),
    # The Berlin arms, which are the primary block: REAL non-overlapping pairs of the
    # discriminating library, with no construction anywhere. They carry no known distance
    # - nothing in that water has a true answer attached - so they answer the mode
    # question (does a joint rule clear what neither mate clears) rather than the
    # distance question, and they answer it on the material the production tool runs on.
    # The class of each mate comes from the committed whole-library label table, so the
    # piles are the ones results/residue_census/ defined and not a new split.
    "berlin_identified": dict(kind="berlin", classes=("identified", "identified"),
                              role="known material: both mates match this sample's own "
                                   "identified contigs. The comparison arm, and where "
                                   "every mode's direction is measured"),
    "berlin_residue": dict(kind="berlin", classes=("none", "none"),
                           role="the residue: neither mate matches SILVA or any own "
                                "contig. The pile the production tool is trying to "
                                "shrink"),
    "rung17_sim": dict(acc="OQ722335.1", nt_div=17.09, prot_id=84.459, kind="simulated",
                       role="a nearer blind, uniform sampling, where single-mate "
                            "recognition is partly alive and a joint rule has room to "
                            "show a gain"),
}
# Inserts, in nucleotides. `near` starts one base past the overlap limit for 145 nt reads.
# `near` starts at 301 rather than 291 so that a pair stays non-overlapping for every
# read length this library carries (max 150), without the simulator having to nudge the
# insert upward per pair - a nudge would make the realised band depend on mate length.
BANDS = {"near": (301, 430), "wide": (600, 1000)}

# The levels, as (program, database stem). `prot_dmnd` is the protein level the five-seed
# run uses and `prot` is BLASTX, kept as a cross-check on one seed rather than dropped.
# WHY: BLASTX of 145 nt queries against 703,330 proteins runs at ~0.7 s per query here,
# which is 14 hours for the five-seed grid, while DIAMOND --ultra-sensitive is minutes.
# results/reach/ measured eleven search methods on one panel and found DIAMOND at three
# sensitivities, BLASTP and full Smith-Waterman returning the SAME nearest relative on all
# 16 contigs - the constraint is the significance statistics, not the seeding - so the
# substitution is one this repo has already calibrated. It is still a different instrument
# with a different E-value scale, so the two are never mixed inside one comparison: every
# mode comparison happens within a level.
LEVELS = {"nt": ("blastn", "db_holdout_nt"),
          "prot_dmnd": ("diamond_blastx", "db_holdout_prot_dmnd")}
BLASTX_LEVEL = {"prot": ("blastx", "db_holdout_prot")}
# Swiss-Prot, for the Berlin block only. The viral protein release answers "is there a
# viral relative", which is the right question for the holdout arms and the wrong
# denominator for a wastewater residue pile that is mostly not viral: a mode that
# converts nothing against 703k viral proteins might convert something against a
# reviewed index covering all life. It is the same `sprot` index the cascade already
# uses, so the operating point is one this repo has elsewhere.
SPROT_LEVEL = {"sprot_dmnd": ("diamond_blastx", "db_sprot_dmnd")}
SPACER = "N" * 10


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load_genomes(path: Path) -> dict[str, str]:
    """Accession -> sequence from the committed candidate FASTA.

    Raises on a missing accession rather than skipping an arm: an arm that silently does
    not run leaves a directory that looks complete, which is the bug family this repo
    keeps finding (a default where a refusal belongs).
    """
    if not path.exists():
        raise SystemExit(f"{path} does not exist; this is the committed candidate FASTA "
                         f"and nothing here can be derived without it")
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def _cdf_at(x: float) -> float:
    """The measured insert distribution's CDF at x, piecewise-linear between its quantiles.

    Five quantiles are all that was measured, so anything between them is interpolation
    and is treated as such: the alternative, snapping an insert to the nearest measured
    quantile, puts every simulated molecule at one of five lengths and makes a band of
    the distribution into a spike. The first pass of this experiment did exactly that and
    every pair in the `near` band came out at 420 nt.
    """
    knots = [(0.0, 2 * 50)] + [(q, v) for q, v in INSERT_QUANTILES]
    for i in range(1, len(knots)):
        q0, v0 = knots[i - 1]
        q1, v1 = knots[i]
        if x <= v1:
            if v1 == v0:
                return q1
            return q0 + (q1 - q0) * (x - v0) / (v1 - v0)
    return 1.0


def conditional_ladder(lo: int, hi: int, steps: int = 40) -> tuple[list, dict]:
    """The measured distribution conditioned on [lo, hi], as a ladder simulate_pairs takes.

    Returns the ladder and a note saying which of two things it is, because the two are
    not the same claim:

    - `conditioned`: the band carries real mass, so the ladder is the measured shape
      renormalised onto it. The band's share of real molecules is reported with it, since
      requiring a non-overlapping pair THROWS AWAY that much of the library and the
      discard is part of any paired mode's cost.
    - `counterfactual_uniform`: the band lies beyond the measured maximum, so there is no
      shape to condition and the ladder is uniform over the band. That arm then answers
      "what would a wider-insert protocol give", not "what does this library give", and
      it must be read as the protocol question it is.
    """
    c_lo, c_hi = _cdf_at(lo), _cdf_at(hi)
    mass = c_hi - c_lo
    if mass <= 0.005:
        ladder = [((i + 1) / steps, int(lo + (hi - lo) * (i + 1) / steps))
                  for i in range(steps)]
        return ladder, {"mode": "counterfactual_uniform", "band_mass_of_library": 0.0,
                        "note": "band is beyond the measured insert maximum (600 nt)"}
    ladder = []
    for i in range(steps):
        u = (i + 1) / steps
        target = c_lo + mass * u
        # invert the piecewise-linear CDF by bisection on length; 1 nt resolution
        a, b = lo, hi
        while b - a > 1:
            m = (a + b) // 2
            if _cdf_at(m) < target:
                a = m
            else:
                b = m
        ladder.append((u, b))
    return ladder, {"mode": "conditioned", "band_mass_of_library": round(mass, 4)}


def map_library(fastq: Path, candfa: Path, wanted: dict[str, str], cache: Path) -> dict:
    """Real library reads whose primary alignment is one of the arm genomes, with where.

    Selection is `results/holdout_run/`'s, unchanged, so the `candidate` pool here is the
    same pool whose single-mate rate that directory published: a read's arm is the
    deposited sequence its primary minimap2 `sr` alignment goes to, and a read can
    therefore be in only one arm.

    WHY THE REAL READ IS NEEDED AT ALL, which is the thing the first pass of this
    experiment got wrong. Pairs simulated off the deposited genome by uniform fragment
    sampling identify at ~47% per mate against the held-out nucleotide index, where the
    same genome's REAL reads identify at 6.1%. The genome's recognisable stretches are a
    small part of it and the library's coverage is not uniform over them, so a uniformly
    sampled read is a far easier query than a real one. Measuring a paired mode on
    uniform reads would therefore be measuring it at a distance the library does not
    have. Mate 1 is consequently a real read, at its real position, and only its partner
    is constructed.
    """
    if not fastq.exists():
        raise SystemExit(f"{fastq} does not exist: the real-mate arms cannot be built "
                         f"without the library, and a simulated-only run is a different "
                         f"experiment (pass --arms with the _sim arms only to run it)")
    if cache.exists():
        return json.loads(cache.read_text())
    import mappy
    idx = mappy.Aligner(str(candfa), preset="sr", n_threads=4)
    if not idx:
        raise SystemExit(f"mappy could not index {candfa}")
    pools = {acc: [] for acc in wanted.values()}
    n = 0
    for name, seq, _q in mappy.fastx_read(str(fastq)):
        n += 1
        best = None
        for h in idx.map(seq):
            if h.is_primary and (best is None or h.mlen > best.mlen):
                best = h
        if best is not None and best.ctg in pools:
            pools[best.ctg].append([name, seq, best.r_st, best.r_en, best.strand])
    out = {"reads_scanned": n, "pools": pools,
           "pool_sizes": {k: len(v) for k, v in pools.items()}}
    cache.write_text(json.dumps(out))
    log(f"  mapped {n} library reads; " +
        ", ".join(f"{k}={len(v)}" for k, v in pools.items()))
    return out


BERLIN_BLOCK = 1_500_000      # the stride results/subtract_known/ uses for this library
BERLIN_TRIM_R2 = 14           # the non-template 5' prefix on mate 2, as every arm trims


def berlin_pairs(r1: Path, r2: Path, classes_npz: Path, want: tuple[str, str],
                 seed: int, n_pairs: int, block: int = BERLIN_BLOCK) -> tuple[list, dict]:
    """Real non-overlapping pairs of one labelled pile, from one block of the library.

    Four things here are the repo's existing conventions rather than choices made for this
    experiment, and each one has a measured reason:

    - **The pile comes from the committed label table**, keyed by spot, so `identified`
      and `none` mean what `results/residue_census/` made them mean.
    - **Five seeds are five disjoint contiguous blocks** of one flow cell at a 1,500,000
      stride, not five random subsamples. Blocks are how this library has always been
      sampled here, and a statistic estimated from the blocks' own spread reverses sign
      against the full series (`results/block_position_bias/`), so the spread below is
      reported as an observed range and nothing is estimated from it.
    - **Mate 2 loses its first 14 nt**, the non-template prefix. Without the trim every
      mate 2 donates a constant 14 nt of unmatched sequence.
    - **Overlapping pairs are excluded with the committed merger**, `merge_real_pair` from
      `experiments/profile_reach.py`, so "non-overlapping" is the complement of the
      criterion that directory measured rather than a new test. A pair it cannot merge is
      either genuinely non-overlapping or overlapping with more mismatches than the
      criterion allows, so this is a lower bound on separation, not a guarantee.
    """
    import numpy as np
    from profile_reach import merge_real_pair
    z = np.load(classes_npz)
    codes = list(z["class_codes"])
    c1, c2 = z["class1"], z["class2"]
    w1, w2 = codes.index(want[0]), codes.index(want[1])
    start = (seed - 1) * block
    out, census = [], dict(block_start=start, scanned=0, in_pile=0, overlapping=0,
                           too_short=0)
    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        i = -1
        while len(out) < n_pairs:
            h1 = f1.readline()
            s1 = f1.readline().strip()
            f1.readline(); f1.readline()
            h2 = f2.readline()
            s2 = f2.readline().strip()
            f2.readline(); f2.readline()
            if not h1 or not h2:
                break
            i += 1
            if i < start:
                continue
            if i >= len(c1):
                break
            census["scanned"] += 1
            if c1[i] != w1 or c2[i] != w2:
                continue
            census["in_pile"] += 1
            m1 = s1.upper()
            m2 = s2.upper()[BERLIN_TRIM_R2:]
            if len(m1) < 50 or len(m2) < 50:
                census["too_short"] += 1
                continue
            merged, _ov = merge_real_pair(m1, m2)
            if merged is not None:
                census["overlapping"] += 1
                continue
            out.append((f"spot{i + 1}", m1, m2))
    census["pairs"] = len(out)
    census["nonoverlap_share_of_pile"] = (
        round(len(out) / census["in_pile"], 4) if census["in_pile"] else None)
    return out, census


def partner_for(read: list, genome: str, rng: random.Random, quant, read_len: int,
                err: float):
    """The other mate of a real read: same molecule, the far end of a real-length insert.

    The real read is mate 1 at its real position and real error profile. The insert is
    drawn from Berlin's measured distribution conditioned on the non-overlapping band,
    laid downstream of the real read in the direction the read's own alignment implies,
    and the partner is the far end of that fragment reverse-complemented, with the
    measured R2 substitution rate applied.

    DIRECTION OF THE BIAS, since there is no control that removes it: the partner carries
    no library artefact, no adapter, no quality collapse and no within-population
    variation beyond the substitution rate, so it is a CLEANER mate than the library would
    have given. Any gain a joint rule shows here is therefore an upper bound on the gain
    it would show on real pairs, which is the safe direction for a negative result and
    must be stated beside a positive one.
    """
    name, seq, r_st, r_en, strand = read
    L = len(seq)
    u = rng.random()
    ins = quant[-1][1]
    for q, v in quant:
        if u <= q:
            ins = v
            break
    ins = max(ins, 2 * L + 1)
    if strand > 0:
        f_st, f_en = r_st, r_st + ins
    else:
        f_en, f_st = r_en, r_en - ins
    if f_st < 0 or f_en > len(genome):
        return None, None          # fragment runs off the genome; counted by the caller
    frag = genome[f_st:f_en]
    far = revcomp(frag)[:L] if strand > 0 else frag[:L]
    alt = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
    far = "".join(rng.choice(alt[b]) if b in alt and rng.random() < err else b
                  for b in far)
    return far, dict(insert=ins, frag_start=f_st, frag_end=f_en, strand=strand,
                     mate1_len=L, separation=ins - 2 * L)


def stage_reads(out: Path, genomes: dict[str, str], seeds, band: str, n_pairs: int,
                read_len: int, arms: list, fastq: Path | None, candfa: Path,
                berlin: tuple | None = None, classes_npz: Path | None = None) -> dict:
    """One FASTA per arm per seed per query form, plus the shuffled arm and its census."""
    lo, hi = BANDS[band]
    if lo <= 2 * read_len:
        raise SystemExit(f"band {band} starts at {lo}, which is not past the overlap "
                         f"limit of {2 * read_len} for {read_len} nt reads")
    quant, qnote = conditional_ladder(lo, hi)
    census = {"band": band, "insert_band": [lo, hi], "insert_ladder": quant,
              "insert_ladder_note": qnote, "read_len": read_len,
              "pairs_per_arm": n_pairs, "arms": {}}
    log(f"  band {band} = {lo}-{hi} nt: {qnote}")
    need_real = [a for a in arms if ARMS.get(a, {}).get("kind") == "real_mate1"]
    lib = None
    if need_real:
        lib = map_library(fastq, candfa,
                          {a: ARMS[a]["acc"] for a in need_real},
                          out / "library_pools.json")
        census["library"] = {"reads_scanned": lib["reads_scanned"],
                             "pool_sizes": lib["pool_sizes"]}
    geom_rows = []
    for seed in seeds:
        shuffled_src = {}
        for arm in arms:
            if arm == "shuffled":
                continue
            spec = ARMS[arm]
            g = None
            if spec["kind"] != "berlin":
                g = genomes.get(spec["acc"])
                if g is None:
                    raise SystemExit(f"arm {arm}: {spec['acc']} absent from the "
                                     f"candidate FASTA")
            rng = random.Random(f"pair_search/{band}/{arm}/{seed}")
            if spec["kind"] == "berlin":
                if not (berlin and berlin[0].exists() and berlin[1].exists()
                        and classes_npz and classes_npz.exists()):
                    raise SystemExit("the Berlin arms need --berlin-r1/--berlin-r2 and "
                                     "--label-classes; refusing rather than running a "
                                     "arm that is not what its name says")
                pairs, cen = berlin_pairs(berlin[0], berlin[1], classes_npz,
                                          spec["classes"], seed, n_pairs)
                if len(pairs) < n_pairs:
                    raise SystemExit(f"arm {arm} seed {seed}: only {len(pairs)} "
                                     f"non-overlapping pairs of that pile in the block; "
                                     f"a short arm is a different experiment")
                census["arms"].setdefault(arm, {})[f"s{seed}"] = cen
            elif spec["kind"] == "simulated":
                pairs = simulate_pairs(g, n_pairs, rng, read_len, quant, ERR1, ERR2)
                census["arms"].setdefault(arm, {})[f"s{seed}"] = dict(pairs=len(pairs))
            else:
                pool = lib["pools"][spec["acc"]]
                if len(pool) < n_pairs:
                    raise SystemExit(f"arm {arm}: pool has {len(pool)} real reads, "
                                     f"fewer than the {n_pairs} asked for")
                # Draw until the arm is full rather than dropping to whatever survives.
                # A read whose fragment would run off the genome end is a read near the
                # end of a short genome, and simply losing those would make the arm
                # smaller AND biased towards the middle of the genome, with nothing in
                # the output saying by how much. The draw order is still the seeded
                # permutation, so the arm is reproducible, and the discard is counted.
                order = rng.sample(range(len(pool)), len(pool))
                pairs, off_end, used = [], 0, 0
                for idx in order:
                    if len(pairs) >= n_pairs:
                        break
                    used += 1
                    rd = pool[idx]
                    far, geo = partner_for(rd, g, rng, quant, read_len, ERR2)
                    if far is None:
                        off_end += 1
                        continue
                    nm = f"{arm}{len(pairs)}"
                    pairs.append((nm, rd[1], far))
                    geom_rows.append(dict(arm=arm, seed=seed, pair=nm, **geo))
                lens = sorted(len(p[1]) for p in pairs)
                census["arms"].setdefault(arm, {})[f"s{seed}"] = dict(
                    pairs=len(pairs), real_reads_offered=used,
                    mate1_len_median=(lens[len(lens) // 2] if lens else None),
                    mate1_len_min=(lens[0] if lens else None),
                    mate1_len_max=(lens[-1] if lens else None))
            write_forms(out, arm, seed, band, pairs)
            if arm in ("candidate", "leak", "candidate_sim", "berlin_identified",
                       "berlin_residue"):
                for nm, r1, r2 in pairs:
                    shuffled_src[f"{arm}:{nm}/1"] = r1
                    shuffled_src[f"{arm}:{nm}/2"] = r2
        # The null shuffles the SAME reads the test arms use, mate by mate, so a
        # difference between test and null cannot come from a different draw. Both mates
        # of a pair are shuffled independently, which is the right null for a joint rule:
        # it destroys the sequence while keeping the pairing, so any gain the null shows
        # is the joint statistic's own false-positive rate rather than lost linkage.
        rng = random.Random(f"pair_search/{band}/shuffled/{seed}")
        res = shuffle_arm(shuffled_src, rng)
        sh, cen = res["shuffled"], res["census"]
        pairs = []
        for key in sorted({k.rsplit("/", 1)[0] for k in sh}):
            r1, r2 = sh.get(f"{key}/1"), sh.get(f"{key}/2")
            if r1 and r2:                      # a pair needs both mates shuffled
                pairs.append((key.replace(":", "_"), r1, r2))
        write_forms(out, "shuffled", seed, band, pairs)
        cen["pairs_with_both_mates"] = len(pairs)
        census["arms"].setdefault("shuffled", {})[f"s{seed}"] = cen
        log(f"  {band} s{seed}: shuffled arm {len(pairs)} pairs "
            f"(offered {cen['offered']}, dropped identical {cen['dropped_identical']})")
    if geom_rows:
        with open(out / f"insert_geometry_{band}.tsv", "w") as w:
            keys = list(geom_rows[0])
            w.write("\t".join(keys) + "\n")
            for r in geom_rows:
                w.write("\t".join(str(r[k]) for k in keys) + "\n")
    (out / f"reads_census_{band}.json").write_text(json.dumps(census, indent=1))
    return census


def write_forms(out: Path, arm: str, seed: int, band: str, pairs) -> None:
    """The four query forms of one arm: each mate alone, and the two concatenations."""
    d = out / "queries"
    d.mkdir(parents=True, exist_ok=True)
    with open(d / f"{band}_{arm}_s{seed}.fa", "w") as w:
        for nm, r1, r2 in pairs:
            w.write(f">{nm}/1\n{r1}\n>{nm}/2\n{r2}\n")
    with open(d / f"{band}_{arm}_s{seed}_cat.fa", "w") as w:
        for nm, r1, r2 in pairs:
            w.write(f">{nm}/cat\n{r1}{r2}\n")
    with open(d / f"{band}_{arm}_s{seed}_catn.fa", "w") as w:
        for nm, r1, r2 in pairs:
            w.write(f">{nm}/catn\n{r1}{SPACER}{r2}\n")
    # the junction offsets are needed to classify a concatenated HSP, and they are a
    # property of the file, so they are written beside it rather than recomputed later
    with open(d / f"{band}_{arm}_s{seed}_geom.tsv", "w") as w:
        w.write("name\tlen1\tlen2\tcat_junction\tcatn_junction\tcatn_spacer\n")
        for nm, r1, r2 in pairs:
            w.write(f"{nm}\t{len(r1)}\t{len(r2)}\t{len(r1)}\t{len(r1)}\t{len(SPACER)}\n")


def stage_positional(out: Path, blastbin: Path, dbdir: Path, genomes: dict, acc: str,
                     read_len: int, step: int = 10) -> dict:
    """Where on the genome a read is recognisable, against where the library's reads are.

    This exists because of a 20-fold discrepancy this experiment ran into and had to
    explain: pairs simulated by uniform sampling off the deposited genome identify at
    65-68% per pair at nucleotide level, while REAL reads of the same genome identify at
    2.4-3.5%. Two candidate causes - the genome's recognisable sequence being concentrated,
    and the library's coverage not being uniform over it - are separable, and this stage
    separates them: tile the genome into read-length windows and search every window, so
    recognisability is a function of position, then weigh those positions by the real
    reads' own mapped starts.

    The weighted figure is what a real library would give under this positional
    distribution, and the gap between it and the flat average is the simulation's error.
    """
    g = genomes[acc]
    q = out / "positional_windows.fa"
    with open(q, "w") as w:
        for i in range(0, max(len(g) - read_len, 1), step):
            w.write(f">pos{i}\n{g[i:i + read_len]}\n")
    rows = {}
    for level, (prog, db) in (("nt", ("blastn", "db_holdout_nt")),
                              ("prot_dmnd", ("diamond_blastx", "db_holdout_prot_dmnd"))):
        o = out / f"positional_{level}.tsv"
        if not o.exists():
            if prog == "blastn":
                cmd = [str(blastbin / "blastn"), "-task", "blastn", "-query", str(q),
                       "-db", str(dbdir / db), "-evalue", "10", "-max_target_seqs", "50",
                       "-num_threads", "4", "-out", str(o)]
                cmd += outfmt_args(DEFAULT_COLUMNS, diamond=False)
            else:
                cmd = [str(blastbin / "diamond"), "blastx", "-q", str(q), "-d",
                       str(dbdir / db), "-e", "10", "--max-target-seqs", "50", "-p", "4",
                       "-o", str(o), "--very-sensitive", "--quiet"]
                cmd += outfmt_args(DEFAULT_COLUMNS, diamond=True)
            subprocess.run(cmd, check=True, capture_output=True)
        best = {}
        for r in read_hits(o, program=prog):
            pos = int(r["qseqid"][3:])
            if pos not in best or r["evalue"] < best[pos]:
                best[pos] = r["evalue"]
        rows[level] = best
    pools = json.loads((out / "library_pools.json").read_text())["pools"].get(acc, [])
    starts = [p[2] for p in pools]
    nwin = len(range(0, max(len(g) - read_len, 1), step))
    res = {"accession": acc, "genome_len": len(g), "windows": nwin, "step": step,
           "read_len": read_len, "real_reads": len(starts), "levels": {}}
    with open(out / "positional_profile.tsv", "w") as w:
        w.write("position\tnt_evalue\tprot_evalue\treal_reads_starting_here\n")
        cov = Counter((s // step) * step for s in starts)
        for i in range(0, max(len(g) - read_len, 1), step):
            w.write(f"{i}\t{rows['nt'].get(i, '')}\t{rows['prot_dmnd'].get(i, '')}"
                    f"\t{cov.get(i, 0)}\n")
    for level, best in rows.items():
        sig = {i for i, e in best.items() if e < 1e-5}
        flat = len(sig) / nwin if nwin else 0.0
        weighted = (sum(1 for s in starts if ((s // step) * step) in sig) / len(starts)
                    if starts else 0.0)
        res["levels"][level] = {
            "windows_significant_pct": round(100 * flat, 3),
            "coverage_weighted_pct": round(100 * weighted, 3),
            "note": ("flat is what uniform sampling sees; weighted is the same windows "
                     "weighted by where the library's reads actually start")}
    (out / "positional.json").write_text(json.dumps(res, indent=1))
    log(f"  positional {acc}: " + "; ".join(
        f"{k} flat {v['windows_significant_pct']}% weighted {v['coverage_weighted_pct']}%"
        for k, v in res["levels"].items()))
    return res


def stage_search(out: Path, blastbin: Path, dbdir: Path, seeds, band: str, arms,
                 threads: int, evalue: str, levels: dict,
                 sensitivity: str = "very-sensitive") -> dict:
    """Search every query form of every arm at a PERMISSIVE floor.

    The floor is deliberately loose (`--evalue 10` by default) because every mode here is
    a threshold applied afterwards, and a mode that combines two sub-threshold mates
    cannot be computed from output that has already thrown the sub-threshold rows away.
    The thresholds live in the report stage, where they are visible, and the search is run
    once.
    """
    timings = {}
    for level, (prog, db) in levels.items():
        dbp = dbdir / db
        diamond = prog.startswith("diamond")
        exists = (dbp.with_suffix(".dmnd").exists() if diamond else
                  dbp.with_suffix(".ndb").exists() or dbp.with_suffix(".pdb").exists())
        if not exists:
            raise SystemExit(f"{dbp} is not a {'DIAMOND' if diamond else 'BLAST'} "
                             f"database: that is 'not searched', not 'zero hits'")
        for seed in seeds:
            for arm in arms:
                for form in ("", "_cat", "_catn"):
                    q = out / "queries" / f"{band}_{arm}_s{seed}{form}.fa"
                    o = out / "hits" / f"{band}_{arm}_s{seed}{form or '_mates'}_{level}.tsv"
                    o.parent.mkdir(parents=True, exist_ok=True)
                    if o.exists() and o.stat().st_size >= 0 and (o.with_suffix(".done")).exists():
                        continue
                    if not q.exists():
                        raise SystemExit(f"{q} missing; run --stages reads first")
                    if diamond:
                        # SENSITIVITY, and why it is not the maximum. results/reach/ ran
                        # eleven search methods over one panel and DIAMOND at sensitive,
                        # very-sensitive and ultra-sensitive, BLASTP and full
                        # Smith-Waterman all returned the SAME nearest relative on all 16
                        # contigs: the binding constraint is the significance statistics,
                        # not the seeding. --ultra-sensitive measured at 182 s per 2,000
                        # queries here against ~40 s for --very-sensitive, which is 3.8
                        # hours against 50 minutes over this grid, for a difference that
                        # directory could not detect. --max-target-seqs matches the BLAST
                        # side so the per-pair subject count, which the Fisher correction
                        # uses, is not set by a tool default.
                        cmd = [str(blastbin / "diamond"), "blastx", "-q", str(q),
                               "-d", str(dbp), "-e", evalue, "--max-target-seqs", "50",
                               "-p", str(threads), "-o", str(o), f"--{sensitivity}",
                               "--quiet"]
                        cmd += outfmt_args(DEFAULT_COLUMNS, diamond=True)
                    else:
                        cmd = [str(blastbin / prog), "-query", str(q), "-db", str(dbp),
                               "-evalue", evalue, "-max_target_seqs", "50",
                               "-num_threads", str(threads), "-out", str(o)]
                        if prog == "blastn":
                            cmd += ["-task", "blastn"]  # most sensitive nucleotide setting
                        cmd += outfmt_args(DEFAULT_COLUMNS, diamond=False)
                    t = time.perf_counter()
                    subprocess.run(cmd, check=True, capture_output=True)
                    o.with_suffix(".done").write_text("ok\n")
                    dt = round(time.perf_counter() - t, 1)
                    timings[o.name] = dt
                    log(f"  searched {o.name} in {dt}s")
    return timings


# --------------------------------------------------------------------- statistics
def e_to_p(e: float) -> float:
    """An E-value read as a p-value: p = 1 - exp(-E), exact under the Poisson model.

    For the small E-values that matter this is E itself; the exact form is used so that a
    mate with E = 3 does not contribute a "p-value" above 1 and quietly break Fisher.
    """
    if e <= 0:
        return 1e-300
    return 1.0 - math.exp(-min(e, 700.0))


def fisher4(p1: float, p2: float) -> float:
    """Fisher's combined p for two independent tests: chi-square, 4 df, closed form.

    sf(x; 4) = exp(-x/2) * (1 + x/2). Even degrees of freedom have an exact elementary
    survival function, so no scipy and no table.
    """
    p1 = min(max(p1, 1e-300), 1.0)
    p2 = min(max(p2, 1e-300), 1.0)
    x = -2.0 * (math.log(p1) + math.log(p2))
    return math.exp(-x / 2.0) * (1.0 + x / 2.0)


def pair_modes(mate_rows, cat_rows, catn_rows, geom) -> dict:
    """Every mode's best statistic for one pair, plus the junction census.

    Returns a dict of mode -> (statistic, subject), where the statistic is on the E-value
    scale for every mode so the same threshold means the same thing across modes. The
    `fisher` statistic is a combined p multiplied by the number of subjects the PAIR
    tried, which is the pair's own multiple-testing burden; a looser correction would make
    this mode look better for a reason that has nothing to do with the pair.
    """
    out = {}
    # 1. best_mate: the production unit
    best = min((r["evalue"], r["sseqid"], r["qseqid"]) for r in mate_rows) \
        if mate_rows else None
    out["best_mate"] = (best[0], best[1]) if best else (float("inf"), "")

    # 2/3. the concatenations, with the junction classified
    junction = Counter()
    for label, rows in (("concat", cat_rows), ("concat_n", catn_rows)):
        j = geom["cat_junction"] if label == "concat" else geom["catn_junction"]
        off = 0 if label == "concat" else geom["catn_spacer"]
        b = None
        for r in rows:
            qs, qe = sorted((r["qstart"], r["qend"]))
            if qe <= j:
                where = "mate1"
            elif qs > j + off:
                where = "mate2"
            else:
                where = "crosses"
            junction[f"{label}:{where}"] += 1
            if b is None or r["evalue"] < b[0]:
                b = (r["evalue"], r["sseqid"], where)
        out[label] = (b[0], b[1]) if b else (float("inf"), "")
        if b:
            junction[f"{label}:best_{b[2]}"] += 1

    # 4. fisher: per subject, the two mates' own best HSPs against that subject
    per_sub = defaultdict(dict)
    for r in mate_rows:
        mate = r["qseqid"].rsplit("/", 1)[-1]
        cur = per_sub[r["sseqid"]].get(mate)
        if cur is None or r["evalue"] < cur:
            per_sub[r["sseqid"]][mate] = r["evalue"]
    both = {s: d for s, d in per_sub.items() if len(d) == 2}
    ntried = max(len(per_sub), 1)
    fb = (float("inf"), "")
    for s, d in both.items():
        e1, e2 = d["1"], d["2"]
        comb = fisher4(e_to_p(e1), e_to_p(e2)) * ntried
        if comb < fb[0]:
            fb = (comb, s)
    out["fisher"] = fb

    # 5. agree: both mates' own best hit names the same subject. The statistic reported is
    # the WORSE of the two mates' E-values, so a threshold on it is "each mate at least
    # this good", which is how a consistency rule is actually applied.
    b1 = min((r["evalue"], r["sseqid"]) for r in mate_rows
             if r["qseqid"].endswith("/1")) if any(r["qseqid"].endswith("/1")
                                                   for r in mate_rows) else None
    b2 = min((r["evalue"], r["sseqid"]) for r in mate_rows
             if r["qseqid"].endswith("/2")) if any(r["qseqid"].endswith("/2")
                                                   for r in mate_rows) else None
    if b1 and b2 and b1[1] == b2[1]:
        out["agree"] = (max(b1[0], b2[0]), b1[1])
    else:
        out["agree"] = (float("inf"), "")
    # 6. best_or_fisher: the only form in which a joint statistic is safe to deploy.
    # `fisher` on its own DILUTES a strong single mate with a useless partner (a mate at
    # E = 1e-9 beside one at E = 5 comes out at 1.5e-4), so as a replacement it loses more
    # than it gains. Taken as the better of the two it can only add, and the question
    # becomes what it costs in the null rather than what it costs in sensitivity.
    out["best_or_fisher"] = (min(out["best_mate"], out["fisher"], key=lambda t: t[0]))

    out["_census"] = dict(junction=dict(junction), subjects_tried=ntried,
                          subjects_both_mates=len(both),
                          mate1_hits=sum(1 for r in mate_rows if r["qseqid"].endswith("/1")),
                          mate2_hits=sum(1 for r in mate_rows if r["qseqid"].endswith("/2")))
    return out


MODES = ("best_mate", "concat", "concat_n", "fisher", "agree", "best_or_fisher")


def stage_report(out: Path, seeds, band: str, arms, fixed: float, levels: dict) -> dict:
    """Per-pair statistics, then each mode's rate at the fixed and matched-null points."""
    rows_path = out / f"pairs_{band}.tsv.gz"
    stats = {level: defaultdict(lambda: defaultdict(list)) for level in levels}
    junction = {level: Counter() for level in levels}
    with gzip.open(rows_path, "wt") as w:
        w.write("band\tlevel\tarm\tseed\tpair\t" + "\t".join(
            f"{m}_e\t{m}_subject" for m in MODES) + "\tsubjects_both\tsubjects_tried\n")
        for level, (prog, _) in levels.items():
            for seed in seeds:
                for arm in arms:
                    geom = {}
                    gp = out / "queries" / f"{band}_{arm}_s{seed}_geom.tsv"
                    if not gp.exists():
                        raise SystemExit(f"{gp} missing; the geometry is written with the "
                                         f"queries and a concatenated HSP cannot be "
                                         f"placed without it")
                    with open(gp) as fh:
                        next(fh)
                        for line in fh:
                            f = line.rstrip("\n").split("\t")
                            geom[f[0]] = dict(len1=int(f[1]), len2=int(f[2]),
                                              cat_junction=int(f[3]),
                                              catn_junction=int(f[4]),
                                              catn_spacer=int(f[5]))
                    def grouped(form):
                        p = out / "hits" / f"{band}_{arm}_s{seed}{form}_{level}.tsv"
                        g = defaultdict(list)
                        for r in read_hits(p, program=prog):
                            g[r["qseqid"].rsplit("/", 1)[0]].append(r)
                        return g
                    gm, gc, gn = grouped("_mates"), grouped("_cat"), grouped("_catn")
                    for pair in geom:
                        m = pair_modes(gm.get(pair, []), gc.get(pair, []),
                                       gn.get(pair, []), geom[pair])
                        junction[level].update(m["_census"]["junction"])
                        cells = []
                        for mode in MODES:
                            e, s = m[mode]
                            cells += [f"{e:.3g}" if e != float("inf") else "inf", s]
                            stats[level][arm][mode].append(e)
                        w.write(f"{band}\t{level}\t{arm}\ts{seed}\t{pair}\t"
                                + "\t".join(cells)
                                + f"\t{m['_census']['subjects_both_mates']}"
                                  f"\t{m['_census']['subjects_tried']}\n")

    summary = {"band": band, "fixed_threshold": fixed, "levels": {}}
    for level in levels:
        per_arm = stats[level]
        # The matched-null point: the strictest threshold at which the shuffled arm is
        # still at zero. Taken as the shuffled arm's own best statistic, so a mode whose
        # null reaches further gets a stricter threshold and cannot buy sensitivity with
        # it. This is the whole guard against a joint statistic that gains by being read
        # at a threshold nobody calibrated.
        null = per_arm.get("shuffled", {})
        matched = {}
        for mode in MODES:
            vals = [v for v in null.get(mode, []) if v != float("inf")]
            matched[mode] = min(vals) if vals else None
        lv = {"matched_null_threshold": {k: (f"{v:.3g}" if v else None)
                                         for k, v in matched.items()},
              "null_pairs": len(null.get("best_mate", [])),
              "junction_census": dict(junction[level]), "arms": {}}
        for arm, modes in per_arm.items():
            n = len(modes["best_mate"])
            a = {"pairs": n, "modes": {}}
            for mode in MODES:
                v = modes[mode]
                at_fixed = sum(1 for x in v if x < fixed)
                thr = matched[mode]
                at_matched = sum(1 for x in v if thr is not None and x < thr)
                a["modes"][mode] = {
                    "pct_at_fixed": round(100.0 * at_fixed / n, 3) if n else None,
                    "pct_at_matched_null": (round(100.0 * at_matched / n, 3)
                                            if n and thr is not None else None),
                }
            lv["arms"][arm] = a
        summary["levels"][level] = lv
    tag = "_".join(sorted(levels))
    (out / f"summary_{band}.json").write_text(json.dumps(summary, indent=1))
    (out / f"summary_{band}_{tag}.json").write_text(json.dumps(summary, indent=1))
    return summary


def per_seed_table(out: Path, band: str, seeds, arms, fixed: float) -> None:
    """The same rates per seed, so the report can quote a range rather than a mean.

    Five seeds here are five independent simulations off ONE genome, so the spread is
    simulation noise and nothing else - not library and not biological variation. It is
    narrower than either by an unknown factor, which is the same caveat
    results/holdout_run/ carries.
    """
    import csv
    rows = []
    with gzip.open(out / f"pairs_{band}.tsv.gz", "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        acc = defaultdict(lambda: defaultdict(list))
        for r in rd:
            key = (r["level"], r["arm"], r["seed"])
            for mode in MODES:
                e = r[f"{mode}_e"]
                acc[key][mode].append(float("inf") if e == "inf" else float(e))
    for (level, arm, seed), modes in sorted(acc.items()):
        n = len(modes["best_mate"])
        row = dict(band=band, level=level, arm=arm, seed=seed, pairs=n)
        for mode in MODES:
            row[mode] = round(100.0 * sum(1 for x in modes[mode] if x < fixed) / n, 3)
        rows.append(row)
    with open(out / f"per_seed_{band}.tsv", "w") as w:
        w.write("\t".join(rows[0].keys()) + "\n")
        for r in rows:
            w.write("\t".join(str(v) for v in r.values()) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, nargs="+", default=[1, 2, 3, 4, 5], dest="seeds")
    ap.add_argument("--band", default="near", choices=sorted(BANDS))
    ap.add_argument("--stages", default="reads,search,report")
    ap.add_argument("--out", type=Path, default=ROOT / "results/pair_search")
    ap.add_argument("--blastbin", type=Path, required=True)
    ap.add_argument("--dbdir", type=Path, required=True,
                    help="directory holding db_holdout_nt and db_holdout_prot, built by "
                         "experiments/holdout_indexes.py")
    ap.add_argument("--candfa", type=Path,
                    default=ROOT / "results/holdout_candidates/candidates.fa")
    ap.add_argument("--pairs", type=int, default=1000)
    ap.add_argument("--read-len", type=int, default=145,
                    help="the holdout library's own read length, so the best_mate rate "
                         "is readable against results/holdout_run/")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--fastq", type=Path, default=None,
                    help="the holdout candidate's containing library (SRR23421544), from "
                         "which the real_mate1 arms take mate 1")
    ap.add_argument("--berlin-r1", type=Path, default=None)
    ap.add_argument("--berlin-r2", type=Path, default=None)
    ap.add_argument("--label-classes", type=Path,
                    default=Path("/mnt/project-files/shared/resasm/"
                                 "library_label_classes.npz"),
                    help="the committed whole-library label table, by spot")
    ap.add_argument("--arms", nargs="+", default=sorted(ARMS) + ["shuffled"],
                    choices=sorted(ARMS) + ["shuffled"])
    ap.add_argument("--search-evalue", default="10",
                    help="permissive floor for the search itself; every mode's threshold "
                         "is applied in the report stage")
    ap.add_argument("--levels", nargs="+", default=sorted(LEVELS),
                    choices=sorted(LEVELS) + sorted(BLASTX_LEVEL) + sorted(SPROT_LEVEL),
                    help="which search levels to run. The default pair is the five-seed "
                         "grid; `prot` adds the BLASTX cross-check, which is ~0.7 s per "
                         "query and is meant for one seed")
    ap.add_argument("--diamond-sensitivity", default="very-sensitive",
                    choices=["sensitive", "very-sensitive", "ultra-sensitive"])
    ap.add_argument("--fixed", type=float, default=1e-5,
                    help="this repo's reporting operating point")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    stages = a.stages.split(",")
    arms = list(a.arms)
    allv = dict(LEVELS, **BLASTX_LEVEL, **SPROT_LEVEL)
    levels = {k: allv[k] for k in a.levels}
    timings, t0 = {}, time.perf_counter()

    if "reads" in stages:
        t = time.perf_counter()
        stage_reads(a.out, load_genomes(a.candfa), a.seeds, a.band, a.pairs, a.read_len,
                    arms, a.fastq, a.candfa,
                    (a.berlin_r1, a.berlin_r2) if a.berlin_r1 else None,
                    a.label_classes)
        timings["reads"] = round(time.perf_counter() - t, 1)
    if "search" in stages:
        t = time.perf_counter()
        timings["search_per_file"] = stage_search(a.out, a.blastbin, a.dbdir, a.seeds,
                                                  a.band, arms, a.threads,
                                                  a.search_evalue, levels,
                                                  a.diamond_sensitivity)
        timings["search"] = round(time.perf_counter() - t, 1)
    if "positional" in stages:
        t = time.perf_counter()
        stage_positional(a.out, a.blastbin, a.dbdir, load_genomes(a.candfa),
                         ARMS["candidate"]["acc"], a.read_len)
        timings["positional"] = round(time.perf_counter() - t, 1)
    if "report" in stages:
        t = time.perf_counter()
        s = stage_report(a.out, a.seeds, a.band, arms, a.fixed, levels)
        per_seed_table(a.out, a.band, a.seeds, arms, a.fixed)
        timings["report"] = round(time.perf_counter() - t, 1)
        for level, lv in s["levels"].items():
            for arm, d in lv["arms"].items():
                log(f"{a.band} {level} {arm}: " + " ".join(
                    f"{m}={d['modes'][m]['pct_at_fixed']}" for m in MODES))
    timings["total"] = round(time.perf_counter() - t0, 1)
    man = dict(provenance(__file__), band=a.band, seeds=a.seeds, pairs=a.pairs,
               read_len=a.read_len, bands=BANDS, arms={k: ARMS[k] for k in arms
                                                        if k in ARMS},
               library=str(a.fastq) if a.fastq else None, insert_source="PR #10",
               error_rates=dict(r1=ERR1, r2=ERR2), search_evalue=a.search_evalue,
               fixed=a.fixed, levels={k: v[0] for k, v in levels.items()},
               diamond_sensitivity=a.diamond_sensitivity,
               timings=timings)
    (a.out / f"manifest_{a.band}.json").write_text(json.dumps(man, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
