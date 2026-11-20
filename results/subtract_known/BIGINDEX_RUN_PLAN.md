# The big-index run: what it does, what it writes, what it costs

Draft plan for the one EC2 run Imo has offered hardware for. **Nothing launches before he
approves a figure**, and the figure's search component is blocked on
`experiments/diamond_scaling.py`, which is queued. Everything else here is measured or
quoted from his launch guide.

## 1. Why a box at all

Four things this container cannot do, in order of how much they are holding up:

1. **Requirement 1.iii** - "the sequences with no match are really out at <80% identity" -
   is currently answered against indexes that fit on a 15 GB disk, plus a 60-read-per-class
   sample sent to the public NCBI service. That sample is rate-limited, loses batches to
   malformed replies, and runs against a database that moves under us with no version to
   cite.
2. **The far end of the removal curve is extrapolated.** The protein side was still buying
   about 20 points per tenfold step at the largest index reachable here, so the curve has no
   measured flattening point.
3. **The production index** the spec's section 1 asks for does not exist yet.
4. A sibling thread's two open numbers both turn on how clean the residue is, and it cannot
   get cleaner here. See section 4.

## 2. What gets staged, and the part that is free

NCBI's preformatted BLAST databases are in `s3://ncbi-blast-databases`, us-east-1, as
finished volume files rather than archives, so the listed size is the on-disk size: copy and
search, no build step, and BLAST memory-maps them so the 64 GB RAM ceiling does not bind.
Measured from the pinned prefix `2026-10-22-01-05-01` (the value of the bucket's
`latest-dir` key on 2026-11-20, and the version this run will cite):

| database | size | why |
|---|---|---|
| `core_nt` | 318.8 GB | the exact database the remote ground truth queries, now ours and pinned |
| `nr` | 928.0 GB | settles UniRef100-versus-nr with data instead of an argument |

The protein production index still needs a build, because searching millions of reads wants
DIAMOND rather than BLAST: UniRef100 is 227,086,028 clusters (UniProt's API, 2026-11-20),
about 100 GB as a DIAMOND database by the 0.44 GB/M measured on the built 10% subset.
UniRef50 (41.2 M, ~18 GB) is the cheaper fallback if the benchmark says UniRef100's search
cost is out of range.

Disk: 1.6 TB covers core_nt, nr and a UniRef100 build with headroom, inside the 2 TB
ceiling.

## 3. What the run measures

1. **The identity distribution of residue reads against core_nt and against nr**, reported
   as a distribution with intervals rather than a hit rate, since "nearly everything hits
   something" is already known and useless. Thousands of reads per class across 5 seeds,
   replacing 60.
2. **UniRef100 against nr on the same reads.** If the identity distributions agree, nr buys
   nothing here and the cheaper index is right; if they differ, nr is the index. The
   distributions ship either way so the choice can be remade when the databases move.
3. **The far end of the removal curve**, measured rather than extrapolated.
4. **The spike ladder re-run against the big indexes**, which is the only thing that can
   say whether the ~10-20% divergence recognition wall is a property of divergence or of
   our index size.

## 4. What the run writes for the other thread, and why it is nearly free

The thread measuring whether the residue is structured has its exact 31-mer coverage of all
99,745 scorable pairs of block 1 committed, and every readout in its report re-runs locally
in minutes against any residue definition. So it needs no compute here, only **read-name
lists**: for each of the five blocks, the pairs still unidentified after core_nt, after nr,
and after both, one name per line, same format as the existing `residue_names` files. The
cascade computes those sets anyway, so emitting them is a file write; all five blocks cost
the same as one.

### 4.1 The exact specification, fixed before launch

Fixed once by the thread that will consume them, because changing the format after the box
has shut down means another launch and another approval:

- **No scorable-pairs filter.** Survivors out of the full 100,000, not out of its 99,745.
  The 255 pairs it drops for having under 40 scorable k-mers are its own to drop, and a
  pre-filtered list makes a pair dropped for that reason indistinguishable from one a rung
  removed. It filters after joining.
- **Cumulative, not deltas**, and each file states what it contains: the pairs still
  unidentified after the existing tier-1 cascade PLUS that rung. Deltas are derivable from
  cumulative sets; the reverse needs an anchor, and the anchor is what drifts.
- **Fifteen files**, one per block per rung: `residue_names_<rung>_s<block>.txt`, rung in
  `{core_nt, nr, core_nt_nr}`. One name per line exactly as the FASTQ carries it
  (`SRR23973193.<n>`), no mate suffix, no header, LF, deduplicated, order irrelevant.
- **Plus the rung-zero list** - current cascade, no big index - from the same run, even
  though it already exists. If it differs from the committed `residue_names_s*.txt` by even
  a few names, the blocks drifted, and that is better found in a file than in a result.
- **A JSON beside each file with the operating point**, which a name list cannot carry:
  identity or bitscore threshold, aligner and mode, database version and date, and the
  count. Names without a threshold are as unreadable as a percentage without one.
- **A sixteenth file if it is free**: block 1 at a MID threshold as well as the strict one.
  Every ladder in this project has shown thresholds dominating database size, and one
  operating point on the biggest index is where that measurement would stop.
- **Confirmed convention**: the big run reads the library through
  `berlin_io.read_pairs` - skip (seed-1) x 1,500,000, 100,000 pairs, mate 2 trimmed 14 nt -
  so the names join its coverage table. If throughput ever forces a different reader the
  names stay correct but the joined set changes, and the rung-zero list is what exposes it.

Two measurements turn on it, both currently conditional on a database known to be
incomplete:

- **Does the 1-in-100,000 neighbourhood wall move?** A bigger index makes the residue cloud
  smaller and cleaner. If the wall is a property of a genome's abundance it stays; if it is
  partly clutter from known material, it falls. The spec's downstream design depends on
  which.
- **Does `cov_step`'s power reach zero or only approach it?** The in-band ladder took it
  from 0.317 to 0.494 at the sample's own contigs, chance to within the interval but on
  1,592 pairs. A core_nt rung says whether the witnessed material is exhausted or merely
  mostly gone, with more pairs surviving to measure on.

## 5. Costs

Measured spot prices, us-east-1, 2026-11-20: `m6i.2xlarge` $0.10/h, `c6i.2xlarge` $0.13/h,
`r6i.2xlarge` $0.29/h. gp3 is $0.08/GB-month, so 1.6 TB for a day is about $4.30. Staging
318.8 GB in-region at gp3's 125 MB/s baseline is about 45 minutes; nr is about 2 hours.

The search component is the number that matters and it is not quoted here on purpose. The
only throughput measurement available is from 8,000-query runs, which are dominated by the
one-off database scan, so multiplying it to a sample's ~2 M residue pairs would overstate
it badly - the same two-points-are-not-a-curve error this project has already made once.
`experiments/diamond_scaling.py` measures query scaling and database scaling separately and
reports the production figure as the linear extrapolation it is, with the extrapolation
factor beside it.

## 6. Compliance checklist, from the launch guide

- us-east-1, spot one-time, no `MaxPrice`, `InstanceInitiatedShutdownBehavior="terminate"`,
  every tag present, `DryRun=True` first expecting `DryRunOperation`,
  `terminate_instances` in a `finally`.
- `shutdown -h +N &` as the first line of UserData; heartbeat log to S3 every few minutes,
  because console output is unreliable; one exit helper on every path that uploads the log
  and an exit-reason object then `shutdown -c; shutdown -h now`.
- `dnf install -y python3.12` and call `python3.12` by name. `test -f` every input, because
  `aws s3 cp --recursive` over an empty prefix exits 0.
- Writes go only under `projects/hapax/{results,logs,scratch}/`. **The bucket has no
  delete**, so the keys are fixed before launch and no probe writes under `results/`.
- S3 clients pass `endpoint_url="https://s3.us-east-1.amazonaws.com"`.
- Code is pulled from S3 at a pinned commit from a clean tree. **This is a property of
  the run, not an implementation detail**: because the tree is clean, every manifest the
  run writes carries a resolved SHA rather than a dirty one, so its numbers are
  reproducible from the start instead of being retrofitted later. That matters because
  almost nothing else in this project can currently make that claim - the repo rule gives
  reproducibility as a design property (a committed script plus a seed), which is weaker
  than having re-run from clean and got the same numbers. The sibling thread has done the
  latter once, by accident, when a provenance fix forced a re-run. The run plan therefore
  pins a commit rather than rsyncing the working tree, even though the working tree is
  what is in hand.
- Copy `tools/ec2_detect_census.py` and `src/amazebio/ec2.py` from `ramics/evo-training-check`
  rather than writing the launcher fresh.
- First chunk of work is checked against a measured baseline in the S3 log before the run is
  trusted, not only at the end.

ElasticBLAST was considered and rejected: its IAM policy needs `iam:CreateRole`,
`iam:CreateInstanceProfile` and `iam:PassRole` plus CloudFormation and AWS Batch, all denied
here, and it defaults to a 32 vCPU machine against an 8 vCPU ceiling. It exists to fan one
search across many nodes; this is one box and a sample.
