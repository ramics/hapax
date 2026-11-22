# A blind positive control for the cascade: one real virus, held out by deposition date

**Run:** `experiments/holdout_run.py`, seeds 1-5, `script_sha256`
`e976a8d39cf2767e6cec476dd9a5862fce0dfe6625225a4edf04b1ca8e4bcf25`, manifest
`manifest_s1.json`. The manifest's `git_sha` reads `888e4ed-dirty`; the dirt was the
manifest files it was writing, not the code, and the content hash above matches committed
`HEAD:experiments/holdout_run.py` byte for byte. One manifest covers all five seeds
because one invocation ran them all.

## What this answers, and what it does not

This is the project's first positive control that nobody designed the distance of. The
candidate is a real sequence deposited after a public wastewater RNA library that
contains it, so the indexes can be rolled back to the day before and asked what they
make of reads from a virus they have genuinely never seen. The alternative - a spike
mutated uniformly to a chosen divergence - sets the answer in advance and degrades
nucleotide and protein identity together, which natural divergence does not do.

It answers one half of SPEC: **what the identification-and-distance readout reports on a
single read, at a known true distance.** It does not run the cascade's removal tiers
(SILVA, Kraken2, the sample's own assembly are not staged here), so no removal percentage
appears below.

**Operating point for every number here:** BLAST+ 2.17.0, `blastn -task blastn` and
`blastx`, E < 1e-5, against RefSeq viral nucleotide (19,650 records) and the whole RefSeq
viral protein release (722,827). Reported identities are per read, best HSP by bitscore,
significant hits only. Units are **reads, not pairs**: this library is single-end in
practice (filed `PAIRED`, no second mate present), so mate merging, mate disagreement and
mate linkage - the three things SPEC §4.3 relies on - are all unavailable and untested
here.

## The four arms

| arm | what it is | what it calibrates |
|---|---|---|
| `leak` | reads from a sequence that is in the reference and predates the cutoff | the readout at zero distance; and a leak alarm - if the candidate read like this, the holdout leaked |
| `candidate` | reads from the held-out virus, 32,575 of them at 452x | the measurement |
| `background` | reads from the same library matching none of the study's deposits | what the readout says on ordinary material, the "already seen" baseline |
| `shuffled` | dinucleotide-preserving shuffles of the candidate reads themselves | what the readout says on sequence carrying no information |

Holdout: 439 of 19,650 nucleotide records and 19,493 of 722,827 proteins dropped for
being deposited after 2022-05-10, keyed on the **original GenBank accession** recovered
from each RefSeq record's provenance line, not on `createdate`. Nothing was dropped for
being undated. The held-out organism has no RefSeq entry of its own, so the holdout is
purely the date filter.

## The ladder

Identification rate and reported identity, five seeds, min-max, full index:

| arm | true distance | reads identified, nt | reads identified, protein | reported identity, nt | reported identity, protein |
|---|---|---|---|---|---|
| `leak` | 0.08% divergent | **56.65-59.30%** | 48.00-51.10% | 100.0 (truth 99.92) | 100.0 (truth 99.81) |
| `background` | unknown, ordinary | 10.85-14.10% | 5.05-5.35% | 80.00-81.38 | 55.81-58.70 |
| `candidate` | **30.3% divergent** | **5.25-7.00%** | 2.45-2.80% | 82.56-83.33 (truth 69.72) | 67.35-73.91 (truth 66.45) |
| `shuffled` | no information | **0.00%** | **0.00%** | - | - |

Three things in that table are worth more than the rest.

**The null is clean.** Zero of 2,000 shuffled reads clear E < 1e-5 at either level, on all
five seeds, at both indexes. The project's standing figure - "92.3% of shuffled reads get
a hit at all" - was taken without an E-value floor; with one, the floor is zero. So the
candidate's 5-7% is signal, not the tail of a null.

**A single 145 nt read survives 30% divergence about one time in sixteen.** Against
roughly one time in two at zero distance. That is the fourth wall measured on a real
virus rather than on a mutated spike, and it is the number to quote at anyone who expects
read-level identification to find novel material.

**Nucleotide beats protein at read level, at every distance including this one.** 5.25-7.00%
against 2.45-2.80% on the candidate. This is the reverse of the contig-level result
(protein identifies 26% of contigs, nucleotide 6%) and the reverse of what "the protein
tier is the backbone past 20% divergence" predicts. The protein tier's advantage is a
length effect: at 145 nt a blastx alignment is ~48 residues and there is not enough of it
to clear the threshold. Nothing here contradicts the contig-level finding - it says the
two live at different lengths, and SPEC should not carry the protein claim down to read
level.

## The bias, its direction, and the control that removes it

Because the true distance is known, the bias in a per-read identity readout can be
measured instead of argued about. It is large and it points the dangerous way: **a novel
thing looks closer to the reference than it is.**

| level | reported (median) | truth | bias | median query coverage |
|---|---|---|---|---|
| nucleotide | 82.56-83.33% | 69.72% | **+12.8 to +13.6 points** | 61-64% |
| protein | 67.35-73.91% | 66.45% | +0.9 to +7.5 points | 88-96% |

At zero distance there is no bias at all (`leak` reads 100.0 against truth 99.92, at 100%
coverage). The bias is created by distance, and coverage is the control that removes most
of it. Candidate reads, nucleotide level, seed 1, binned by coverage:

| query coverage | reads | median identity | bias |
|---|---|---|---|
| under 50% | 36 | 87.0% | **+17.2** |
| 50-70% | 53 | 82.8% | +13.1 |
| 70-90% | 23 | 82.2% | +12.5 |
| over 90% | 10 | 76.8% | **+7.1** |

The mechanism is plain: a read from a 30%-divergent genome aligns over its most conserved
stretch and reports that stretch's identity. The p90 of candidate nucleotide identity
reaches **92.00-96.88%** - the same genome yields individual reads that report up to ~97%
identity to a known relative while being 69.7% identical to it.

**So a per-read identity is not a distance.** It is the distance to the most conserved
fraction of that read, and without coverage beside it, it is unreadable. This is the
ground-truth version of the claim this repo already retracted once ("residue aligns 36-48
nt at ~89%"): the retraction was right, and here is the size of the error it avoided.
Even at over 90% coverage the readout still flatters by 7 points, so coverage narrows the
bias but does not close it - a distance belongs to a molecule, not a read.

## Did the holdout do anything?

Almost nothing, and that is the finding. Full index against held-out index differs by at
most 0.15 points of identification rate in any arm, and the candidate's reported identity
is identical to three decimal places. Hits to post-cutoff subjects are 0.00% everywhere
except `background` nucleotide, where one seed reads 1.38% on the full index and
**0.00% on the held-out one** - the only place the holdout bites, and it confirms the
mechanism works.

Read in the one-sided direction that is safe: measuring a candidate against today's
reference can only make it look *less* novel, so 30.3% divergence is a floor. Rolling the
reference back four and a half years did not move it, which means this candidate's
distance is not an artefact of recent deposition, and the `leak` arm reading 100% while
the candidate reads 83% means nothing leaked.

## What the first pass got wrong

The first complete pass of all 80 jobs was discarded. Its shuffled control read 3.35-5.45%
identified at 80.8-83.6% median identity, against the candidate's 5.25-7.00% at
82.56-83.33% - a null that reads like the test, which would have been published as "the
candidate is indistinguishable from noise" and is in fact a bug. 70.4% of the shuffled arm
was unshuffled candidate reads: the Altschul-Erikson walk dead-ends unless each vertex's
reserved last edge leads to the terminal vertex, that condition was missing, and the
fallback returned the input sequence. The function now returns `None` on failure and the
caller counts the drops into `reads_manifest.json` (0 failures, 0 returned inputs, all
five seeds), and the shuffled arm now shuffles the *same* subsample the candidate arm uses
rather than an independent draw. All 80 jobs were re-run, not the 20 shuffled ones, because
the script hash changed.

The tell was the one this repo has a rule for - **the control moved with the test** - and
it was visible in the summary table before any of it was understood.

## Caveats, in the order they would bite

1. **Five seeds here are five subsamples of one read pool**, not five blocks of a flow cell
   and not five libraries. The min-max above is subsampling spread and nothing else; it
   says nothing about library or biological variation, and it is narrower than the real
   thing by an unknown factor.
2. **One candidate, one library, one distance.** 30.3% divergence is one rung. The shape
   of the curve between `leak` and `candidate` is not measured here.
3. **No mates.** Everything SPEC §4.3 gains from pairs is absent, which means this run
   measures a weaker instrument than the production tool. A genuinely paired candidate
   exists but costs ~30 GB of staging per library.
4. **2,000 reads per arm per seed** out of 32,575 candidate reads available. The
   identification rates are proportions on that sample.
5. Reads, not pairs, throughout.

## Cost

Zero. Staged libraries, local BLAST, in-container. 92 minutes of search for the re-run
(`timings_s.search` 5,499 s) plus 22 s to extract the arms.

## Re-running this

The hit tables are committed gzipped (`hits_*.tsv.gz`, 2.2 MB against 20 MB raw). The
script reads `hits_*.tsv`, so `gunzip results/holdout_run/hits_*.tsv.gz` before
`--stages report`, or run `--stages reads,search,report` and let it rebuild them. The
`.done` resume markers have been removed for the same reason: a marker beside a gzipped
table would make the script skip a search whose output it then cannot read. The two
held-out BLAST databases and the arm FASTAs are gitignored and rebuilt by
`--stages holdout,reads`; the inputs they need are the three staged SRA libraries and the
RefSeq viral nucleotide and protein releases.
