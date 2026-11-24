"""Every stage a run has, declared: what it needs, what it writes, what reads it.

This is the registry's content. It describes the pipeline that exists today, including
the parts of it that are missing: a file the report reads with no stage behind it is
declared here as an orphan rather than left to be discovered by a reader of their own
water, which is how it was discovered the first time.

HOW THE ORDER IS DERIVED. Each stage names the stages it needs, and `Registry.order`
computes the sequence. The old `STAGE_ORDER` was a hand-kept list typed in three places,
and the launcher's copy is the one that ran: `panel` and `calibration` left production in
two copies and not the third, so a validation stage ran in a production pass. A list
cannot disagree with a derived order, because there is nothing to disagree with.

WHAT IS NOT HERE. Paths, databases, thread counts and seeds. Those are properties of a
machine and an input, they live in the run configuration, and a stage declaration that
carried one would be the shape of the defect this repo has unpicked twenty times.
"""
from __future__ import annotations

from .registry import Output, Registry, Stage

#: The report's sections, as the generator's one literal lists them
#: (experiments/make_web.py, `const SUMMARY`). Identical on every sample by design.
WATCHLIST = "Watchlist screen"
WATER = "What this water is"
DISTANCE = "How far the assembled viruses sit from anything published"
REMAINDER = "The unnamed remainder"
METHODS = "Methods"
SECTIONS = [WATCHLIST, WATER, DISTANCE, REMAINDER, METHODS]

R = Registry()

R.add(Stage(
    name="cascade",
    why="Subtract everything identifiable. This is the subtrahend and the whole first "
        "half of the thesis: the sample's own assembly removes 71.9% of one library "
        "where 16 GB of Kraken2 removes 0.66%, so which index is consulted beats how "
        "fast it is consulted by orders of magnitude.",
    script="subtract_known.py",
    outputs=(Output("subtract_s{seed}.json", "cascade"),
             Output("residue_rows_s{seed}.jsonl", "residue_rows"),
             Output("tiers.json", "tiers")),
))

R.add(Stage(
    name="contigs",
    why="The contigs that assembled and could not be identified: tier 2, about 8.4% of "
        "pairs at median depth 6. Kept as its own stage because subtracting "
        "unidentified contigs inflated removal to 85.9% against the true 77.0% - the "
        "subtrahend is identification, not assembly.",
    needs=("cascade",),
    outputs=(Output("contigs.jsonl", "contigs"),),
))

R.add(Stage(
    name="reads",
    why="Characterise the residual reads field by field, each beside its shuffled null "
        "and its known-positive arm, so a field that reads the same on both is reported "
        "with that fact rather than hidden.",
    script="residue_report.py",
    needs=("cascade",),
    outputs=(Output("reads.jsonl", "reads"),
             Output("fields.json", "fields")),
))

R.add(Stage(
    name="linkage",
    why="Whether a residual read's mate is also unplaced. Both-mates-unwitnessed beats "
        "exactly-one at matched depth, which is what makes the residue missing molecules "
        "rather than errors.",
    script="mate_friends.py",
    needs=("reads",),
    outputs=(Output("linkage.json", "linkage"),),
))

R.add(Stage(
    name="census",
    why="An explicit census of what was searched and found nothing. An empty table and a "
        "search that never ran look identical, and that cost a withdrawn claim.",
    script="linkage_census.py",
    needs=("reads",),
    outputs=(Output("census.json", "census", section=WATCHLIST),),
))

R.add(Stage(
    name="intervals",
    why="The cascade's tiers with intervals that account for the reads being a "
        "contiguous block of one flow cell rather than a random sample. Wilson alone "
        "understates it by the measured design effect.",
    needs=("cascade",),
    outputs=(Output("tier_table.json", "tier_table", section=REMAINDER),),
))

R.add(Stage(
    name="layers",
    why="Which layer each residual read sits in, and the fields that decided it, so a "
        "layer whose fields were only ever computed on the residue cannot be compared "
        "against the known pile and says so.",
    needs=("reads",),
    outputs=(Output("layers.json", "layers", section=WATCHLIST),
             Output("layers.jsonl", "layer_rows")),
))

R.add(Stage(
    name="strains",
    why="Novel strains of known viruses: the middle output class, where the nearest "
        "relative has a name and the distance is real.",
    needs=("reads",),
    outputs=(Output("strains.json", "strains", section=DISTANCE),),
))

R.add(Stage(
    name="divergent",
    why="The divergent-but-related contigs, ranked on the conserved gene rather than on "
        "a whole-contig identity. One contig's own genes differ by a median 17.4 points "
        "against a near-reference control's 2.35, so a whole-contig number reports "
        "whichever gene won on bitscore.",
    script="divergent_known.py",
    needs=("contigs",),
    outputs=(Output("divergent_known_s{seed}.json", "divergent"),),
))

R.add(Stage(
    name="shortlist",
    why="The removed pile as cells with counts on them rather than a list ordered by "
        "identity. Ordering a shortlist by identity is actively wrong: the holdout "
        "positive control sits at the 93.46th identity percentile of the pile that "
        "removed it, with 7,104 contigs ahead of it, and coverage is what works.",
    script="removal_shortlist.py",
    needs=("cascade",),
    outputs=(Output("removal_shortlist.json", "shortlist"),),
))

R.add(Stage(
    name="crossplant",
    why="Whether the residue recurs in other plants' background. Depth decides and the "
        "pile label almost nothing: at depth 1 residue survives 79.0% of a 512-library "
        "panel and identified reads 70.5%, and that 8.5-point gap is the whole "
        "known-against-unknown difference at matched depth.",
    script="crossplant_stage.py",
    needs=("reads",),
    outputs=(Output("crossplant_stage.json", "crossplant"),),
))

R.add(Stage(
    name="profile",
    why="The watchlist screen and the composition. Adopted into the run on the branch "
        "that found the page reading an absent file for it (PLAN.md 1.3).",
    script="public_health_profile.py",
    needs=("cascade",),
    outputs=(Output("public_health_profile_s{seed}.json", "profile", section=WATCHLIST),
             Output("public_health_profile_s{seed}.json", "profile", section=WATER)),
))

R.add(Stage(
    name="taxa",
    why="The distance line on the page: how far each assembled virus sits from the "
        "nearest thing published. Distance without a nearest-relative identity is half "
        "an answer, and the combined signal is new AND like X.",
    script="known_strain_divergence.py",
    needs=("contigs",),
    outputs=(Output("known_strain_divergence_s{seed}.json", "taxa", section=DISTANCE),),
))

R.add(Stage(
    name="ranking",
    why="The unnamed remainder, ranked. Ranked on occurrence and shape, never on a "
        "search score: an index can say whether something is already known and cannot "
        "say whether it is new, because score ordering tracks database proximity.",
    script="concern_ranking.py",
    needs=("reads", "strains"),
    outputs=(Output("ranking.json", "ranking", section=REMAINDER),),
))


# ------------------------------------------------- the enrichments, adopted in step 6
#
# These three were the orphans that writing this file down found: files the report reads
# which no stage wrote. They are NOT the same kind of gap as `profile`, `taxa` and
# `ranking`, and the difference decides what the report does about them. Those three had
# no producer and their sections rendered a census of an absent file, which tells a
# reader of their own water that something is missing from it. These three have no
# producer and the sections that read them render WITHOUT them, saying what that costs.
#
# So they are declared, they run, and their outputs are optional with the consequence of
# each absence stated. A run with no Pfam database still produces a report; it produces a
# report that says the distance panel ranked on the whole-contig hit.

R.add(Stage(
    name="workup",
    in_driver=False,
    why="The per-gene workup of the divergent end. A whole-contig identity is not a "
        "property of the contig: it reports whichever gene won on bitscore, and one "
        "contig's own genes differ by a median 17.4 points against a near-reference "
        "control's 2.35. This is what lets the distance panel rank on the conserved gene "
        "instead.",
    script="divergent_contig_workup.py",
    # `taxa`, not `strains`: it selects contigs out of the per-contig identity table
    # `known_strain_divergence.py` writes. The `strains` stage's own output answers to the
    # same word and holds a read-level subject table, which this script cannot read.
    needs=("taxa",),
    outputs=(Output("divergent_contig_workup_s{seed}.json", "workup", section=DISTANCE,
                    optional=True,
                    absence_means="the distance panel ranks on the whole-contig best hit "
                                  "and says so, which is the reading SPEC 2.2 calls not a "
                                  "property of the contig. It needs a Pfam HMM library "
                                  "and a Swiss-Prot protein index, which not every "
                                  "machine has staged."),),
))

R.add(Stage(
    name="architecture",
    in_driver=False,
    why="What the page is allowed to say an assembled thing IS, in plain English, from "
        "experiments/virus_vocabulary.json and the evidence that actually fired. It "
        "declines with a reason rather than ever recording 'not a phage', which is why "
        "the renderer has no vocabulary of its own and a wrong word is wrong in one "
        "place.",
    script="architecture_call.py",
    needs=("workup",),
    outputs=(Output("architecture_s{seed}.json", "architecture", section=DISTANCE,
                    optional=True,
                    absence_means="the page names no architecture for an assembled thing "
                                  "and makes no claim about what it is. It reads the "
                                  "per-gene workup, so it is absent wherever that is."),),
))

R.add(Stage(
    name="names",
    in_driver=False,
    why="Accession to organism name, read from the same protein FASTA the search index "
        "was built from rather than from a second source that could disagree with it. "
        "Imo's reading preference is names, not accessions.",
    script="subject_names.py",
    needs=("strains", "ranking"),
    outputs=(Output("subject_names.json", "names", optional=True,
                    absence_means="the strains panel and the sampled read table show the "
                                  "accession alone rather than inventing a label. It "
                                  "needs the FASTA the protein index was built from, "
                                  "which a run given only a .dmnd does not have."),),
))

# ---------------------------------------------------------------- validation, not output
#
# Imo, 2026-11-23: "why would we, with unknown wastewater, have a genome with a distance
# fixed before the run? by default we know nothing. we always know nothing." Both stages
# below score a run against sequence whose distance was fixed before it. Both are real and
# both are how this tool is known to work at all, and neither has an answer on water
# nobody has deposited from, which is every sample production is for. So they are declared
# as validation and are not in the default set, and no report section names their output.

R.add(Stage(
    name="calibration",
    why="The per-sample divergence ladder: a known positive at a chosen distance. The "
        "only positive control the discriminating dataset has, since no depositor truth "
        "panel exists for that water.",
    script="spike_holdout.py",
    needs=("cascade",),
    validation=True,
    outputs=(Output("calibration_ladder.json", "calibration"),),
))

R.add(Stage(
    name="panel",
    why="Recovery against the depositors' own assemblies out of the same water. 25 "
        "genomes at 0.4-700x with at least six near-identical pairs, which is also why "
        "purity is scored against the near-identical group and not the genome.",
    script="truth_panel.py",
    needs=("cascade",),
    validation=True,
    outputs=(Output("truth_panel.json", "panel"),),
))

R.add(Stage(
    name="report",
    why="The page itself. It is a stage because its inputs are declared like any other "
        "stage's, which is what lets a section with no producer be a build error instead "
        "of a page that tells a reader something is missing from their sample.",
    script="make_web.py",
    # The package's, not the driver's: the driver never built the page, which is why a run
    # could finish with every stage green and no report. Marked rather than special-cased
    # at each use, because `--stages report` to the driver is an argparse error.
    in_driver=False,
    needs=("intervals", "layers", "census", "strains", "ranking", "profile", "taxa",
           "workup", "architecture", "names"),
    outputs=(Output("report.html", "report", section=METHODS),
             Output("render_check.json", "render_check")),
))


#: Files the report reads out of a run directory, and the stage that writes each.
#:
#: WHY THIS IS WRITTEN DOWN. A report section with no stage behind it is a hole in the
#: tool and not in the water, and the way that hole gets made is a renderer that reads a
#: file nobody writes: `load()` returns None for an absent file by design, so the page
#: prints "absent" and the build succeeds. PLAN.md 1.3 named three of these. Declaring
#: the whole set found three more, below.
#:
#: `tests/test_registry.py` checks this map against what `make_web.py` actually reads, so
#: a new read added to the renderer without a producer fails the suite rather than
#: rendering as an absence.
REPORT_INPUTS = {
    "manifest.json": "the run itself",
    "blind_spots.json": "the run itself",
    "tier_table.json": "intervals",
    "layers.json": "layers",
    "layers.jsonl": "layers",
    "fields.json": "reads",
    "reads.jsonl": "reads",
    "linkage.json": "linkage",
    "strains.json": "strains",
    "census.json": "census",
    "ranking.json": "ranking",
    "public_health_profile_s*.json": "profile",
    "known_strain_divergence_s*.json": "taxa",
    "calibration_ladder.json": "calibration",
    "divergent_contig_workup_s*.json": "workup",
    "architecture_s*.json": "architecture",
    "subject_names.json": "names",
}

#: Files the report reads that no stage writes. Empty since PLAN.md step 6.
#:
#: It held six. Three - `public_health_profile_s*.json`, `known_strain_divergence_s*.json`
#: and `ranking.json` - were sections rendering a census of an absent file, and are stages
#: now. Three more, found by writing `REPORT_INPUTS` down, were enrichments whose sections
#: render without them: they are stages now too, and their outputs are declared optional
#: with the consequence of each absence stated, which is a different claim and has to look
#: like one.
#:
#: `tests/test_registry.py` asserts this is empty, so a report input added without a
#: producer fails the suite rather than rendering as an absence on somebody's page.
ORPHAN_REPORT_INPUTS: dict[str, str] = {}
