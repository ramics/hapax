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
    needs=("intervals", "layers", "census", "strains", "ranking", "profile", "taxa"),
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
}

#: Files the report reads that NO stage writes, with the script that would write each.
#:
#: Each is a section or a panel that renders an absence on every run, and none of them is
#: an absence in the water: the script exists in `experiments/` and was never adopted
#: into `hapax run`, which is the same fault as the three PLAN.md 1.3 named and was found
#: the same way, by writing the map above down.
#:
#: PLAN.md step 6 empties this. Until then `tests/test_registry.py` asserts it EXACTLY,
#: so an orphan that is adopted has to be removed from here and a new one fails the
#: suite. A census, not a hidden absence.
ORPHAN_REPORT_INPUTS = {
    "architecture_s*.json": "experiments/architecture_call.py, which decides what the "
                            "page is allowed to say an assembled thing is, from "
                            "virus_vocabulary.json and the evidence that fired",
    "subject_names.json": "experiments/subject_names.py, which maps an accession to an "
                          "organism name from the same FASTA the DIAMOND database was "
                          "built from. Absent, the page shows the accession alone",
    "divergent_contig_workup_s*.json": "experiments/divergent_contig_queries.py and "
                                       "conserved_positions.py, the per-gene workup the "
                                       "distance panel ranks on. Absent, the panel falls "
                                       "back to the whole-contig hit, which is the "
                                       "reading SPEC 2.2 says is not a property of the "
                                       "contig",
}
