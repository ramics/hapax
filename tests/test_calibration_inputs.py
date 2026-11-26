"""The divergence ladder's two new inputs, as refusals and as measurements.

The ladder is a production stage rather than validation, and what made that possible is
that neither of the two things it needs from a person is needed any more. The spike genome
was named by hand - chosen once, as absent from one sample's assembly, and then carried to
whatever library the ladder was pointed at - and the cross-sample reference was required,
which refused the whole positive control on any library with no sibling samples.

So two properties are asserted here. A genome is CHOSEN against the assembly the run just
produced, because absence from one water is not a fact a reference set can carry; and a
candidate already resident in the sample is rejected with the count that rejected it,
because a ladder run against resident sequence is removed by the sample's own contigs at
every rung and reads as the tool working when it is the sample working.
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

mappy = pytest.importorskip("mappy", reason="the absence test maps windows with minimap2")

import spike_holdout  # noqa: E402


def sequence(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def fasta(path: Path, records: dict[str, str]) -> Path:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in records.items()))
    return path


def test_the_ladder_genome_is_chosen_as_the_one_absent_from_this_sample(tmp_path):
    resident = sequence(1, 4_000)
    clean = sequence(2, 4_000)
    pool = fasta(tmp_path / "pool.fa", {"resident": resident, "clean": clean})
    assembly = fasta(tmp_path / "asm.fa", {"contig_1": resident})
    index = mappy.Aligner(str(assembly), preset="sr")

    chosen, record = spike_holdout.pick_genome(pool, {"selfasm_identified": index},
                                               min_mlen=40)

    assert list(chosen) == ["clean"], record
    rejected = [c for c in record["candidates"] if c["genome"] == "resident"]
    assert rejected and rejected[0]["windows_mapping_this_sample"] > 0
    assert "rejected_because" in rejected[0]
    # The count that rejected it travels, so the record says how clean the chosen one is
    # rather than asserting it.
    taken = [c for c in record["candidates"] if c["genome"] == "clean"]
    assert taken and taken[0]["windows_mapping_this_sample"] == 0


def test_three_segments_of_one_virus_are_not_three_genomes(tmp_path):
    """Asked for three out of the real release in order, the first answer was segments L,
    M and S of one Pseudomonas phage: three accessions, one virus. The point of taking
    several is that one genome's recognition decay is that genome's, which three segments
    of the same organism cannot answer while looking as though they do.
    """
    pool = tmp_path / "pool.fa"
    pool.write_text(
        ">NC_000001.1 Example phage phiX segment L, complete sequence\n"
        + sequence(11, 3_000) + "\n"
        ">NC_000002.1 Example phage phiX segment M, complete sequence\n"
        + sequence(12, 3_000) + "\n"
        ">NC_000003.1 Another virus 1, complete genome\n"
        + sequence(13, 3_000) + "\n")
    assembly = fasta(tmp_path / "asm.fa", {"contig_1": sequence(14, 2_000)})
    index = mappy.Aligner(str(assembly), preset="sr")

    chosen, record = spike_holdout.pick_genome(pool, {"selfasm_identified": index},
                                               min_mlen=40, want=3)

    assert list(chosen) == ["NC_000001.1", "NC_000003.1"], record["candidates"]
    second = [c for c in record["candidates"] if c["genome"] == "NC_000002.1"][0]
    assert "another segment" in second["rejected_because"]
    assert second["windows_mapping_this_sample"] == 0, (
        "the second segment was rejected for being the same organism, not for mapping")


def test_the_organism_key_strips_the_segment_and_the_completeness():
    assert (spike_holdout.organism_of(
        "NC_042068.1 Pseudomonas phage phiNN segment L, complete sequence")
        == spike_holdout.organism_of(
            "NC_042070.1 Pseudomonas phage phiNN segment S, complete sequence"))
    assert (spike_holdout.organism_of("NC_001477.1 Dengue virus 1, complete genome")
            != spike_holdout.organism_of("NC_042068.1 Pseudomonas phage phiNN segment L"))


def test_a_pool_with_nothing_clean_is_a_refusal_with_its_reason(tmp_path):
    resident = sequence(3, 4_000)
    pool = fasta(tmp_path / "pool.fa", {"resident": resident})
    assembly = fasta(tmp_path / "asm.fa", {"contig_1": resident})
    index = mappy.Aligner(str(assembly), preset="sr")

    chosen, record = spike_holdout.pick_genome(pool, {"selfasm_identified": index},
                                               min_mlen=40)

    assert chosen == {}
    assert "failed_because" in record and len(record["failed_because"]) > 40


def test_a_genome_outside_the_length_band_is_not_considered(tmp_path):
    pool = fasta(tmp_path / "pool.fa", {"tiny": sequence(4, 300),
                                        "right": sequence(5, 3_000)})
    assembly = fasta(tmp_path / "asm.fa", {"contig_1": sequence(6, 2_000)})
    index = mappy.Aligner(str(assembly), preset="sr")

    chosen, record = spike_holdout.pick_genome(pool, {"selfasm_identified": index},
                                               min_mlen=40)

    assert list(chosen) == ["right"]
    assert [c["genome"] for c in record["candidates"]] == ["right"]


def test_the_two_mates_of_a_simulated_pair_get_their_own_lengths():
    """A pair-symmetric instrument cannot measure an asymmetric one.

    After the non-template prefix comes off, SRR23973193's mate 2 is 93 nt against mate
    1's 107, which is 62 scorable 31-mers against 76. Simulating both at one length
    removes from the ladder the asymmetry it exists to measure.
    """
    genome = sequence(7, 5_000)
    quantiles = [(0.5, 194), (0.75, 246), (0.95, 363), (0.99, 485), (1.0, 600)]
    pairs = spike_holdout.simulate_pairs(genome, 20, random.Random(1), 107, quantiles,
                                         0.0, 0.0, read_len2=93)
    assert {len(a) for _n, a, _b in pairs} == {107}
    assert {len(b) for _n, _a, b in pairs} == {93}
    # And the default is still one length for both, which is what the function did before
    # the argument existed.
    same = spike_holdout.simulate_pairs(genome, 20, random.Random(1), 107, quantiles,
                                        0.0, 0.0)
    assert {len(b) for _n, _a, b in same} == {107}
