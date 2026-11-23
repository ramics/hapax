# hapax

Finds the viruses in wastewater that sit far from anything published.

![The stages of a run, the three output classes, and the distance scale they share](docs/hapax-overview.svg)

hapax removes the reads a library's reference databases already explain, assembles the
remainder, and reports every contig as a distance to its closest published relative: which
relative, over how much of the contig, and at what significance.

| Class | Contents |
|---|---|
| **A** | Contigs at real distance from anything published |
| **B** | Novel strains of named viruses |
| **C** | Unplaced sequence, reported with the depth it was seen at |

Known viruses are recognised out to about 10% divergence, fading by 20%. Detection needs
abundance above roughly one in a million.

## Install

```
python3 -m pip install -r requirements.txt
python3 -m pytest
```

The tests run on the Python packages alone. A real run also drives BLAST, DIAMOND, Kraken2
and foldseek, pinned in [`TOOLS.md`](TOOLS.md).

## Run

```
python3 hapax.py run SAMPLE_1.fastq.gz SAMPLE_2.fastq.gz \
    --seed N --out RUNDIR --library ACCESSION --trim-r2 K \
    --contigs IDENTIFIED_RRNA.fa IDENTIFIED_OTHER.fa \
    --unidentified-contigs UNIDENTIFIED.fa \
    --identified-contigs IDENTIFIED_OTHER.fa \
    --rrna-contigs IDENTIFIED_RRNA.fa \
    --assembly ASSEMBLY.fa \
    --cross-sample-ref OTHER_SAMPLES.fa \
    --spacers SPACERS.fa

python3 experiments/make_web.py results/web/pages/SAMPLE.json   # run viewer
```

Each sample's own files are passed explicitly. Add `--plan` to print the stages and exit;
`python3 hapax.py run --help` lists every option. hapax checks all inputs, indexes and tools
up front and exits 2 if one is missing.

## More

- [`SPEC.md`](SPEC.md): the spec, with the bounds on every claim
- [`results/README.md`](results/README.md): every experiment, failures included, plus the
  two archived earlier directions (codon-aware assembly, strain separation)
- [`CLAUDE.md`](CLAUDE.md): how to work in the repo
