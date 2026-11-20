"""Which row of the peptide-absence headroom table applies to the real Berlin data.

`results/peptide_absence_noise/REPORT.md` swept the substitution rate because
SRR23973193's own quality profile was not on disk, and §7 recorded that as the gap which
decides the answer: at 0.001 the one-per-million budget needs every k-mer in a read to be
absent, and at 0.002 it is not reachable at all. The sweep is only useful once we know
which row we are on.

This measures that from the run's own quality strings. Reads are fetched from ENA over
HTTP, which needs no SRA toolkit, and only a prefix of each file is read because a quality
profile converges long before a million reads.

Two things are measured, and they answer different questions:

1. The mean error probability per base, converted from Phred. This is the number that
   indexes the headroom table. It is reported both over all bases and over the leading
   portion of the read, because Illumina error rises sharply toward the 3' end and a
   pipeline that quality-trims is operating at the lower figure.
2. The profile by position, so the shape can be compared with the flat and ramped models
   used in the noise-floor experiment rather than assumed to match either.

A caveat that travels with every number here: Phred scores are the instrument's estimate
of its own error rate, and they are known to be miscalibrated, usually optimistically.
These are therefore a lower bound on the true substitution rate, which matters because the
headroom table gets worse as the rate rises, not better.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENA = ("https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/093/SRR23973193/"
       "SRR23973193_{r}.fastq.gz")


def phred_to_p(q: int) -> float:
    return 10 ** (-q / 10)


def profile(path: Path, max_reads: int) -> dict:
    """Per-position error probability and read-length distribution."""
    pos_sum: list[float] = []
    pos_n: list[int] = []
    lens: dict[int, int] = {}
    n = 0
    with path.open() as fh:
        while n < max_reads:
            if not fh.readline():
                break
            fh.readline()
            fh.readline()
            qual = fh.readline().strip()
            if not qual:
                break
            lens[len(qual)] = lens.get(len(qual), 0) + 1
            if len(qual) > len(pos_sum):
                pos_sum.extend([0.0] * (len(qual) - len(pos_sum)))
                pos_n.extend([0] * (len(qual) - len(pos_n)))
            for i, ch in enumerate(qual):
                pos_sum[i] += phred_to_p(ord(ch) - 33)
                pos_n[i] += 1
            n += 1
    per_pos = [pos_sum[i] / pos_n[i] if pos_n[i] else 0.0 for i in range(len(pos_sum))]
    total_p = sum(pos_sum)
    total_n = sum(pos_n)
    modal_len = max(lens, key=lambda k: lens[k]) if lens else 0
    first_two_thirds = per_pos[:max(1, (modal_len * 2) // 3)]
    return {
        "reads": n,
        "modal_read_length": modal_len,
        "read_length_counts": dict(sorted(lens.items())[:10]),
        "mean_error_probability": total_p / total_n if total_n else 0.0,
        "mean_error_probability_first_two_thirds": (
            sum(first_two_thirds) / len(first_two_thirds) if first_two_thirds else 0.0),
        "per_position_error_probability": [round(x, 6) for x in per_pos],
        "first_position": per_pos[0] if per_pos else 0.0,
        "last_position": per_pos[-1] if per_pos else 0.0,
    }


def fetch_prefix(r: int, dest: Path, lines: int) -> Path:
    """Stream the first `lines` lines of a gzipped FASTQ without downloading it whole."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    import gzip
    import io
    req = urllib.request.Request(ENA.format(r=r),
                                 headers={"User-Agent": "hapax-experiment/1.0"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        gz = gzip.GzipFile(fileobj=io.BufferedReader(resp))  # type: ignore[arg-type]
        with dest.open("w") as out:
            for i, line in enumerate(gz):
                if i >= lines:
                    break
                out.write(line.decode())
    return dest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="unused: this is a deterministic read of a file prefix, but the "
                         "flag is kept so every experiment takes --seed as the repo requires")
    ap.add_argument("--max-reads", type=int, default=400000)
    ap.add_argument("--data", type=Path, default=ROOT / "data" / "berlin")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "berlin_quality")
    args = ap.parse_args()

    t0 = time.perf_counter()
    args.out.mkdir(parents=True, exist_ok=True)
    out = {}
    for r in (1, 2):
        path = fetch_prefix(r, args.data / f"sample_R{r}.fastq", args.max_reads * 4)
        out[f"R{r}"] = profile(path, args.max_reads)

    # The headroom table's rows, so the mapping is explicit rather than eyeballed.
    rows = [0.0005, 0.001, 0.002, 0.005, 0.01]
    for key in ("R1", "R2"):
        e = out[key]["mean_error_probability"]
        out[key]["nearest_headroom_row"] = min(rows, key=lambda x: abs(x - e))
        out[key]["above_highest_row"] = e > rows[-1]

    (args.out / "berlin_quality.json").write_text(json.dumps(out, indent=2))
    (args.out / "manifest.json").write_text(json.dumps({
        "script": "experiments/berlin_quality_profile.py",
        "accession": "SRR23973193",
        "source": "ENA fastq over HTTP",
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "git_sha": _git_sha(),
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "platform": platform.platform(),
    }, indent=2))

    for key in ("R1", "R2"):
        d = out[key]
        print(f"{key}: {d['reads']} reads, modal length {d['modal_read_length']}, "
              f"mean error {d['mean_error_probability']:.5f} "
              f"(first two thirds {d['mean_error_probability_first_two_thirds']:.5f}), "
              f"position 1 {d['first_position']:.5f} -> last {d['last_position']:.5f}, "
              f"nearest table row {d['nearest_headroom_row']}")


def _git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
