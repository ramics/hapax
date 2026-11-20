"""Are mosaic joins really clustered at the first window boundary?

The first run with join features put 5 of 6 joins at position 300, which looked
like either genuine early switching or an unreliable lineage call near contig
ends. Both would matter. But the right null is not uniform over positions: a
mosaic contig of 765 nt has exactly one interior 300 nt boundary, so it can only
ever report a join at 300. The null is uniform over the interior boundaries each
contig actually has, which makes this a Poisson-binomial with p = 1/b per
contig - exact at these counts, so no simulation.

Run it on any rows file holding join observables; it is a lookup over data
already committed, which is the point of persisting it.
"""

from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

WINDOW = 300


def interior_boundaries(contig_len: int, window: int = WINDOW) -> int:
    """Boundaries between consecutive scored windows, which is where a switch can
    be reported. lineage_windows() steps by `window` while a full window fits."""
    n_windows = max(0, (contig_len - window) // window + 1)
    return max(0, n_windows - 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=Path, nargs="+", required=True)
    ap.add_argument("--window", type=int, default=WINDOW)
    args = ap.parse_args()

    joins = []
    for f in args.rows:
        for line in f.open():
            row = json.loads(line)
            for j in row.get("joins", []):
                joins.append((row["pair"], row["seed"], j["contig_len"], j["pos"]))
    if not joins:
        print("no joins recorded in these rows")
        return

    ps = [1.0 / b if (b := interior_boundaries(L, args.window)) else 0.0
          for _, _, L, _ in joins]
    obs = sum(1 for *_, pos in joins if pos == args.window)
    exp = sum(ps)

    if len(ps) <= 22:
        tail = 0.0
        for bits in product([0, 1], repeat=len(ps)):
            if sum(bits) >= obs:
                pr = 1.0
                for bit, p in zip(bits, ps):
                    pr *= p if bit else (1 - p)
                tail += pr
        tail_s = f"{tail:.4f}"
    else:  # normal approximation, adequate once there are this many joins
        import math
        mu = exp
        sd = math.sqrt(sum(p * (1 - p) for p in ps))
        z = (obs - 0.5 - mu) / sd if sd else 0.0
        tail_s = f"{0.5 * math.erfc(z / math.sqrt(2)):.4f} (normal approx)"

    bs = [interior_boundaries(L, args.window) for _, _, L, _ in joins]
    print(f"joins {len(joins)}  mean interior boundaries {sum(bs) / len(bs):.2f}")
    print(f"at the first boundary: observed {obs}, expected {exp:.2f} under "
          f"uniform-over-boundaries")
    print(f"P(>= observed | uniform) = {tail_s}")


if __name__ == "__main__":
    main()
