"""Do any committed five-seed directories have seeds that are empty, identical or short?

WHY. The block stride bug left two of five seeds reading nothing while the directory still
looked complete, and this repo quotes a min-max over five seeds for nearly every number.
An empty or duplicated seed inside a published range is invisible in the range itself.
"""
import argparse, glob, hashlib, json, re, statistics
from collections import defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--seed", type=int, default=1)  # unused; the convention requires it
ap.add_argument("--out", type=Path, default=ROOT / "results/fail_loudly_audit")
args = ap.parse_args()
groups = defaultdict(dict)
for f in glob.glob(str(ROOT/"results/**/*"), recursive=True):
    p = Path(f)
    if not p.is_file():
        continue
    m = re.search(r"^(.*?)_?s(\d+)(\..+)?$", p.name)
    if not m or not m.group(1):
        continue
    key = (str(p.parent.relative_to(ROOT)), m.group(1) + (m.group(3) or ""))
    groups[key][int(m.group(2))] = p
# hapax-style run_s1/ directories: key on the per-seed manifest
for f in glob.glob(str(ROOT/"results/**/run_s*/*.json"), recursive=True):
    p = Path(f)
    sd = int(re.search(r"run_s(\d+)", str(p)).group(1))
    groups[(str(p.parent.parent.relative_to(ROOT)), "run_s*/"+p.name)][sd] = p
out = []
for (d, stem), per in sorted(groups.items()):
    if len(per) < 2:
        continue
    sizes = {s: p.stat().st_size for s, p in per.items()}
    hashes = {s: hashlib.sha256(p.read_bytes()).hexdigest()[:12] for s, p in per.items()}
    dup = defaultdict(list)
    for s, h in hashes.items():
        dup[h].append(s)
    med = statistics.median(sizes.values())
    flags = []
    empty = sorted(s for s, b in sizes.items() if b == 0)
    short = sorted(s for s, b in sizes.items() if b and med and b < 0.5 * med)
    ident = [sorted(v) for v in dup.values() if len(v) > 1]
    if empty: flags.append(f"empty seeds {empty}")
    if short: flags.append(f"under half median size: seeds {short}")
    if ident: flags.append(f"byte-identical: {ident}")
    if flags:
        out.append({"dir": d, "files": stem, "seeds": sorted(per),
                    "sizes": sizes, "flags": flags})
args.out.mkdir(parents=True, exist_ok=True)
(args.out / "seed_degeneracy.json").write_text(json.dumps(
    {"script": "experiments/seed_degeneracy_audit.py", "seed": args.seed,
     "groups_checked": len(groups), "flagged": out}, indent=2))
print(len(groups), "multi-seed file groups checked;", len(out), "flagged")
for o in out:
    print(" ", o["dir"], o["files"], o["flags"])
