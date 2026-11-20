"""Do incumbent assemblers keep two near-relative lineages apart? A de-collinearised ladder.

Why this exists. This is the baseline the strain-separation direction has to beat, and
the benchmark half of the project's deliverable. The archived two-lineage study
(experiments/codon_assembler/results/abundance_sweep, §5) measured pooling at a fixed
sub-threshold budget, with one genome pair per rung and a labeller that ignored strand.
This run fixes both design faults and asks the separation question directly:

* three distinct real pairs per rung (experiments/select_pairs.py,
  results/pair_ladder/ladder.json), so between-pair spread inside a rung can be set
  against the trend across rungs, instead of divergence being collinear with pair
  identity;
* the strand-aware window labeller audited in results/labeller_audit (0.0% false
  mosaics on either strand), so mosaic counts mean something;
* per-lineage coverage at 5x, 15x and 40x, so the regime where each lineage alone is
  assemblable is covered, not only the regime where nothing is;
* contigs are kept, so any better labeller is a re-score, not a re-run.

Reads are spike-only (no metagenomic background). That is the cheap arm: collapse
between two lineages happens in their shared k-mers, which a background does not supply.
Whether that shortcut holds is checked separately on a real background, not assumed.

Readouts are continuous, per lineage: breadth recovered by contigs that are labelled
*purely* that lineage ("separated breadth"), breadth by any contig, mosaic mass, longest
pure contig. Single-lineage controls run through MEGAHIT give the labeller's null on real
assembler output at every rung.

Usage: python experiments/pair_ladder.py --seed 1 --n-seeds 5 --tools megahit,metaspades,penguin,haploflow
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from labeller_audit import lineage_windows_stranded, read_fasta  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASM = Path(os.environ.get("HAPAX_TOOLS", "/home/user/mamba/envs/asm/bin"))
AUX = Path(os.environ.get("HAPAX_AUX", "/home/user/mamba/envs/aux2/bin"))
GENOMES = [ROOT / "data/pairs/rabies_200.fa", ROOT / "data/labeller_audit/rabies_pairs.fa"]
READ_LEN, INSERT, INSERT_SD = 150, 400, 40
# HiSeq 2500 at its native 150 nt. The archived runs used MSv3 (a 250 nt MiSeq profile)
# truncated to shorter reads; its 3'-end error rate left PenguiN, whose overlaps demand
# 99% identity, unable to build a single 1 kb contig even from one lineage at 60x, while
# HS25 reads give it the whole genome. A profile under which one tool cannot run at all
# is not a fair benchmark, so all tools here get HS25.
ART_PROFILE = "HS25"
MIN_CONTIG = 500
WINDOW = 300


def sh(cmd, **kw):
    env = dict(os.environ, PATH=f"{ASM}:{AUX}:" + os.environ["PATH"])
    return subprocess.run([str(c) for c in cmd], check=True, capture_output=True, text=True, env=env, **kw)


def simulate(seq_name: str, seq: str, cov: float, seed: int, out: Path) -> tuple[Path, Path]:
    out.mkdir(parents=True, exist_ok=True)
    fa = out / "g.fa"
    fa.write_text(f">{seq_name}\n{seq}\n")
    pref = out / "r"
    sh(["art_illumina", "-ss", ART_PROFILE, "-i", fa, "-p", "-l", READ_LEN, "-f", cov,
        "-m", INSERT, "-s", INSERT_SD, "-rs", seed, "-na", "-o", pref])
    return Path(f"{pref}1.fq"), Path(f"{pref}2.fq")


def cat(paths, out: Path) -> Path:
    with out.open("w") as fh:
        for p in paths:
            fh.write(p.read_text())
    return out


def assemble(tool: str, r1: Path, r2: Path, out: Path, threads: int) -> Path:
    if out.exists():
        shutil.rmtree(out)
    if tool == "megahit":
        sh(["megahit", "-1", r1, "-2", r2, "-o", out, "-t", threads, "--min-contig-len", MIN_CONTIG])
        return out / "final.contigs.fa"
    if tool == "metaspades":
        sh(["metaspades.py", "-1", r1, "-2", r2, "-o", out, "-t", threads, "-m", 6, "--only-assembler"])
        return out / "contigs.fasta"
    if tool in ("penguin", "penguin_guided"):
        out.mkdir(parents=True)
        mode = "guided_nuclassemble" if tool == "penguin_guided" else "nuclassemble"
        sh(["penguin", mode, r1, r2, out / "contigs.fa", out / "tmp", "--threads", threads])
        shutil.rmtree(out / "tmp", ignore_errors=True)
        return out / "contigs.fa"
    if tool == "haploflow":
        out.mkdir(parents=True)
        inter = cat([r1, r2], out / "reads.fq")
        sh(["haploflow", "--read-file", inter, "--out", out / "hf", "--log", out / "hf.log"])
        inter.unlink()
        return out / "hf" / "contigs.fa"
    raise ValueError(tool)


def score(contigs: Path, names: tuple[str, str], seqs: tuple[str, str], work: Path,
          spike_only_filter: bool = False) -> tuple[dict, list]:
    """Per-contig lineage labels and per-lineage breadth, separated and total.

    With a background, contigs that do not align to either lineage over >=80% of their
    length are dropped before labelling (the archived sweep's spike rule): a background
    contig scored against two rabies genomes gets an arbitrary nearest call, which would
    count as a lineage label it never earned.
    """
    cs = {k: v for k, v in read_fasta(contigs).items() if len(v) >= MIN_CONTIG} if contigs.exists() else {}
    if not cs:
        z = {"contigs": 0, "mass": 0}
        for s in "ab":
            z |= {f"breadth_{s}": 0.0, f"sep_breadth_{s}": 0.0, f"longest_pure_{s}": 0}
        return z | {"mosaic_contigs": 0, "mosaic_mass": 0, "unlabelled_mass": 0, "longest": 0}, []
    ref = work / "ref2.fa"
    ref.write_text(f">{names[0]}\n{seqs[0]}\n>{names[1]}\n{seqs[1]}\n")
    q = work / "q.fa"
    q.write_text("".join(f">{k}\n{v}\n" for k, v in cs.items()))
    paf = sh(["minimap2", "-x", "asm20", "-c", "--secondary=no", ref, q]).stdout
    if spike_only_filter:
        aligned: dict[str, int] = {}
        for line in paf.splitlines():
            f = line.split("\t")
            if int(f[10]) and int(f[9]) / int(f[10]) >= 0.85:
                aligned[f[0]] = aligned.get(f[0], 0) + int(f[3]) - int(f[2])
        cs = {k: v for k, v in cs.items() if aligned.get(k, 0) >= 0.8 * len(v)}
        if not cs:
            return score(work / "none.fa", names, seqs, work)
    per = []
    labels = {}
    for k, v in cs.items():
        wa, wb, sw, calls = lineage_windows_stranded(v, seqs[0], seqs[1], WINDOW)
        if wa and wb:
            lab = "mosaic"
        elif wa:
            lab = "a"
        elif wb:
            lab = "b"
        else:
            lab = "unlabelled"
        labels[k] = lab
        per.append({"contig": k, "len": len(v), "label": lab, "calls": "".join(calls), "switches": sw})
    cov = {"a": set(), "b": set()}
    sep = {"a": set(), "b": set()}
    for line in paf.splitlines():
        f = line.split("\t")
        if int(f[10]) == 0 or int(f[9]) / int(f[10]) < 0.85:
            continue
        side = "a" if f[5] == names[0] else "b"
        span = range(int(f[7]), int(f[8]))
        cov[side].update(span)
        if labels.get(f[0]) == side:
            sep[side].update(span)
    row = {"contigs": len(cs), "mass": sum(map(len, cs.values())),
           "longest": max(map(len, cs.values()))}
    for i, s in enumerate("ab"):
        row[f"breadth_{s}"] = round(len(cov[s]) / len(seqs[i]), 4)
        row[f"sep_breadth_{s}"] = round(len(sep[s]) / len(seqs[i]), 4)
        row[f"longest_pure_{s}"] = max((p["len"] for p in per if p["label"] == s), default=0)
    row["mosaic_contigs"] = sum(p["label"] == "mosaic" for p in per)
    row["mosaic_mass"] = sum(p["len"] for p in per if p["label"] == "mosaic")
    row["unlabelled_mass"] = sum(p["len"] for p in per if p["label"] == "unlabelled")
    return row, per


def one_job(job: dict, genomes: dict[str, str], work_root: Path, keep_root: Path, threads: int,
            background: tuple[Path, Path] | None = None) -> dict:
    pair, mode, cov, seed, tool = job["pair"], job["mode"], job["cov"], job["seed"], job["tool"]
    a, b = pair["a"], pair["b"]
    tag = f"{pair['tag']}_{a}_{b}_{mode}_c{cov}_s{seed}"
    # one directory per tool: parallel jobs must not share read files. ART is seeded,
    # so every tool still sees byte-identical reads for the same (pair, coverage, seed).
    work = work_root / f"{tag}_{tool}"
    work.mkdir(parents=True, exist_ok=True)
    # Reads depend on (pair, coverage, seed) only, so every tool sees identical reads.
    ra = simulate(a, genomes[a], cov, seed, work / "ra") if mode in ("mix", "a_alone") else None
    rb = simulate(b, genomes[b], cov, seed + 1000, work / "rb") if mode in ("mix", "b_alone") else None
    parts = [r for r in (ra, rb) if r]
    if background:
        parts.append(background)
    r1 = cat([p[0] for p in parts], work / "r1.fq") if not (work / "r1.fq").exists() else work / "r1.fq"
    r2 = cat([p[1] for p in parts], work / "r2.fq") if not (work / "r2.fq").exists() else work / "r2.fq"
    t0 = time.perf_counter()
    err = None
    try:
        contigs = assemble(tool, r1, r2, work / tool, threads)
    except subprocess.CalledProcessError as e:
        contigs, err = work / "none.fa", (e.stderr or "")[-400:]
    secs = time.perf_counter() - t0
    row, per = score(contigs, (a, b), (genomes[a], genomes[b]), work, spike_only_filter=bool(background))
    keep = keep_root / tag
    keep.mkdir(parents=True, exist_ok=True)
    if contigs.exists() and background:
        # a background assembly is ~60k contigs; keep only the spike contigs that were scored
        mine = {p["contig"] for p in per}
        with gzip.open(keep / f"{tool}.contigs.fa.gz", "wt") as fo:
            for k, v in read_fasta(contigs).items():
                if k in mine:
                    fo.write(f">{k}\n{v}\n")
    elif contigs.exists():
        with open(contigs, "rb") as fi, gzip.open(keep / f"{tool}.contigs.fa.gz", "wb") as fo:
            shutil.copyfileobj(fi, fo)
    # Per-contig labels go beside the kept contigs, not into the committed rows:
    # PenguiN's unreduced output alone is ~900 contigs per run.
    with gzip.open(keep / f"{tool}.labels.json.gz", "wt") as fo:
        json.dump(per, fo)
    shutil.rmtree(work, ignore_errors=True)
    return {"rung": pair["tag"], "a": a, "b": b, "divergence": pair["divergence"],
            "archived_pair": pair["archived"], "mode": mode, "cov_per_lineage": cov,
            "seed": seed, "tool": tool, "seconds": round(secs, 1), "error": err, **row}


def versions() -> dict:
    v = {"python": platform.python_version()}
    for name, cmd in {"megahit": ["megahit", "--version"], "spades": ["spades.py", "--version"],
                      "penguin": ["penguin", "version"], "minimap2": ["minimap2", "--version"],
                      "art": ["art_illumina", "--help"]}.items():
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               env=dict(os.environ, PATH=f"{ASM}:{AUX}:" + os.environ["PATH"]))
            v[name] = (r.stdout or r.stderr).strip().splitlines()[0][:80] if (r.stdout or r.stderr) else "?"
        except Exception as e:  # noqa: BLE001
            v[name] = f"unavailable: {e}"
    return v


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--covs", default="5,15,40")
    ap.add_argument("--tools", default="megahit,metaspades,penguin,penguin_guided,haploflow")
    ap.add_argument("--control-tools", default="megahit")
    ap.add_argument("--rungs", default="", help="comma list of rung tags; empty = all")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--work", type=Path, default=Path("/home/user/ladderwork"))
    ap.add_argument("--keep", type=Path, default=Path("/mnt/project-files/strain_ladder/contigs"))
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "pair_ladder")
    ap.add_argument("--rows-name", default="rows.jsonl")
    ap.add_argument("--background", default="",
                    help="r1.fq,r2.fq of a fixed real background added to every run (uncompressed)")
    args = ap.parse_args()

    genomes = {}
    for g in GENOMES:
        genomes |= read_fasta(g)
    ladder = json.loads((ROOT / "results/pair_ladder/ladder.json").read_text())
    want = set(filter(None, args.rungs.split(",")))
    jobs = []
    for rung in ladder:
        if want and rung["tag"] not in want:
            continue
        for pair in rung["pairs"]:
            for cov in [float(c) for c in args.covs.split(",")]:
                for seed in range(args.seed, args.seed + args.n_seeds):
                    for tool in args.tools.split(","):
                        jobs.append({"pair": pair, "mode": "mix", "cov": cov, "seed": seed, "tool": tool})
                    for tool in filter(None, args.control_tools.split(",")):
                        for mode in ("a_alone", "b_alone"):
                            jobs.append({"pair": pair, "mode": mode, "cov": cov, "seed": seed, "tool": tool})
    args.out.mkdir(parents=True, exist_ok=True)
    rows_path = args.out / args.rows_name
    done = set()
    if rows_path.exists():
        for line in rows_path.open():
            r = json.loads(line)
            done.add((r["rung"], r["a"], r["b"], r["mode"], r["cov_per_lineage"], r["seed"], r["tool"]))
    todo = [j for j in jobs if (j["pair"]["tag"], j["pair"]["a"], j["pair"]["b"], j["mode"],
                                 j["cov"], j["seed"], j["tool"]) not in done]
    print(f"{len(jobs)} jobs, {len(todo)} to run", flush=True)
    (args.out / f"manifest_{args.rows_name.split('.')[0]}.json").write_text(json.dumps({
        "script": "experiments/pair_ladder.py", "argv": sys.argv[1:],
        "git_sha": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                                  text=True).stdout.strip(),
        "read_len": READ_LEN, "insert": INSERT, "insert_sd": INSERT_SD, "art_profile": ART_PROFILE,
        "min_contig": MIN_CONTIG, "window": WINDOW, "versions": versions(),
        "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "jobs": len(jobs),
    }, indent=2))
    bg = tuple(Path(x) for x in args.background.split(",")) if args.background else None
    with ProcessPoolExecutor(args.jobs) as ex, rows_path.open("a") as fh:
        futs = [ex.submit(one_job, j, genomes, args.work, args.keep, args.threads, bg) for j in todo]
        for n, f in enumerate(as_completed(futs), 1):
            try:
                r = f.result()
            except Exception as e:  # noqa: BLE001
                print("job failed:", repr(e)[:300], flush=True)
                continue
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            if n % 10 == 0 or n == len(futs):
                print(f"{n}/{len(futs)} {r['tool']} {r['rung']} {r['mode']} c{r['cov_per_lineage']} "
                      f"{r['seconds']}s", flush=True)


if __name__ == "__main__":
    main()
