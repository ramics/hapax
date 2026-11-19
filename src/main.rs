//! hapax command line.
//!
//! hapax assemble -i reads.fa -o outdir [options]
//! hapax orfs -i genomes.fa [--orf-caller scored|union] [--min-orf INT]
//!     prints the ORFs and per-base region runs the annotation would use
//!
//! Options (defaults in brackets):
//!   --mode codon|nucleotide    overlap scoring [codon]; nucleotide is the ablation baseline
//!   --seed contiguous|spaced   seed shape [contiguous]; spaced = codon-wobble-masked
//!   -k INT  -w INT             seed informative positions / minimizer window [15 / 5]
//!   --tau FLOAT                tier-2 accept: min log-odds per aligned column [0.8]
//!   --min-ident FLOAT          tier-2 accept: min nucleotide identity [0.70]
//!   --easy-div / --hard-div    tier-1 accept / reject ungapped divergence [0.03 / 0.30]
//!   --div FLOAT                model pairwise divergence [0.10]
//!   --omega FLOAT              model dN/dS [0.23]
//!   --sub-err FLOAT            per-base substitution error [0.005]
//!   --max-partners INT         candidate partners per read [16]
//!   --bucket-cap INT           per-bucket, per-strand entries paired with each read [48]; 0 = all
//!   --orf-caller scored|union  region annotation ORF caller [scored]
//!   --min-orf INT              shortest ORF annotated [300]
//!   --min-contig INT           shortest contig written [300]
//!   --best-k INT               stage-1 clusters from each read's top k edges [6]
//!   --lazy-t2 0|1              skip the HMM for pairs tier 1 already joined [1]
//!   --stage1 bestk|tier1       first-stage clusters for grouping [bestk]
//!   --min-bridge INT           edges needed to merge two clusters [5]
//!   --bridge-frac FLOAT        ...and as a fraction of the smaller cluster [0.1]
//!   --dump-edges 0|1           write edges.tsv for calibration [0]
//!   --split-cut FLOAT          split a read group at a cut whose crossing overlaps
//!                              are below this fraction of the smaller side's own
//!                              overlaps (0 = off) [0]; two swarms of one virus are
//!                              two dense clusters joined by a thin cut
//!   --split-cut-depth INT      how many times a group may be cut [1]
//!   --split-phase 0|1          split a contig into two populations when a
//!                              two-component fit over its variant sites separates
//!                              them: two swarms of the same virus join through
//!                              conserved windows and no per-pair score can refuse
//!                              those joins (src/phase.rs) [0]
//!   --split-sep FLOAT          |p0-p1| for a site to discriminate two populations [0.8]
//!   --split-min-af FLOAT       minor-allele frequency a site needs to carry phase [0.10]
//!   --split-sites INT          sites the chain needs [20]
//!   --split-span FLOAT         fraction of the contig the chain must span [0.5]
//!   --split-depth INT          how many times a contig may be split [2]
//!   -t INT                     threads [all]

use hapax::consensus::{self, ConsensusParams, Region};
use hapax::index::{Index, SeedParams, SeedShape};
use hapax::layout;
use hapax::model::{ModelParams, Scorer};
use hapax::overlap::{self, DispatchStats, Mode, OverlapParams};
use hapax::phase::{self, PhaseParams};
use hapax::lineage::{self, LineageParams};
use hapax::seq;
use std::io::{BufWriter, Write};
use std::time::Instant;

/// A contig whose variants form one linked chain across it holds two
/// populations that overlap scoring cannot separate (see src/phase.rs); rebuild
/// each side as its own contig, and try again on each side up to `max_depth`.
fn split_recursive(
    c: consensus::Contig,
    reads: &[Vec<u8>],
    cp: &ConsensusParams,
    pp: &PhaseParams,
    depth: usize,
    counter: &std::sync::atomic::AtomicUsize,
) -> Vec<consensus::Contig> {
    if depth >= pp.max_depth {
        return vec![c];
    }
    match phase::find_split(&c, pp) {
        None => vec![c],
        Some(s) => {
            let (a, b) = phase::split_placements(&c, &s);
            if a.len() < 2 || b.len() < 2 {
                return vec![c];
            }
            counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let ca = consensus::build(reads, &a, cp);
            let cb = consensus::build(reads, &b, cp);
            let mut out = split_recursive(ca, reads, cp, pp, depth + 1, counter);
            out.extend(split_recursive(cb, reads, cp, pp, depth + 1, counter));
            out
        }
    }
}

/// Split a contig into lineages (src/lineage.rs) and rebuild each; repeat on each
/// piece, since a rebuilt consensus is cleaner and can expose a lineage the
/// mixed one hid.
fn split_lineages(
    c: consensus::Contig,
    reads: &[Vec<u8>],
    cp: &ConsensusParams,
    lp: &LineageParams,
    depth: usize,
    counter: &std::sync::atomic::AtomicUsize,
) -> Vec<consensus::Contig> {
    if depth >= lp.max_depth {
        return vec![c];
    }
    let sets = lineage::split(&c, lp);
    if sets.is_empty() {
        return vec![c];
    }
    counter.fetch_add(sets.len() - 1, std::sync::atomic::Ordering::Relaxed);
    sets.iter().flat_map(|s| split_lineages(consensus::build(reads, s, cp), reads, cp, lp, depth + 1, counter)).collect()
}

fn usage() -> ! {
    eprintln!("usage: hapax assemble -i reads.fa -o outdir [options]  (see src/main.rs header)");
    std::process::exit(2);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() >= 2 && args[1] == "orfs" {
        return orfs_cmd(&args[2..]);
    }
    if args.len() < 2 || args[1] != "assemble" {
        usage();
    }
    let mut input = None;
    let mut out = None;
    let mut sp = SeedParams::default();
    let mut op = OverlapParams::default();
    let mut mp = ModelParams::default();
    let mut cp = ConsensusParams::default();
    let mut min_contig = 300usize;
    let mut threads = 0usize;
    let mut dump_edges = false;
    let mut split_phase = false;
    let mut split_lineage = false;
    let mut lp = LineageParams::default();
    let mut split_cut = 0.0f32;
    let mut split_cut_depth = 1usize;
    let mut pp = PhaseParams::default();
    let mut best_k = 6usize;
    let mut min_bridge = 5usize;
    let mut stage1_tier1 = false;
    let mut bridge_frac = 0.1f32;
    let mut i = 2;
    while i < args.len() {
        let a = args[i].as_str();
        let v = || args.get(i + 1).cloned().unwrap_or_else(|| usage());
        match a {
            "-i" => input = Some(v()),
            "-o" => out = Some(v()),
            "--mode" => op.mode = if v() == "nucleotide" { Mode::Nucleotide } else { Mode::Codon },
            "--seed" => sp.shape = if v() == "spaced" { SeedShape::CodonSpaced } else { SeedShape::Contiguous },
            "-k" => sp.k = v().parse().unwrap(),
            "-w" => sp.w = v().parse().unwrap(),
            "--tau" => op.tau = v().parse().unwrap(),
            "--min-ident" => op.min_ident = v().parse().unwrap(),
            "--easy-div" => op.easy_div = v().parse().unwrap(),
            "--hard-div" => op.hard_div = v().parse().unwrap(),
            "--div" => mp.divergence = v().parse().unwrap(),
            "--omega" => mp.omega = v().parse().unwrap(),
            "--sub-err" => mp.sub_err = v().parse().unwrap(),
            "--max-partners" => op.max_partners = v().parse().unwrap(),
            "--bucket-cap" => op.bucket_cap = v().parse().unwrap(),
            "--orf-caller" => cp.orf_caller = parse_caller(&v()),
            "--min-orf" => cp.min_orf = v().parse().unwrap(),
            "--rounds" => cp.rounds = v().parse().unwrap(),
            "--min-contig" => min_contig = v().parse().unwrap(),
            "-t" => threads = v().parse().unwrap(),
            "--best-k" => best_k = v().parse().unwrap(),
            "--lazy-t2" => op.lazy_t2 = v() == "1",
            "--stage1" => stage1_tier1 = v() == "tier1",
            "--min-bridge" => min_bridge = v().parse().unwrap(),
            "--bridge-frac" => bridge_frac = v().parse().unwrap(),
            "--dump-edges" => dump_edges = v() == "1",
            "--split-cut" => split_cut = v().parse().unwrap(),
            "--split-cut-depth" => split_cut_depth = v().parse().unwrap(),
            "--split-phase" => split_phase = v() == "1",
            "--split-lineage" => split_lineage = v() == "1",
            "--lineage-div" => lp.max_div = v().parse().unwrap(),
            "--lineage-min-reads" => lp.min_lineage_reads = v().parse().unwrap(),
            "--lineage-depth" => lp.max_depth = v().parse().unwrap(),
            "--lineage-split-div" => lp.min_split_div = v().parse().unwrap(),
            "--lineage-conflict-frac" => lp.max_conflict_frac = v().parse().unwrap(),
            "--split-sep" => pp.min_sep = v().parse().unwrap(),
            "--split-min-af" => pp.min_af = v().parse().unwrap(),
            "--split-sites" => pp.min_sites = v().parse().unwrap(),
            "--split-span" => pp.span_frac = v().parse().unwrap(),
            "--split-depth" => pp.max_depth = v().parse().unwrap(),
            _ => usage(),
        }
        i += 2;
    }
    let (input, out) = (input.unwrap_or_else(|| usage()), out.unwrap_or_else(|| usage()));
    if threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    }
    std::fs::create_dir_all(&out).unwrap();
    let t0 = Instant::now();
    let recs = seq::read_seqs(&input).expect("reading input");
    let names: Vec<String> = recs.iter().map(|r| r.name.clone()).collect();
    let reads: Vec<Vec<u8>> = recs.into_iter().map(|r| r.seq).collect();
    let t_read = t0.elapsed().as_secs_f64();

    let t = Instant::now();
    let idx = Index::build(&reads, sp);
    let t_index = t.elapsed().as_secs_f64();

    let t = Instant::now();
    let cands = overlap::candidates(&reads, &idx, &op);
    let t_cand = t.elapsed().as_secs_f64();

    let t = Instant::now();
    let scorer = Scorer::new(mp.clone());
    let st = DispatchStats::default();
    let disp = overlap::dispatch(&reads, &cands, &scorer, &op, &st);
    let edges = disp.edges;
    let t_overlap = t.elapsed().as_secs_f64();

    if dump_edges {
        let mut f = BufWriter::new(std::fs::File::create(format!("{}/edges.tsv", out)).unwrap());
        writeln!(f, "a\tb\trev\toffset\ttier\tscore\tper_col\tidentity").unwrap();
        for e in &edges {
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.4}\t{:.4}", names[e.a as usize], names[e.b as usize], e.rev as u8, e.offset, e.tier, e.score, e.per_col, e.identity).unwrap();
        }
    }
    let t = Instant::now();
    let lens: Vec<usize> = reads.iter().map(|r| r.len()).collect();
    let stage1 = if stage1_tier1 { Some(&disp.tier1_comp[..]) } else { None };
    let mut groups = layout::group(reads.len(), &edges, best_k, min_bridge, bridge_frac, stage1);
    if split_cut > 0.0 {
        groups = layout::split_sparse(reads.len(), &groups, &edges, split_cut, 0.15, split_cut_depth, 50);
    }
    let comps = layout::layout(&lens, &edges, Some(&groups));
    let n_groups = { let mut g = groups.clone(); g.sort(); g.dedup(); g.len() };
    let t_layout = t.elapsed().as_secs_f64();

    let t = Instant::now();
    use rayon::prelude::*;
    let splits = std::sync::atomic::AtomicUsize::new(0);
    let contigs: Vec<consensus::Contig> = comps
        .par_iter()
        .filter(|c| c.len() >= 2)
        .flat_map(|c| {
            let built = consensus::build(&reads, c, &cp);
            if split_lineage {
                return split_lineages(built, &reads, &cp, &lp, 0, &splits);
            }
            if split_phase {
                split_recursive(built, &reads, &cp, &pp, 0, &splits)
            } else {
                vec![built]
            }
        })
        .filter(|c| c.seq.len() >= min_contig)
        .collect();
    let t_cons = t.elapsed().as_secs_f64();
    let n_splits = splits.load(std::sync::atomic::Ordering::Relaxed);

    // ---- outputs
    let mut fa = BufWriter::new(std::fs::File::create(format!("{}/contigs.fa", out)).unwrap());
    let mut lay = BufWriter::new(std::fs::File::create(format!("{}/layout.tsv", out)).unwrap());
    let mut var = BufWriter::new(std::fs::File::create(format!("{}/variants.tsv", out)).unwrap());
    let mut lnk = BufWriter::new(std::fs::File::create(format!("{}/linkage.tsv", out)).unwrap());
    let mut reg = BufWriter::new(std::fs::File::create(format!("{}/regions.tsv", out)).unwrap());
    writeln!(lay, "contig\tread\tstart\tend\tstrand").unwrap();
    writeln!(var, "contig\tpos\tdepth\tcons\talt\talt_freq\tA\tC\tG\tT\tdel\tregion\teffects").unwrap();
    writeln!(lnk, "contig\tpos1\tpos2\tr2\thaplotypes").unwrap();
    writeln!(reg, "contig\tstart\tend\tregion\torf_frames").unwrap();
    let bases = |b: u8| (b"ACGT-"[b.min(4) as usize]) as char;
    let mut n_var = 0;
    let mut total_len = 0;
    for (ci, c) in contigs.iter().enumerate() {
        let name = format!("contig_{}", ci + 1);
        total_len += c.seq.len();
        seq::write_fasta(&mut fa, &format!("{} len={} reads={} variants={}", name, c.seq.len(), c.alns.len(), c.variants.len()), &c.seq).unwrap();
        for a in &c.alns {
            writeln!(lay, "{}\t{}\t{}\t{}\t{}", name, names[a.read as usize], a.ref_start, a.ref_end, if a.rev { '-' } else { '+' }).unwrap();
        }
        for v in &c.variants {
            n_var += 1;
            let eff: Vec<String> = v.effects.iter().map(|(f, s)| format!("f{}:{}", f, if *s { "syn" } else { "nonsyn" })).collect();
            writeln!(
                var,
                "{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                name,
                v.pos + 1,
                v.depth,
                bases(v.cons),
                bases(v.alt),
                v.alt_freq,
                v.counts[0],
                v.counts[1],
                v.counts[2],
                v.counts[3],
                v.counts[4],
                v.region.name(),
                if eff.is_empty() { "-".into() } else { eff.join(",") }
            )
            .unwrap();
        }
        for (s1, s2, h, r2) in consensus::linkage(c) {
            let hs: Vec<String> = h.iter().map(|((a, b), n)| format!("{}{}:{}", bases(*a), bases(*b), n)).collect();
            writeln!(lnk, "{}\t{}\t{}\t{:.3}\t{}", name, s1 + 1, s2 + 1, r2, hs.join(",")).unwrap();
        }
        // run-length regions
        let mut s = 0;
        while s < c.region.len() {
            let mut e = s;
            while e < c.region.len() && c.region[e] == c.region[s] {
                e += 1;
            }
            let frames: Vec<String> = c.orfs.iter().filter(|o| o.start < e && o.end > s).map(|o| o.frame.to_string()).collect();
            let frames = if c.region[s] == Region::Noncoding { "-".to_string() } else { frames.join(",") };
            writeln!(reg, "{}\t{}\t{}\t{}\t{}", name, s + 1, e, c.region[s].name(), frames).unwrap();
            s = e;
        }
    }
    let g = |x: &std::sync::atomic::AtomicU64| st.get(x);
    let placed: usize = contigs.iter().map(|c| c.alns.len()).sum();
    let stats = format!(
        concat!(
            "{{\n  \"reads\": {},\n  \"mode\": \"{:?}\",\n  \"seed_shape\": \"{:?}\",\n  \"k\": {},\n  \"w\": {},\n  \"tau\": {},\n",
            "  \"candidates\": {},\n  \"tier1_accept\": {},\n  \"tier1_reject\": {},\n  \"tier2\": {},\n  \"tier2_accept\": {},\n",
            "  \"tier2_dual_fallback\": {},\n  \"tier2_skipped_same_cluster\": {},\n  \"hmm_fraction\": {:.4},\n  \"tier2_codon_scored_nt\": {},\n  \"tier2_nucleotide_scored_nt\": {},\n",
            "  \"tier2_syn_mismatch\": {},\n  \"tier2_nonsyn_mismatch\": {},\n  \"phase_wins\": [{}, {}, {}],\n  \"edges\": {},\n  \"read_groups\": {},\n",
            "  \"components\": {},\n  \"phase_splits\": {},\n  \"contigs\": {},\n  \"contig_bases\": {},\n  \"placed_reads\": {},\n  \"variants\": {},\n",
            "  \"seconds\": {{\"read\": {:.3}, \"index\": {:.3}, \"candidates\": {:.3}, \"overlap\": {:.3}, \"layout\": {:.3}, \"consensus\": {:.3}, \"total\": {:.3}}}\n}}\n"
        ),
        reads.len(),
        op.mode,
        sp.shape,
        sp.k,
        sp.w,
        op.tau,
        g(&st.candidates),
        g(&st.t1_accept),
        g(&st.t1_reject),
        g(&st.t2),
        g(&st.t2_accept),
        g(&st.t2_dual),
        g(&st.t2_skipped),
        st.hmm_fraction(),
        g(&st.t2_codon_nt),
        g(&st.t2_nuc_nt),
        g(&st.t2_syn),
        g(&st.t2_nonsyn),
        g(&st.phase_wins[0]),
        g(&st.phase_wins[1]),
        g(&st.phase_wins[2]),
        edges.len(),
        n_groups,
        comps.len(),
        n_splits,
        contigs.len(),
        total_len,
        placed,
        n_var,
        t_read,
        t_index,
        t_cand,
        t_overlap,
        t_layout,
        t_cons,
        t0.elapsed().as_secs_f64()
    );
    std::fs::write(format!("{}/stats.json", out), &stats).unwrap();
    eprint!("{}", stats);
}

fn parse_caller(s: &str) -> hapax::orf::OrfCaller {
    if s == "union" {
        hapax::orf::OrfCaller::Union
    } else if s == "atg" {
        hapax::orf::OrfCaller::Atg
    } else {
        hapax::orf::OrfCaller::Scored
    }
}

fn orfs_cmd(args: &[String]) {
    let mut input = None;
    let mut caller = hapax::orf::OrfCaller::Scored;
    let mut min_orf = ConsensusParams::default().min_orf;
    let mut i = 0;
    while i < args.len() {
        let v = args.get(i + 1).cloned().unwrap_or_else(|| usage());
        match args[i].as_str() {
            "-i" => input = Some(v),
            "--orf-caller" => caller = parse_caller(&v),
            "--min-orf" => min_orf = v.parse().unwrap(),
            _ => usage(),
        }
        i += 2;
    }
    let recs = seq::read_seqs(&input.unwrap_or_else(|| usage())).expect("reading input");
    println!("seq\tkind\tstart\tend\tvalue");
    for r in recs {
        let orfs = hapax::orf::call_orfs(&r.seq, min_orf, caller);
        for o in &orfs {
            println!("{}\torf\t{}\t{}\t{}", r.name, o.start + 1, o.end, o.frame);
        }
        let reg = consensus::annotate_regions(r.seq.len(), &orfs);
        let mut s = 0;
        while s < reg.len() {
            let mut e = s;
            while e < reg.len() && reg[e] == reg[s] {
                e += 1;
            }
            println!("{}\tregion\t{}\t{}\t{}", r.name, s + 1, e, reg[s].name());
            s = e;
        }
    }
}
