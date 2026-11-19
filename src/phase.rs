//! Splitting a merged contig into co-existing populations.
//!
//! Why this exists. On a metagenome holding two swarms of the same virus from
//! two different hosts, hapax returns both as one contig
//! (`results/swarm_metagenome/`). The joins responsible are not mistakes any
//! per-pair score could refuse: over a 150-nt window two hosts' sequences
//! often differ at one or two sites, which is inside sequencing error, so the
//! overlap really is a nucleotide-identical overlap.
//! `experiments/sister_pairs.py` shows the merge happens at 30x as readily as
//! at 300x and that no setting of the grouping thresholds prevents it. The
//! information that separates the hosts is not in any pair of reads; it is in
//! the population structure, and the place to use it is after consensus, where
//! the variant sites and the reads covering them are both in hand.
//!
//! What was tried first, and why it failed. A sister population should show a
//! chain of neighbouring variant sites in near-perfect linkage running the
//! length of the contig. It does not: on the merged HIV *env* contig only 97 of
//! 677 neighbouring site pairs reach r^2 >= 0.9 and the longest unbroken run
//! is 5 sites. Each host is *itself* a swarm with 60 haplotypes at ~5%
//! diversity, so between-host sites are interleaved with far more within-host
//! sites, and any chain through consecutive sites breaks immediately. Linkage
//! between *neighbours* is the wrong statistic.
//!
//! What works is to fit the two populations to all the sites at once: a
//! two-component Bernoulli mixture over per-read alleles, by EM. A site that
//! separates the hosts gets p0 near 0 and p1 near 1 (or the reverse); a site
//! that varies within a host gets p0 ~ p1 and contributes nothing either way,
//! which is exactly the asymmetry that makes a swarm different from two
//! swarms. The contig is split only when the fitted model has many sites
//! discriminating across most of its length, so a single swarm is left whole:
//! splitting a real population would be worse than leaving two merged, and the
//! chimaera number reports either outcome.

use crate::consensus::{Contig, Variant};
use crate::layout::Placement;

#[derive(Clone, Debug)]
pub struct PhaseParams {
    /// minor-allele frequency a site needs before it can carry phase
    pub min_af: f32,
    /// per-read sites needed before a read is used to fit, or assigned
    pub min_sites_per_read: usize,
    /// |p0 - p1| for a site to count as discriminating the two populations
    pub min_sep: f32,
    /// discriminating sites the fit needs
    pub min_sites: usize,
    /// fraction of the contig the discriminating sites must span
    pub span_frac: f32,
    /// reads each side needs, as a fraction of the contig's assigned reads
    pub min_side_frac: f32,
    /// how many times a contig may be split (2 = up to four populations)
    pub max_depth: usize,
    /// EM restarts, and iterations per restart
    pub restarts: usize,
    pub iters: usize,
    /// rounds of sequence-level reassignment after the mixture fit
    pub refine_rounds: usize,
}

impl Default for PhaseParams {
    fn default() -> Self {
        // min_af 0.10: a between-host site sits near the hosts' abundance
        // ratio, so this is the lowest-abundance population that can be split
        // out. min_sep 0.8 with min_sites 20 over half the contig is a
        // deliberately high bar: see the module note on which error is worse.
        PhaseParams { min_af: 0.10, min_sites_per_read: 3, min_sep: 0.8, min_sites: 20,
                      span_frac: 0.5, min_side_frac: 0.15, max_depth: 2, restarts: 6, iters: 60,
                      refine_rounds: 5 }
    }
}

/// Deterministic generator: a split must be reproducible from the reads alone.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// Per-read (site index, allele) pairs, consensus = 0 and alt = 1. Bases that
/// are neither are dropped: at a biallelic site a third base is an error.
fn observations(c: &Contig, sites: &[&Variant]) -> Vec<Vec<(u32, u8)>> {
    let mut by_read: Vec<Vec<(u32, u8)>> = vec![Vec::new(); c.alns.len()];
    for (si, v) in sites.iter().enumerate() {
        for (ri, a) in c.alns.iter().enumerate() {
            if v.pos < a.ref_start || v.pos >= a.ref_end {
                continue;
            }
            let b = a.cols[v.pos - a.ref_start];
            if b == v.cons {
                by_read[ri].push((si as u32, 0));
            } else if b == v.alt {
                by_read[ri].push((si as u32, 1));
            }
        }
    }
    by_read
}

struct Fit {
    ll: f64,
    p: Vec<[f64; 2]>,
    post: Vec<f64>,
}

fn em(obs: &[Vec<(u32, u8)>], n_sites: usize, used: &[usize], p: &PhaseParams, rng: &mut Lcg) -> Fit {
    let mut post: Vec<f64> = obs.iter().map(|_| 0.5).collect();
    for &r in used {
        post[r] = rng.next_f64();
    }
    let mut prob = vec![[0.5f64; 2]; n_sites];
    let mut pi = 0.5f64;
    let mut ll = f64::NEG_INFINITY;
    for _ in 0..p.iters {
        // M step: allele frequency per site per component, with a Beta(1,1) prior
        let mut num = vec![[1.0f64; 2]; n_sites];
        let mut den = vec![[2.0f64; 2]; n_sites];
        let mut wsum = 0.0;
        for &r in used {
            let w1 = post[r];
            let w0 = 1.0 - w1;
            wsum += w1;
            for &(si, a) in &obs[r] {
                let s = si as usize;
                den[s][0] += w0;
                den[s][1] += w1;
                if a == 1 {
                    num[s][0] += w0;
                    num[s][1] += w1;
                }
            }
        }
        for s in 0..n_sites {
            prob[s] = [num[s][0] / den[s][0], num[s][1] / den[s][1]];
        }
        pi = (wsum / used.len() as f64).clamp(0.05, 0.95);
        // E step
        let mut new_ll = 0.0;
        for &r in used {
            let mut l0 = (1.0 - pi).ln();
            let mut l1 = pi.ln();
            for &(si, a) in &obs[r] {
                let s = si as usize;
                let (q0, q1) = (prob[s][0], prob[s][1]);
                l0 += if a == 1 { q0.ln() } else { (1.0 - q0).ln() };
                l1 += if a == 1 { q1.ln() } else { (1.0 - q1).ln() };
            }
            let m = l0.max(l1);
            let z = (l0 - m).exp() + (l1 - m).exp();
            new_ll += m + z.ln();
            post[r] = (l1 - m).exp() / z;
        }
        if (new_ll - ll).abs() < 1e-6 {
            ll = new_ll;
            break;
        }
        ll = new_ll;
    }
    Fit { ll, p: prob.iter().map(|x| [x[0], x[1]]).collect(), post }
}

pub struct Split {
    /// read index (into the contig's `alns`) -> side, or None for "both sides"
    pub side: Vec<Option<bool>>,
    pub sites: usize,
    pub span: f32,
}

/// Fit two populations to the contig's variant sites; return the split if the
/// fit is strong and spread out enough to be population structure.
pub fn find_split(c: &Contig, p: &PhaseParams) -> Option<Split> {
    if c.seq.is_empty() || c.alns.len() < 20 {
        return None;
    }
    let sites: Vec<&Variant> = c.variants.iter().filter(|v| v.alt_freq >= p.min_af).collect();
    if sites.len() < p.min_sites {
        return None;
    }
    let obs = observations(c, &sites);
    let used: Vec<usize> = (0..obs.len()).filter(|&r| obs[r].len() >= p.min_sites_per_read).collect();
    if used.len() < 20 {
        return None;
    }
    let mut rng = Lcg(0x5EED ^ (c.seq.len() as u64) << 16 ^ sites.len() as u64);
    let mut best: Option<Fit> = None;
    for _ in 0..p.restarts {
        let f = em(&obs, sites.len(), &used, p, &mut rng);
        if best.as_ref().map_or(true, |b| f.ll > b.ll) {
            best = Some(f);
        }
    }
    let f = best?;
    let disc: Vec<usize> = (0..sites.len()).filter(|&s| (f.p[s][0] - f.p[s][1]).abs() >= p.min_sep as f64).collect();
    if disc.len() < p.min_sites {
        return None;
    }
    let span = (sites[*disc.last().unwrap()].pos - sites[disc[0]].pos) as f32 / c.seq.len() as f32;
    if span < p.span_frac {
        return None;
    }
    // Confident reads first: the mixture posterior decides them.
    let mut side = vec![None; c.alns.len()];
    let (mut n0, mut n1) = (0usize, 0usize);
    for r in 0..obs.len() {
        if obs[r].len() < p.min_sites_per_read {
            continue;
        }
        if f.post[r] >= 0.9 {
            side[r] = Some(true);
            n1 += 1;
        } else if f.post[r] <= 0.1 {
            side[r] = Some(false);
            n0 += 1;
        }
    }
    let need = ((p.min_side_frac * (n0 + n1) as f32) as usize).max(10);
    if n0.min(n1) < need {
        return None;
    }
    // The mixture posterior is a start, not an answer: on the real HIV *env*
    // pair it put each host at about 70/30 in its own side, and a 70/30 side
    // gives a mosaic consensus. So sharpen it by sequence. Each side's column
    // consensus is computed from the reads currently on it (in the parent
    // contig's coordinates, which both sides share), every read is then moved
    // to the side it matches at the columns where the two sides differ, and
    // that repeats. Once the two consensuses are distinct, a 150-nt read
    // carries about twenty differences between two hosts, so assignment by
    // sequence is far sharper than assignment by a handful of phased sites.
    // A read matching both equally (conserved sequence, where the populations
    // agree) stays on both: both populations contain it, and withholding it
    // would truncate one consensus.
    let n = c.seq.len();
    for _ in 0..p.refine_rounds {
        let mut cols = [vec![[0u32; 5]; n], vec![[0u32; 5]; n]];
        for (ri, a) in c.alns.iter().enumerate() {
            let s = match side[ri] {
                Some(s) => s as usize,
                None => continue,
            };
            for (k, &b) in a.cols.iter().enumerate() {
                let pos = a.ref_start + k;
                if pos < n {
                    cols[s][pos][b.min(4) as usize] += 1;
                }
            }
        }
        let maj: Vec<[u8; 2]> = (0..n)
            .map(|i| {
                [0usize, 1].map(|s| {
                    let cc = &cols[s][i];
                    let mut best = 255u8;
                    let mut bn = 0u32;
                    for b in 0..5u8 {
                        if cc[b as usize] > bn {
                            bn = cc[b as usize];
                            best = b;
                        }
                    }
                    // a column one side barely covers cannot discriminate
                    if bn < 3 { 255 } else { best }
                })
            })
            .collect();
        let mut moved = 0usize;
        for (ri, a) in c.alns.iter().enumerate() {
            let (mut m0, mut m1) = (0i32, 0i32);
            for (k, &b) in a.cols.iter().enumerate() {
                let pos = a.ref_start + k;
                if pos >= n {
                    break;
                }
                let (e0, e1) = (maj[pos][0], maj[pos][1]);
                if e0 == 255 || e1 == 255 || e0 == e1 {
                    continue;
                }
                if b == e0 {
                    m0 += 1;
                } else if b == e1 {
                    m1 += 1;
                }
            }
            let new_side = if m0 > m1 { Some(false) } else if m1 > m0 { Some(true) } else { None };
            if new_side != side[ri] {
                moved += 1;
            }
            side[ri] = new_side;
        }
        if moved == 0 {
            break;
        }
    }
    let n1 = side.iter().filter(|s| **s == Some(true)).count();
    let n0 = side.iter().filter(|s| **s == Some(false)).count();
    let need = ((p.min_side_frac * (n0 + n1) as f32) as usize).max(10);
    if n0.min(n1) < need {
        return None;
    }
    Some(Split { side, sites: disc.len(), span })
}

/// The read sets a split implies, as placements in the contig's coordinates.
/// Reads assigned to neither side appear in both.
pub fn split_placements(c: &Contig, s: &Split) -> (Vec<Placement>, Vec<Placement>) {
    let mut a = Vec::new();
    let mut b = Vec::new();
    for (ri, aln) in c.alns.iter().enumerate() {
        let pl = Placement { read: aln.read, rev: aln.rev, x: aln.ref_start as i64 };
        match s.side[ri] {
            Some(false) => a.push(pl),
            Some(true) => b.push(pl),
            None => {
                a.push(pl);
                b.push(pl);
            }
        }
    }
    (a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two populations 10% apart, each with its own within-population variation,
    /// must split; one population with the same amount of variation must not.
    fn synthetic(two: bool) -> (Vec<Vec<u8>>, Vec<Placement>) {
        let n = 2000usize;
        let mut rng = Lcg(7);
        let root: Vec<u8> = (0..n).map(|_| (rng.next_f64() * 4.0) as u8).collect();
        let mut pops = vec![root.clone()];
        if two {
            let mut b = root.clone();
            for i in 0..n {
                if rng.next_f64() < 0.10 {
                    b[i] = (b[i] + 1) % 4;
                }
            }
            pops.push(b);
        }
        let mut reads = Vec::new();
        let mut pl = Vec::new();
        let mut x = 0usize;
        for i in 0..1600 {
            let base = &pops[i % pops.len()];
            let start = ((rng.next_f64() * (n - 150) as f64) as usize).min(n - 150);
            let mut r: Vec<u8> = base[start..start + 150].to_vec();
            // within-population variation: a few sites per read, not shared
            for k in 0..150 {
                if rng.next_f64() < 0.01 {
                    r[k] = (r[k] + 1) % 4;
                }
            }
            reads.push(r);
            pl.push(Placement { read: i as u32, rev: false, x: start as i64 });
            x = x.max(start);
        }
        let _ = x;
        (reads, pl)
    }

    #[test]
    fn splits_two_populations_and_keeps_one_whole() {
        let p = PhaseParams::default();
        let cp = crate::consensus::ConsensusParams::default();
        let (reads, pl) = synthetic(true);
        let c = crate::consensus::build(&reads, &pl, &cp);
        let s = find_split(&c, &p).expect("two populations 10% apart should split");
        assert!(s.sites >= p.min_sites);
        let (a, b) = split_placements(&c, &s);
        assert!(a.len() > 100 && b.len() > 100);

        let (reads1, pl1) = synthetic(false);
        let c1 = crate::consensus::build(&reads1, &pl1, &cp);
        assert!(find_split(&c1, &p).is_none(), "one population must stay whole");
    }
}
