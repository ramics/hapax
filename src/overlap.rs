//! Candidate generation and two-tier overlap dispatch.
//!
//! Tier 1 is cheap nucleotide work: minimizer anchors are grouped by read pair
//! and diagonal, and an ungapped comparison along the best diagonal settles the
//! easy majority. A pair whose anchors all sit on one diagonal (no indel) and
//! whose ungapped mismatch rate is below `easy_div` is accepted outright -- a
//! nucleotide assembler would join it too, and the codon model has nothing to
//! add. A pair whose ungapped mismatch rate is above `hard_div` is rejected
//! outright: no sensible within-population join looks like that.
//!
//! Everything in between is ambiguous and goes to tier 2, the codon HMM
//! (`hmm::score_pair`). This is the zone where the synonymous/non-synonymous
//! character of the differences is what decides, and it is the only place the
//! expensive machinery runs. The fraction of candidates reaching tier 2 is the
//! kill-criterion number of METHOD.md 5 and is reported on every run.
//!
//! Joins always demand shared nucleotide identity: a pair is only a candidate
//! if it shares exact seeds on a consistent diagonal, and an accepted tier-2
//! overlap must also clear a nucleotide-identity floor (`min_ident`). The codon
//! model changes how differences are *weighed*, never whether homology is
//! demanded.

use crate::hmm::{self, HmmParams};
use crate::index::Index;
use crate::model::Scorer;
use crate::seq::revcomp;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// codon scoring with per-region fallback (the method)
    Codon,
    /// nucleotide scoring only (internal ablation baseline)
    Nucleotide,
}

#[derive(Clone, Debug)]
pub struct OverlapParams {
    pub mode: Mode,
    pub min_anchors: usize,
    pub min_overlap: usize,
    /// max candidate partners kept per read (by anchor support)
    pub max_partners: usize,
    /// diagonal tolerance when clustering anchors
    pub diag_tol: i64,
    pub easy_div: f32,
    pub hard_div: f32,
    pub min_ident: f32,
    /// minimum log-odds per aligned column for a tier-2 accept
    pub tau: f32,
    pub dual_margin: f32,
    pub hmm: HmmParams,
}

impl Default for OverlapParams {
    fn default() -> Self {
        OverlapParams {
            mode: Mode::Codon,
            min_anchors: 2,
            min_overlap: 50,
            max_partners: 16,
            diag_tol: 10,
            easy_div: 0.03,
            hard_div: 0.30,
            min_ident: 0.70,
            tau: 0.8,
            dual_margin: 8.0,
            hmm: HmmParams::default(),
        }
    }
}

/// A candidate pair from seeding: B (oriented by `rev`) is expected to start at
/// A coordinate `diag`.
#[derive(Clone, Copy, Debug)]
pub struct Candidate {
    pub a: u32,
    pub b: u32,
    pub rev: bool,
    pub diag: i64,
    pub anchors: u32,
    pub spread: i64,
}

#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub a: u32,
    pub b: u32,
    pub rev: bool,
    /// start of oriented B in A coordinates
    pub offset: i64,
    pub score: f32,
    pub per_col: f32,
    pub identity: f32,
    pub tier: u8,
}

#[derive(Default, Debug)]
pub struct DispatchStats {
    pub candidates: AtomicU64,
    pub t1_accept: AtomicU64,
    pub t1_reject: AtomicU64,
    pub t2: AtomicU64,
    pub t2_accept: AtomicU64,
    pub t2_dual: AtomicU64,
    pub t2_codon_nt: AtomicU64,
    pub t2_nuc_nt: AtomicU64,
    pub t2_syn: AtomicU64,
    pub t2_nonsyn: AtomicU64,
    pub phase_wins: [AtomicU64; 3],
}

impl DispatchStats {
    pub fn get(&self, x: &AtomicU64) -> u64 {
        x.load(Ordering::Relaxed)
    }
    pub fn hmm_fraction(&self) -> f64 {
        let c = self.get(&self.candidates);
        if c == 0 {
            0.0
        } else {
            self.get(&self.t2) as f64 / c as f64
        }
    }
}

/// Convert a candidate seen from read `a` into the frame of read `b`.
fn flip(c: Candidate, la: i64, lb: i64) -> Candidate {
    let diag = if c.rev { lb - la + c.diag } else { -c.diag };
    Candidate { a: c.b, b: c.a, rev: c.rev, diag, anchors: c.anchors, spread: c.spread }
}

pub fn candidates(reads: &[Vec<u8>], idx: &Index, p: &OverlapParams) -> Vec<Candidate> {
    let span = idx.params.span() as i64;
    let per_read: Vec<Vec<Candidate>> = (0..reads.len())
        .into_par_iter()
        .map(|i| {
            let la = reads[i].len() as i64;
            let mut anchors: Vec<(u32, bool, i64)> = Vec::new();
            for s in &idx.per_read[i] {
                for &(_, j, pb, rb) in idx.lookup(s.hash) {
                    if j as usize == i {
                        continue;
                    }
                    let lb = reads[j as usize].len() as i64;
                    let rev = rb != s.rev;
                    let d = if !rev { s.pos as i64 - pb as i64 } else { s.pos as i64 - (lb - pb as i64 - span) };
                    anchors.push((j, rev, d));
                }
            }
            anchors.sort_unstable();
            let mut out: Vec<Candidate> = Vec::new();
            let mut g = 0;
            while g < anchors.len() {
                let mut h = g;
                while h < anchors.len() && anchors[h].0 == anchors[g].0 && anchors[h].1 == anchors[g].1 {
                    h += 1;
                }
                // best diagonal cluster within the group (diags sorted)
                let grp = &anchors[g..h];
                let (mut best_n, mut best_lo, mut best_hi) = (0usize, 0usize, 0usize);
                let mut lo = 0;
                for hi in 0..grp.len() {
                    while grp[hi].2 - grp[lo].2 > p.diag_tol {
                        lo += 1;
                    }
                    if hi - lo + 1 > best_n {
                        best_n = hi - lo + 1;
                        best_lo = lo;
                        best_hi = hi;
                    }
                }
                if best_n >= p.min_anchors {
                    let cl = &grp[best_lo..=best_hi];
                    let med = cl[cl.len() / 2].2;
                    let lb = reads[grp[0].0 as usize].len() as i64;
                    let ov = (la.min(med + lb) - med.max(0)) as usize;
                    if ov >= p.min_overlap {
                        out.push(Candidate {
                            a: i as u32,
                            b: grp[0].0,
                            rev: grp[0].1,
                            diag: med,
                            anchors: best_n as u32,
                            spread: cl[cl.len() - 1].2 - cl[0].2,
                        });
                    }
                }
                g = h;
            }
            out.sort_unstable_by(|x, y| y.anchors.cmp(&x.anchors));
            out.truncate(p.max_partners);
            out
        })
        .collect();
    // canonicalise to a < b and dedupe
    let mut all: Vec<Candidate> = per_read
        .into_iter()
        .flatten()
        .map(|c| if c.a < c.b { c } else { flip(c, reads[c.a as usize].len() as i64, reads[c.b as usize].len() as i64) })
        .collect();
    all.par_sort_unstable_by_key(|c| (c.a, c.b, c.rev, std::cmp::Reverse(c.anchors)));
    all.dedup_by_key(|c| (c.a, c.b, c.rev));
    all
}

/// Run the two tiers over all candidates.
pub fn dispatch(reads: &[Vec<u8>], cands: &[Candidate], sc: &Scorer, p: &OverlapParams, st: &DispatchStats) -> Vec<Edge> {
    st.candidates.fetch_add(cands.len() as u64, Ordering::Relaxed);
    cands
        .par_iter()
        .filter_map(|c| {
            let a = &reads[c.a as usize];
            let braw = &reads[c.b as usize];
            let bo;
            let b: &[u8] = if c.rev {
                bo = revcomp(braw);
                &bo
            } else {
                braw
            };
            // ---- tier 1: ungapped along the median diagonal
            let d = c.diag;
            let s0 = d.max(0) as usize;
            let e0 = (a.len() as i64).min(d + b.len() as i64).max(0) as usize;
            if e0 <= s0 || e0 - s0 < p.min_overlap {
                st.t1_reject.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            let mut mm = 0usize;
            let mut sc1 = 0.0f32;
            for ia in s0..e0 {
                let x = a[ia];
                let y = b[(ia as i64 - d) as usize];
                if x != y {
                    mm += 1;
                }
                sc1 += sc.nuc(x, y);
            }
            let l = (e0 - s0) as f32;
            let div = mm as f32 / l;
            if c.spread == 0 && div <= p.easy_div {
                st.t1_accept.fetch_add(1, Ordering::Relaxed);
                return Some(Edge { a: c.a, b: c.b, rev: c.rev, offset: d, score: sc1, per_col: sc1 / l, identity: 1.0 - div, tier: 1 });
            }
            if c.spread == 0 && div > p.hard_div {
                st.t1_reject.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            // ---- tier 2: codon HMM over three offsets (or nucleotide-only)
            st.t2.fetch_add(1, Ordering::Relaxed);
            let aln = match p.mode {
                Mode::Nucleotide => hmm::align(a, b, d, None, sc, &p.hmm)?,
                Mode::Codon => {
                    let ps = hmm::score_pair(a, b, d, sc, &p.hmm, p.dual_margin)?;
                    if ps.dual {
                        st.t2_dual.fetch_add(1, Ordering::Relaxed);
                        ps.nucleotide
                    } else {
                        // count a phase win only where codon scoring was actually used
                        // and beat the other grids (identical reads tie on all three)
                        let second = ps.per_phase.iter().filter(|&&x| x < ps.best.score).fold(f32::MIN, |m, &x| m.max(x));
                        if ps.best.codon_nt > 0 && ps.per_phase.iter().filter(|&&x| x >= ps.best.score).count() == 1 && second < ps.best.score {
                            st.phase_wins[ps.best.phase as usize].fetch_add(1, Ordering::Relaxed);
                        }
                        ps.best
                    }
                }
            };
            st.t2_codon_nt.fetch_add(aln.codon_nt as u64, Ordering::Relaxed);
            st.t2_nuc_nt.fetch_add(aln.nuc_nt as u64, Ordering::Relaxed);
            st.t2_syn.fetch_add(aln.syn as u64, Ordering::Relaxed);
            st.t2_nonsyn.fetch_add(aln.nonsyn as u64, Ordering::Relaxed);
            let per_col = aln.score / aln.cols.max(1) as f32;
            if aln.identity() >= p.min_ident && aln.overlap_len() as usize >= p.min_overlap && per_col >= p.tau {
                st.t2_accept.fetch_add(1, Ordering::Relaxed);
                Some(Edge {
                    a: c.a,
                    b: c.b,
                    rev: c.rev,
                    offset: aln.a_start as i64 - aln.b_start as i64,
                    score: aln.score,
                    per_col,
                    identity: aln.identity(),
                    tier: 2,
                })
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn flip_roundtrip() {
        let c = Candidate { a: 1, b: 2, rev: true, diag: 40, anchors: 3, spread: 0 };
        let f = flip(flip(c, 150, 140), 140, 150);
        assert_eq!(f.diag, 40);
        let c = Candidate { rev: false, ..c };
        assert_eq!(flip(c, 150, 140).diag, -40);
    }
}
