//! ORF calling for region annotation (coding / dual-coding / noncoding).
//!
//! What this feeds: `regions.tsv` and the syn/nonsyn effect column of
//! `variants.tsv`. It does **not** feed overlap scoring. The overlap HMM picks
//! its scoring mode per region from the reads themselves (its nucleotide state
//! and the two-offsets-both-win dual fallback in `hmm::score_pair`), so a bad
//! ORF call mislabels output but cannot switch codon scoring off.
//!
//! The first caller (`Union`) took every stop-free stretch of at least
//! `min_len` in all six frames and called any base covered twice "dual". On
//! real astrovirus sequence (NC_001943) that labelled roughly a quarter of the
//! genome dual-coding where the true overlap is 8 nt (found by the wastewater
//! probe). The reason: real coding sequence is full of long chance stop-free
//! runs in its shadow frames, and a longer minimum length does not remove them.
//!
//! `Scored`, the default, does what gene finders do in miniature:
//! 1. **Candidates** start at the first ATG of a stop-free segment, or at the
//!    contig edge when the segment runs off it (contigs are often partial).
//! 2. **Training set**: the longest candidates, taken greedily with no base
//!    overlap. These are almost always real genes.
//! 3. **Codon-usage score**: each candidate's mean per-codon log ratio of the
//!    training set's in-frame codon usage against the codon usage the same
//!    training spans show in their five *other* frames. Real genes read in
//!    their own frame look like the genome's genes; shadow ORFs inherit the
//!    shifted statistics of the frame they shadow, and score below zero.
//! 4. A candidate outside the training set -- one that overlaps a longer ORF
//!    in another frame -- is kept only if it scores at least `min_llr` and is
//!    at least 1.5 x `min_len` long. Overlapping genes need more evidence than
//!    free-standing ones, because the overlap is exactly where chance shadow
//!    runs live. Genuine nested genes of that size (HBV S inside P) survive;
//!    short ones (HIV rev exon 2, vpr) are below what a de novo caller can
//!    separate from chance and are not called.
//!
//! Known remaining error: an ORF's start is its first in-frame ATG, which can
//! sit upstream of the true start and create a false tail-to-head overlap
//! (astrovirus ORF2 starts ~180 nt early, overlapping ORF1b).
//!
//! Validated against GenBank CDS annotations of HIV-1 HXB2, HBV and human
//! astrovirus by `experiments/orf_annotation.py`.

use crate::genetic::{codon_index, is_stop};
use crate::seq::revcomp;

#[derive(Clone, Debug)]
pub struct Orf {
    /// 0..3 forward frames, 3..6 reverse
    pub frame: u8,
    pub start: usize,
    pub end: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrfCaller {
    Union,
    Atg,
    Scored,
}

/// Minimum mean per-codon log-likelihood ratio for a non-training ORF.
/// Set with `experiments/orf_annotation.py` on three genomes (see results);
/// a shadow ORF scores well below this, a real overlapping gene near or above.
pub const DEFAULT_MIN_LLR: f64 = -0.1;

fn to_fwd(strand: usize, f: usize, start: usize, end: usize, len: usize) -> Orf {
    if strand == 0 {
        Orf { frame: f as u8, start, end }
    } else {
        Orf { frame: 3 + f as u8, start: len - end, end: len - start }
    }
}

/// Stop-free segments per frame, in strand-local coordinates:
/// (strand, frame, seg_start, seg_end, open_left) where open_left means the
/// segment runs off the start of the sequence.
fn segments(seq: &[u8]) -> Vec<(usize, usize, usize, usize, bool)> {
    let rc = revcomp(seq);
    let mut out = Vec::new();
    for strand in 0..2 {
        let s = if strand == 0 { seq } else { &rc[..] };
        for f in 0..3 {
            let mut start = f;
            let mut open = true;
            let mut i = f;
            while i + 3 <= s.len() {
                if codon_index(&s[i..i + 3]).map_or(false, is_stop) {
                    out.push((strand, f, start, i, open));
                    start = i + 3;
                    open = false;
                }
                i += 3;
            }
            let end = f + (s.len() - f) / 3 * 3;
            if end > start {
                out.push((strand, f, start, end, open));
            }
        }
    }
    out
}

/// The original caller: every stop-free stretch >= min_len in six frames.
pub fn union_orfs(seq: &[u8], min_len: usize) -> Vec<Orf> {
    segments(seq)
        .into_iter()
        .filter(|&(_, _, s, e, _)| e - s >= min_len)
        .map(|(strand, f, s, e, _)| to_fwd(strand, f, s, e, seq.len()))
        .collect()
}

fn candidates(seq: &[u8], min_len: usize) -> Vec<Orf> {
    let rc = revcomp(seq);
    let mut out = Vec::new();
    for (strand, f, s0, e, open) in segments(seq) {
        let s = if strand == 0 { seq } else { &rc[..] };
        let start = if open {
            Some(s0)
        } else {
            (s0..e).step_by(3).find(|&i| s[i] == 0 && s[i + 1] == 3 && s[i + 2] == 2) // ATG
        };
        if let Some(st) = start {
            if e - st >= min_len {
                out.push(to_fwd(strand, f, st, e, seq.len()));
            }
        }
    }
    out
}

/// Codons of an ORF in its own reading direction.
fn codons_of(seq: &[u8], rc: &[u8], o: &Orf) -> Vec<usize> {
    let (s, a, b) = if o.frame < 3 { (seq, o.start, o.end) } else { (rc, seq.len() - o.end, seq.len() - o.start) };
    (a..b).step_by(3).filter(|&i| i + 3 <= b).filter_map(|i| codon_index(&s[i..i + 3])).collect()
}

/// Codons covering [start, end) read in frame `fr` (0..3 offset, either strand).
fn shadow_codons(seq: &[u8], rc: &[u8], start: usize, end: usize, strand: usize, fr: usize) -> Vec<usize> {
    let (s, a, b) = if strand == 0 { (seq, start, end) } else { (rc, seq.len() - end, seq.len() - start) };
    let first = a + (3 + fr - a % 3) % 3;
    (first..b).step_by(3).filter(|&i| i + 3 <= b).filter_map(|i| codon_index(&s[i..i + 3])).collect()
}

fn overlaps(a: &Orf, b: &Orf) -> bool {
    a.start < b.end && b.start < a.end
}

pub fn scored_orfs(seq: &[u8], min_len: usize, min_llr: f64) -> Vec<Orf> {
    let mut cands = candidates(seq, min_len);
    cands.sort_by(|a, b| (b.end - b.start).cmp(&(a.end - a.start)).then(a.start.cmp(&b.start)));
    let mut train: Vec<usize> = Vec::new();
    for (i, c) in cands.iter().enumerate() {
        if train.iter().all(|&t| !overlaps(&cands[t], c)) {
            train.push(i);
        }
    }
    let rc = revcomp(seq);
    let mut fg = vec![1.0f64; 64];
    let mut bg = vec![1.0f64; 64];
    for &t in &train {
        let o = &cands[t];
        for c in codons_of(seq, &rc, o) {
            fg[c] += 1.0;
        }
        let own_strand = (o.frame >= 3) as usize;
        let own_fr = if own_strand == 0 { o.start % 3 } else { (seq.len() - o.end) % 3 };
        for strand in 0..2 {
            for fr in 0..3 {
                if strand == own_strand && fr == own_fr {
                    continue;
                }
                for c in shadow_codons(seq, &rc, o.start, o.end, strand, fr) {
                    bg[c] += 1.0;
                }
            }
        }
    }
    let (fs, bs): (f64, f64) = (fg.iter().sum(), bg.iter().sum());
    let llr: Vec<f64> = (0..64).map(|c| ((fg[c] / fs) / (bg[c] / bs)).ln()).collect();
    let n_train_codons: f64 = fs - 64.0;
    cands
        .iter()
        .enumerate()
        .filter(|(i, c)| {
            if train.contains(i) {
                return true;
            }
            // too little sequence to learn usage from: keep only the training set
            if n_train_codons < 200.0 {
                return false;
            }
            if c.end - c.start < min_len * 3 / 2 {
                return false;
            }
            let cs = codons_of(seq, &rc, c);
            let m = cs.iter().map(|&x| llr[x]).sum::<f64>() / cs.len().max(1) as f64;
            m >= min_llr
        })
        .map(|(_, c)| c.clone())
        .collect()
}

pub fn call_orfs(seq: &[u8], min_len: usize, caller: OrfCaller) -> Vec<Orf> {
    match caller {
        OrfCaller::Union => union_orfs(seq, min_len),
        OrfCaller::Atg => candidates(seq, min_len),
        OrfCaller::Scored => scored_orfs(seq, min_len, DEFAULT_MIN_LLR),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode;

    #[test]
    fn union_finds_stop_free_run() {
        let s = encode(b"ATGAAACCCGGGTTTAAACCCGGGTTTTAA");
        let orfs = union_orfs(&s, 18);
        assert!(orfs.iter().any(|o| o.frame == 0 && o.start == 0 && o.end == 27));
    }

    #[test]
    fn shadow_orf_inside_gene_not_called_dual() {
        // a long gene whose reverse strand carries a long ATG-started open
        // stretch: the scored caller must not call the gene dual-coding
        let mut x: u64 = 9;
        let mut g = vec![0u8, 3, 2];
        while g.len() < 3000 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let c = ((x >> 33) % 64) as usize;
            // skewed usage, as real genes have: favour codons ending in G/C
            if is_stop(c) || (c % 4 < 2 && (x >> 20) % 3 != 0) {
                continue;
            }
            g.extend_from_slice(&[crate::genetic::base_of(c, 0), crate::genetic::base_of(c, 1), crate::genetic::base_of(c, 2)]);
        }
        g.extend_from_slice(&[3, 0, 0]);
        let scored = scored_orfs(&g, 300, DEFAULT_MIN_LLR);
        let dual_nt: usize = {
            let mut n = vec![0u8; g.len()];
            for o in &scored {
                for v in &mut n[o.start..o.end] {
                    *v += 1;
                }
            }
            n.iter().filter(|&&k| k >= 2).count()
        };
        assert!(dual_nt < 150, "scored caller called {} nt dual", dual_nt);
        assert!(scored.iter().any(|o| o.frame == 0 && o.start == 0));
    }

    #[test]
    fn candidates_need_atg_unless_open() {
        // frame 0: TAA then stop-free run without ATG; must not be a candidate
        let s = encode(b"TAACCCGGGTTTCCCGGGTTTCCCGGGTAA");
        assert!(candidates(&s, 18).iter().all(|o| !(o.frame == 0 && o.start == 3)));
        let s = encode(b"TAACCCATGTTTCCCGGGTTTCCCGGGTTTTAA");
        assert!(candidates(&s, 18).iter().any(|o| o.frame == 0 && o.start == 6));
    }
}
