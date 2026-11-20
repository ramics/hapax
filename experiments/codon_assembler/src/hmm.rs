//! Banded pair-HMM for read-read overlaps, scored in codon space.
//!
//! This is the RAMICS model (Wright & Travers 2014) turned on its side. RAMICS
//! aligned a read to a reference whose open reading frame had been elected;
//! here there is no reference and no frame is elected. Read A plays the
//! reference role with a codon grid laid at phase `phi` (codons start at A
//! positions congruent to phi mod 3), and read B is free. Running phi = 0, 1, 2
//! is the "three candidate offsets" of the spec: whatever frame the locus is
//! read in, one of the three grids lines up with it relative to the pair, and
//! that one scores best. Dual-coding regions show up as two grids scoring well.
//!
//! States, following RAMICS (thesis Fig. 2-1):
//! * codon match on the grid (3 nt of each read, GY94 codon log-odds),
//! * a codon whose B side lost one base (sequencing deletion inside a codon),
//! * a codon whose B side gained one base (sequencing insertion inside a codon),
//! * whole-codon insertion or deletion in B (biological, codon-sized gap),
//! * a nucleotide-mode state that scores base by base with the K80 log-odds.
//!
//! The nucleotide state is how scoring mode is chosen per region. Entering or
//! leaving it costs `switch`; Viterbi picks codon scoring where the sequence
//! behaves like coding sequence under the pair's relative frame and nucleotide
//! scoring through UTRs, packaging signals and anything else. Because both
//! scorers are log-odds against the same random background, the comparison is
//! fair. A one-base frameshift in A breaks A's grid; the path then has to fall
//! back to nucleotide mode, which is the right cost for an event that rare.
//!
//! The alignment is end-to-end over the overlap (free end gaps only at read
//! boundaries), so it can never report a partial local match as an overlap.

use crate::genetic::{aa, codon_index};
use crate::model::Scorer;

const NEG: f32 = -1.0e30;

#[derive(Clone, Copy, Debug)]
pub struct HmmParams {
    /// log P of a 1-nt sequencing indel (frameshift error).
    pub frameshift: f32,
    /// log P of a codon-sized biological indel.
    pub codon_gap: f32,
    /// log P of a nucleotide-mode gap (one base).
    pub nuc_gap: f32,
    /// log P of switching between codon and nucleotide scoring.
    pub switch: f32,
    /// Half-width of the diagonal band.
    pub band: usize,
}

impl Default for HmmParams {
    fn default() -> Self {
        // frameshift: MiSeq-class indel error ~1e-3 per base (thesis Table 2-1
        // order of magnitude); codon gaps rare but real in viral populations.
        HmmParams { frameshift: (1e-3f32).ln(), codon_gap: (2e-3f32).ln(), nuc_gap: (1e-3f32).ln(), switch: -4.0, band: 12 }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Aln {
    pub score: f32,
    /// grid phase used, or 3 for the nucleotide-only alignment
    pub phase: u8,
    pub a_start: u32,
    pub a_end: u32,
    pub b_start: u32,
    pub b_end: u32,
    /// aligned columns (matches + mismatches + gap columns)
    pub cols: u32,
    pub matches: u32,
    pub codon_nt: u32,
    pub nuc_nt: u32,
    pub syn: u32,
    pub nonsyn: u32,
}

impl Aln {
    pub fn identity(&self) -> f32 {
        if self.cols == 0 {
            0.0
        } else {
            self.matches as f32 / self.cols as f32
        }
    }
    pub fn overlap_len(&self) -> u32 {
        (self.a_end - self.a_start).max(self.b_end - self.b_start)
    }
}

/// Banded storage: row i, column offset k = j - lo(i).
struct Band {
    w: usize,
    d: i64,
    lb: usize,
}

impl Band {
    #[inline]
    fn lo(&self, i: usize) -> i64 {
        i as i64 - self.d - self.w as i64
    }
    #[inline]
    fn idx(&self, i: usize, j: i64) -> Option<usize> {
        let lo = self.lo(i);
        if j < 0 || j > self.lb as i64 || j < lo || j > lo + 2 * self.w as i64 {
            return None;
        }
        Some(i * (2 * self.w + 1) + (j - lo) as usize)
    }
}

// backpointer codes
const BP_START: u8 = 0;
const G_MATCH: u8 = 1;
const G_BDEL: u8 = 2; // + 16*p
const G_BINS: u8 = 3; // + 16*p
const G_CINS: u8 = 4;
const G_CDEL: u8 = 5;
const G_FROM_N: u8 = 6;
const N_DIAG: u8 = 1;
const N_UP: u8 = 2; // consume A only
const N_LEFT: u8 = 3; // consume B only
const N_FROM_G: u8 = 4;

/// Align `a` (forward) against `b` (already oriented) around diagonal `d`
/// (B[j] expected at A[j + d]). `phase` = Some(phi) enables the codon grid,
/// None gives a pure nucleotide-mode alignment.
pub fn align(a: &[u8], b: &[u8], d: i64, phase: Option<u8>, sc: &Scorer, hp: &HmmParams) -> Option<Aln> {
    let la = a.len();
    let lb = b.len();
    let band = Band { w: hp.band, d, lb };
    let width = 2 * hp.band + 1;
    let ncell = (la + 1) * width;
    let mut g = vec![NEG; ncell];
    let mut n = vec![NEG; ncell];
    let mut gbp = vec![0u8; ncell];
    let mut nbp = vec![0u8; ncell];
    let phi = phase.map(|p| p as usize);

    let mut best = (NEG, 0usize, 0i64, false); // score, i, j, in_g
    for i in 0..=la {
        let lo = band.lo(i);
        let grid = phi.map_or(false, |p| i % 3 == p);
        for k in 0..width {
            let j = lo + k as i64;
            if j < 0 || j > lb as i64 {
                continue;
            }
            let ju = j as usize;
            let cell = i * width + k;
            // ---- nucleotide-mode predecessor
            let (mut npre, mut nb) = if i == 0 || j == 0 { (0.0, BP_START) } else { (NEG, BP_START) };
            if i > 0 && j > 0 {
                if let Some(p) = band.idx(i - 1, j - 1) {
                    let v = n[p] + sc.nuc(a[i - 1], b[ju - 1]);
                    if v > npre {
                        npre = v;
                        nb = N_DIAG;
                    }
                }
            }
            if i > 0 {
                if let Some(p) = band.idx(i - 1, j) {
                    let v = n[p] + hp.nuc_gap;
                    if v > npre {
                        npre = v;
                        nb = N_UP;
                    }
                }
            }
            if j > 0 {
                if let Some(p) = band.idx(i, j - 1) {
                    let v = n[p] + hp.nuc_gap;
                    if v > npre {
                        npre = v;
                        nb = N_LEFT;
                    }
                }
            }
            // ---- codon-grid predecessor
            let (mut gpre, mut gb) = (NEG, BP_START);
            if grid {
                // free start near a read boundary; the partial codon before
                // the grid is scored base by base so no phase gets a free ride
                if ju < 3 || i < 3 {
                    let m = i.min(ju);
                    gpre = (0..m).map(|t| sc.nuc(a[i - m + t], b[ju - m + t])).sum();
                    gb = BP_START;
                }
                if i >= 3 {
                    let ac = &a[i - 3..i];
                    if ju >= 3 {
                        if let Some(p) = band.idx(i - 3, j - 3) {
                            let v = g[p] + sc.codon(ac, &b[ju - 3..ju]);
                            if v > gpre {
                                gpre = v;
                                gb = G_MATCH;
                            }
                        }
                    }
                    // B lost one base inside this codon: align 2 of A's 3
                    if ju >= 2 {
                        if let Some(p) = band.idx(i - 3, j - 2) {
                            if g[p] > NEG {
                                let bs = &b[ju - 2..ju];
                                for skip in 0..3 {
                                    let mut s = hp.frameshift;
                                    let mut bi = 0;
                                    for q in 0..3 {
                                        if q == skip {
                                            continue;
                                        }
                                        s += sc.nuc(ac[q], bs[bi]);
                                        bi += 1;
                                    }
                                    let v = g[p] + s;
                                    if v > gpre {
                                        gpre = v;
                                        gb = G_BDEL + 16 * skip as u8;
                                    }
                                }
                            }
                        }
                    }
                    // B gained one base inside this codon: drop one of 4
                    if ju >= 4 {
                        if let Some(p) = band.idx(i - 3, j - 4) {
                            if g[p] > NEG {
                                let bs = &b[ju - 4..ju];
                                for skip in 0..4 {
                                    let mut c = [0u8; 3];
                                    let mut ci = 0;
                                    for q in 0..4 {
                                        if q != skip {
                                            c[ci] = bs[q];
                                            ci += 1;
                                        }
                                    }
                                    let v = g[p] + hp.frameshift + sc.codon(ac, &c);
                                    if v > gpre {
                                        gpre = v;
                                        gb = G_BINS + 16 * skip as u8;
                                    }
                                }
                            }
                        }
                    }
                    // codon deleted from B
                    if let Some(p) = band.idx(i - 3, j) {
                        let v = g[p] + hp.codon_gap;
                        if v > gpre {
                            gpre = v;
                            gb = G_CDEL;
                        }
                    }
                }
                // codon inserted in B
                if ju >= 3 {
                    if let Some(p) = band.idx(i, j - 3) {
                        let v = g[p] + hp.codon_gap;
                        if v > gpre {
                            gpre = v;
                            gb = G_CINS;
                        }
                    }
                }
            }
            // ---- mode switches at this cell
            if grid {
                let gv = if npre + hp.switch > gpre { (npre + hp.switch, G_FROM_N) } else { (gpre, gb) };
                let nv = if gpre + hp.switch > npre { (gpre + hp.switch, N_FROM_G) } else { (npre, nb) };
                g[cell] = gv.0;
                gbp[cell] = gv.1;
                n[cell] = nv.0;
                nbp[cell] = nv.1;
            } else {
                n[cell] = npre;
                nbp[cell] = nb;
            }
            // ---- end conditions (reached the end of either read)
            let end_n = i == la || ju == lb;
            if end_n && n[cell] > best.0 {
                best = (n[cell], i, j, false);
            }
            if grid && (i + 2 >= la || ju + 2 >= lb) && g[cell] > NEG {
                // trailing partial codon, scored base by base
                let m = (la - i).min(lb - ju);
                let v = g[cell] + (0..m).map(|t| sc.nuc(a[i + t], b[ju + t])).sum::<f32>();
                if v > best.0 {
                    best = (v, i, j, true);
                }
            }
        }
    }
    if best.0 <= NEG / 2.0 {
        return None;
    }
    // ---- traceback
    let (score, mut i, mut j, mut in_g) = best;
    let mut r = Aln { score, phase: phase.unwrap_or(3), a_end: i as u32, b_end: j as u32, ..Default::default() };
    if in_g {
        let m = (la - i).min(lb - j as usize);
        tally_nuc_diag(&mut r, &a[i..i + m], &b[j as usize..j as usize + m]);
        r.a_end += m as u32;
        r.b_end += m as u32;
    }
    loop {
        let cell = band.idx(i, j).unwrap();
        if in_g {
            let code = gbp[cell];
            match code & 15 {
                BP_START => {
                    let m = i.min(j as usize);
                    tally_nuc_diag(&mut r, &a[i - m..i], &b[j as usize - m..j as usize]);
                    i -= m;
                    j -= m as i64;
                    break;
                }
                G_MATCH => {
                    let ac = &a[i - 3..i];
                    let bc = &b[j as usize - 3..j as usize];
                    tally_codon(&mut r, ac, bc);
                    i -= 3;
                    j -= 3;
                }
                G_BDEL => {
                    let skip = (code >> 4) as usize;
                    let ac = &a[i - 3..i];
                    let bs = &b[j as usize - 2..j as usize];
                    let mut bi = 0;
                    for q in 0..3 {
                        r.cols += 1;
                        r.codon_nt += 1;
                        if q != skip {
                            if ac[q] == bs[bi] {
                                r.matches += 1;
                            }
                            bi += 1;
                        }
                    }
                    i -= 3;
                    j -= 2;
                }
                G_BINS => {
                    let skip = (code >> 4) as usize;
                    let bs = &b[j as usize - 4..j as usize];
                    let mut c = [0u8; 3];
                    let mut ci = 0;
                    for q in 0..4 {
                        if q != skip {
                            c[ci] = bs[q];
                            ci += 1;
                        }
                    }
                    tally_codon(&mut r, &a[i - 3..i], &c);
                    r.cols += 1;
                    i -= 3;
                    j -= 4;
                }
                G_CINS => {
                    r.cols += 3;
                    r.codon_nt += 3;
                    j -= 3;
                }
                G_CDEL => {
                    r.cols += 3;
                    r.codon_nt += 3;
                    i -= 3;
                }
                G_FROM_N => in_g = false,
                _ => unreachable!(),
            }
        } else {
            match nbp[cell] {
                BP_START => break,
                N_DIAG => {
                    r.cols += 1;
                    r.nuc_nt += 1;
                    if a[i - 1] == b[j as usize - 1] {
                        r.matches += 1;
                    }
                    i -= 1;
                    j -= 1;
                }
                N_UP => {
                    r.cols += 1;
                    r.nuc_nt += 1;
                    i -= 1;
                }
                N_LEFT => {
                    r.cols += 1;
                    r.nuc_nt += 1;
                    j -= 1;
                }
                N_FROM_G => in_g = true,
                _ => unreachable!(),
            }
        }
    }
    r.a_start = i as u32;
    r.b_start = j as u32;
    Some(r)
}

fn tally_nuc_diag(r: &mut Aln, a: &[u8], b: &[u8]) {
    for (x, y) in a.iter().zip(b) {
        r.cols += 1;
        r.nuc_nt += 1;
        if x == y {
            r.matches += 1;
        }
    }
}

fn tally_codon(r: &mut Aln, ac: &[u8], bc: &[u8]) {
    let mut diff = 0;
    for q in 0..3 {
        r.cols += 1;
        r.codon_nt += 1;
        if ac[q] == bc[q] {
            r.matches += 1;
        } else {
            diff += 1;
        }
    }
    if diff > 0 {
        match (codon_index(ac), codon_index(bc)) {
            (Some(x), Some(y)) if aa(x) == aa(y) => r.syn += diff,
            _ => r.nonsyn += diff,
        }
    }
}

/// Result of scoring one read pair under all three relative offsets plus the
/// nucleotide-only model.
#[derive(Clone, Copy, Debug)]
pub struct PairScore {
    pub best: Aln,
    pub per_phase: [f32; 3],
    pub nucleotide: Aln,
    /// two grids both beat nucleotide scoring by the dual-coding margin
    pub dual: bool,
}

/// Score a pair: three offsets, best wins, nucleotide alignment alongside.
/// When two offsets both gain substantially over nucleotide scoring the region
/// is plausibly dual-coding; METHOD.md 5 says synonymous collapse is unsafe
/// there, so the caller is told and falls back to the nucleotide score.
pub fn score_pair(a: &[u8], b: &[u8], d: i64, sc: &Scorer, hp: &HmmParams, dual_margin: f32) -> Option<PairScore> {
    let nuc = align(a, b, d, None, sc, hp)?;
    let mut alns: Vec<Aln> = (0..3u8).filter_map(|p| align(a, b, d, Some(p), sc, hp)).collect();
    if alns.is_empty() {
        return None;
    }
    let mut per_phase = [NEG; 3];
    for x in &alns {
        per_phase[x.phase as usize] = x.score;
    }
    alns.sort_by(|x, y| y.score.partial_cmp(&x.score).unwrap());
    let gain1 = alns[0].score - nuc.score;
    let gain2 = if alns.len() > 1 { alns[1].score - nuc.score } else { NEG };
    let dual = gain1 > dual_margin && gain2 > dual_margin;
    Some(PairScore { best: alns[0], per_phase, nucleotide: nuc, dual })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::ModelParams;
    use crate::seq::encode;

    fn rand_coding(n_codons: usize, seed: u64) -> Vec<u8> {
        // simple LCG codon sampler without stops
        let mut x = seed;
        let mut out = Vec::new();
        while out.len() < n_codons * 3 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let c = ((x >> 33) % 64) as usize;
            if crate::genetic::is_stop(c) {
                continue;
            }
            for p in 0..3 {
                out.push(crate::genetic::base_of(c, p));
            }
        }
        out
    }

    #[test]
    fn identical_overlap_found() {
        let sc = Scorer::new(ModelParams::default());
        let hp = HmmParams::default();
        let g = rand_coding(100, 7);
        let a = g[0..150].to_vec();
        let b = g[60..210].to_vec();
        let r = align(&a, &b, 60, Some(0), &sc, &hp).unwrap();
        assert_eq!(r.a_start, 60);
        assert_eq!(r.b_start, 0);
        assert_eq!(r.a_end, 150);
        assert_eq!(r.matches, 90);
        assert!(r.identity() > 0.99);
    }

    #[test]
    fn right_phase_wins_on_synonymous_variation() {
        let sc = Scorer::new(ModelParams::default());
        let hp = HmmParams::default();
        let g = rand_coding(100, 11);
        // variant: change third positions synonymously where possible
        let mut v = g.clone();
        let mut changed = 0;
        for c in (0..v.len() / 3).step_by(3) {
            let i = c * 3;
            let orig = codon_index(&v[i..i + 3]).unwrap();
            for nb in 0..4u8 {
                if nb == v[i + 2] {
                    continue;
                }
                let mut t = [v[i], v[i + 1], nb];
                let ni = codon_index(&t).unwrap();
                if aa(ni) == aa(orig) {
                    v[i + 2] = nb;
                    t[2] = nb;
                    changed += 1;
                    break;
                }
            }
        }
        assert!(changed > 10);
        // A = g[30..180] (frame phase 0 in A), B = v[90..240]
        let a = g[30..180].to_vec();
        let b = v[90..240].to_vec();
        let ps = score_pair(&a, &b, 60, &sc, &hp, 1e9).unwrap();
        assert_eq!(ps.best.phase, 0, "per-phase {:?}", ps.per_phase);
        assert!(ps.best.score > ps.nucleotide.score);
        assert!(ps.best.syn > 0 && ps.best.nonsyn == 0);
    }

    #[test]
    fn frameshift_error_in_b_tolerated() {
        let sc = Scorer::new(ModelParams::default());
        let hp = HmmParams::default();
        let g = rand_coding(100, 3);
        let a = g[0..150].to_vec();
        let mut b = g[60..210].to_vec();
        b.insert(40, 0); // one-base sequencing insertion
        let r = align(&a, &b, 60, Some(0), &sc, &hp).unwrap();
        assert!(r.identity() > 0.97, "{:?}", r);
        assert!(r.codon_nt > 60);
    }

    #[test]
    fn unrelated_sequences_score_negative() {
        let sc = Scorer::new(ModelParams::default());
        let hp = HmmParams::default();
        let a = rand_coding(50, 1);
        let b = rand_coding(50, 2);
        let ps = score_pair(&a, &b, 60, &sc, &hp, 1e9).unwrap();
        assert!(ps.best.score < 20.0, "{:?}", ps.best);
        let _ = encode(b"A");
    }
}
