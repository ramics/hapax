//! Substitution models that turn a pair of aligned bases or codons into a
//! log-odds score.
//!
//! Both scorers answer the same question -- how much more likely is this pair
//! under "two reads from one viral population" than under "two unrelated
//! sequences" -- so their scores are on the same scale and can be compared
//! region by region. That comparability is what lets the overlap HMM choose
//! codon or nucleotide scoring per region rather than once per run.
//!
//! The codon model is Goldman & Yang (1994): single-nucleotide instantaneous
//! changes between sense codons, transitions weighted by kappa, non-synonymous
//! changes weighted by omega. The default omega of 0.23 is the mean backbone
//! dN/dS for intra-patient HIV-1 env (Lemey et al. 2007), the same value the
//! author's thesis used to drive EvolveAGene. RAMICS used the empirical Kosiol
//! (ECM) matrix; a mechanistic matrix is used here because the parameters that
//! matter for synonymous collapse (omega, divergence) are then explicit knobs
//! that experiments can sweep, rather than baked into an empirical table.
//!
//! Sequencing error is convolved into both models: an observed mismatch may be
//! machine error rather than biology, at rate `sub_err` per base per read
//! (thesis default 0.5%).

use crate::genetic::{base_of, is_stop, is_transition, sense_codons};

#[derive(Clone, Debug)]
pub struct ModelParams {
    /// Expected substitutions per nucleotide site separating two reads of one
    /// population (pairwise distance, not distance to the root).
    pub divergence: f64,
    pub kappa: f64,
    pub omega: f64,
    pub sub_err: f64,
}

impl Default for ModelParams {
    fn default() -> Self {
        ModelParams { divergence: 0.10, kappa: 2.0, omega: 0.23, sub_err: 0.005 }
    }
}

/// Precomputed log-odds tables.
#[derive(Clone)]
pub struct Scorer {
    /// codon[a*64+b]: log-odds for codon a in read A aligned to codon b in read B.
    pub codon: Vec<f32>,
    /// nuc[a*4+b]: log-odds for base a aligned to base b.
    pub nuc: [f32; 16],
    pub params: ModelParams,
}

impl Scorer {
    pub fn new(p: ModelParams) -> Scorer {
        let codon = codon_log_odds(&p);
        let nuc = nuc_log_odds(&p);
        Scorer { codon, nuc, params: p }
    }

    #[inline]
    pub fn nuc(&self, a: u8, b: u8) -> f32 {
        if a > 3 || b > 3 {
            return 0.0;
        }
        self.nuc[(a * 4 + b) as usize]
    }

    /// Score a codon pair given as 3-base slices. Codons containing N fall back
    /// to the sum of nucleotide scores, which is neutral (0) for the N positions.
    #[inline]
    pub fn codon(&self, a: &[u8], b: &[u8]) -> f32 {
        match (crate::genetic::codon_index(a), crate::genetic::codon_index(b)) {
            (Some(x), Some(y)) => self.codon[x * 64 + y],
            _ => self.nuc(a[0], b[0]) + self.nuc(a[1], b[1]) + self.nuc(a[2], b[2]),
        }
    }
}

/// GY94 rate matrix over the 61 sense codons with uniform codon frequencies,
/// scaled so that one time unit is one expected substitution per nucleotide.
fn gy94_q(kappa: f64, omega: f64) -> (Vec<usize>, Vec<f64>) {
    let sense = sense_codons();
    let n = sense.len();
    let mut q = vec![0.0; n * n];
    for (i, &a) in sense.iter().enumerate() {
        for (j, &b) in sense.iter().enumerate() {
            if i == j {
                continue;
            }
            let diffs: Vec<usize> = (0..3).filter(|&p| base_of(a, p) != base_of(b, p)).collect();
            if diffs.len() != 1 {
                continue;
            }
            let p = diffs[0];
            let mut r = 1.0 / n as f64;
            if is_transition(base_of(a, p), base_of(b, p)) {
                r *= kappa;
            }
            if crate::genetic::aa(a) != crate::genetic::aa(b) {
                r *= omega;
            }
            q[i * n + j] = r;
        }
    }
    let mut total = 0.0;
    for i in 0..n {
        let s: f64 = (0..n).map(|j| q[i * n + j]).sum();
        q[i * n + i] = -s;
        total += s / n as f64; // uniform stationary distribution
    }
    // total = expected substitutions per codon per unit time; per nucleotide is /3
    let scale = 3.0 / total;
    for v in q.iter_mut() {
        *v *= scale;
    }
    (sense, q)
}

fn matmul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    c
}

/// Matrix exponential by scaling and squaring with a Taylor series. The
/// matrices here are small (61x61) and well conditioned, so this is accurate
/// to well below the precision the scores need.
pub fn expm(m: &[f64], n: usize) -> Vec<f64> {
    let norm = (0..n).map(|i| (0..n).map(|j| m[i * n + j].abs()).sum::<f64>()).fold(0.0, f64::max);
    let mut s = 0;
    while norm / (1u64 << s) as f64 > 0.25 {
        s += 1;
    }
    let scale = 1.0 / (1u64 << s) as f64;
    let a: Vec<f64> = m.iter().map(|v| v * scale).collect();
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        result[i * n + i] = 1.0;
    }
    let mut term = result.clone();
    for k in 1..=14 {
        term = matmul(&term, &a, n);
        for v in term.iter_mut() {
            *v /= k as f64;
        }
        for (r, t) in result.iter_mut().zip(term.iter()) {
            *r += t;
        }
    }
    for _ in 0..s {
        result = matmul(&result, &result, n);
    }
    result
}

/// Per-base error channel applied to one codon position: stay with 1-e,
/// move to each other base with e/3.
fn err_matrix_nuc(e: f64) -> [f64; 16] {
    let mut m = [0.0; 16];
    for a in 0..4 {
        for b in 0..4 {
            m[a * 4 + b] = if a == b { 1.0 - e } else { e / 3.0 };
        }
    }
    m
}

fn codon_log_odds(p: &ModelParams) -> Vec<f32> {
    let (sense, q) = gy94_q(p.kappa, p.omega);
    let n = sense.len();
    // Pairwise distance d between two tips is time d/2 from each to their
    // ancestor, but under a reversible model the joint distribution of the
    // two tips is pi_a P(d)_ab directly.
    let qt: Vec<f64> = q.iter().map(|v| v * p.divergence).collect();
    let pt = expm(&qt, n);
    let mut joint = vec![0.0f64; 64 * 64];
    for (i, &a) in sense.iter().enumerate() {
        for (j, &b) in sense.iter().enumerate() {
            joint[a * 64 + b] = pt[i * n + j] / n as f64;
        }
    }
    // convolve observation error from both reads (combined rate ~2e)
    let e = err_matrix_nuc(2.0 * p.sub_err);
    let mut obs = vec![0.0f64; 64 * 64];
    for a in 0..64 {
        for b in 0..64 {
            let j = joint[a * 64 + b];
            if j == 0.0 {
                continue;
            }
            for c in 0..64 {
                let mut pe = 1.0;
                for pos in 0..3 {
                    pe *= e[(base_of(b, pos) * 4 + base_of(c, pos)) as usize];
                }
                obs[a * 64 + c] += j * pe;
            }
        }
    }
    // also convolve A side symmetrically so stop codons in A get error mass
    let mut obs2 = vec![0.0f64; 64 * 64];
    let e1 = err_matrix_nuc(p.sub_err);
    for a in 0..64 {
        for c in 0..64 {
            let mut pe = 1.0;
            for pos in 0..3 {
                pe *= e1[(base_of(a, pos) * 4 + base_of(c, pos)) as usize];
            }
            if pe < 1e-12 {
                continue;
            }
            for b in 0..64 {
                obs2[c * 64 + b] += pe * obs[a * 64 + b];
            }
        }
    }
    let bg = 1.0 / 64.0;
    obs2.iter().map(|&v| (v.max(1e-30) / (bg * bg)).ln() as f32).collect()
}

fn nuc_log_odds(p: &ModelParams) -> [f32; 16] {
    // K80 with the same kappa; P(transition) and P(transversion) at distance d.
    let k = p.kappa;
    let d = p.divergence;
    // rate parameterisation: alpha (ts), beta (each tv), d = (alpha + 2 beta) t
    let beta = 1.0 / (k + 2.0);
    let alpha = k * beta;
    let pts = 0.25 + 0.25 * (-4.0 * beta * d).exp() - 0.5 * (-2.0 * (alpha + beta) * d).exp();
    let ptv = 0.25 - 0.25 * (-4.0 * beta * d).exp();
    let pid = 1.0 - pts - 2.0 * ptv;
    let mut m = [[0.0f64; 4]; 4];
    for a in 0..4u8 {
        for b in 0..4u8 {
            m[a as usize][b as usize] = if a == b {
                pid
            } else if is_transition(a, b) {
                pts
            } else {
                ptv
            } * 0.25;
        }
    }
    let e = err_matrix_nuc(2.0 * p.sub_err);
    let mut out = [0.0f32; 16];
    for a in 0..4 {
        for c in 0..4 {
            let mut s = 0.0;
            for b in 0..4 {
                s += m[a][b] * e[b * 4 + c];
            }
            out[a * 4 + c] = (s / (0.25 * 0.25)).ln() as f32;
        }
    }
    out
}

/// Is codon index `c` a stop? Re-exported for scorers that need it.
pub fn stop(c: usize) -> bool {
    is_stop(c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode;

    #[test]
    fn expm_rows_sum_to_one() {
        let (_s, q) = gy94_q(2.0, 0.23);
        let n = 61;
        let qt: Vec<f64> = q.iter().map(|v| v * 0.2).collect();
        let p = expm(&qt, n);
        for i in 0..n {
            let s: f64 = (0..n).map(|j| p[i * n + j]).sum();
            assert!((s - 1.0).abs() < 1e-9, "row {} sums to {}", i, s);
        }
        // expected per-nucleotide divergence matches the time scale at small t
        let qt: Vec<f64> = q.iter().map(|v| v * 0.001).collect();
        let p = expm(&qt, n);
        let stay: f64 = (0..n).map(|i| p[i * n + i]).sum::<f64>() / n as f64;
        assert!(((1.0 - stay) / 3.0 - 0.001).abs() < 1e-4);
    }

    #[test]
    fn synonymous_beats_nonsynonymous() {
        let s = Scorer::new(ModelParams::default());
        let c = |x: &[u8]| encode(x);
        let same = s.codon(&c(b"CTG"), &c(b"CTG"));
        let syn = s.codon(&c(b"CTG"), &c(b"CTA")); // Leu -> Leu, transition
        let non = s.codon(&c(b"CTG"), &c(b"CCG")); // Leu -> Pro, transition
        let stop = s.codon(&c(b"TGG"), &c(b"TAG")); // Trp -> stop
        assert!(same > syn && syn > non && non > stop, "{} {} {} {}", same, syn, non, stop);
        // the codon model forgives a synonymous change more than the nucleotide
        // model forgives the same single mismatch: this is the collapse
        let nuc_syn = s.nuc(1, 1) + s.nuc(3, 3) + s.nuc(2, 0); // CTG vs CTA, base by base
        assert!(syn > nuc_syn, "codon syn {} vs nuc {}", syn, nuc_syn);
    }

    #[test]
    fn nucleotide_identity_positive() {
        let s = Scorer::new(ModelParams::default());
        for a in 0..4 {
            assert!(s.nuc(a, a) > 1.0);
            for b in 0..4 {
                if a != b {
                    assert!(s.nuc(a, b) < 0.0);
                }
            }
        }
    }
}
