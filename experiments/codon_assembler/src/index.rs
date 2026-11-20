//! Minimizer seeding index.
//!
//! All-against-all overlap is quadratic in any language; the index is what
//! makes candidate generation near-linear. Each read is reduced to its (w,k)
//! minimizers -- the smallest-hash k-mer in every window of w consecutive
//! k-mers -- and reads sharing a minimizer land in the same bucket. Only pairs
//! that share buckets are ever compared.
//!
//! Two seed shapes are supported:
//!
//! * `Contiguous`: ordinary exact k-mers. This is the spec's tier-1 seed.
//! * `CodonSpaced`: a spaced seed with mask `110110...` that ignores every third
//!   position. When the mask phase lines up with a reading frame, the ignored
//!   positions are the wobble positions where most within-population variation
//!   sits, so two reads from different variants of one swarm still share seeds.
//!   Every position is seeded, so every phase is tried and no frame is elected.
//!   The seed still demands exact nucleotide identity on the positions it reads.
//!
//! Hashes are canonical (min of forward and reverse-complement), so reads from
//! either strand meet in the same bucket.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeedShape {
    Contiguous,
    CodonSpaced,
}

#[derive(Clone, Copy, Debug)]
pub struct SeedParams {
    pub shape: SeedShape,
    /// Number of informative positions in the seed.
    pub k: usize,
    pub w: usize,
    /// Buckets larger than this are dropped (repeats / low complexity).
    pub max_occ: usize,
}

impl Default for SeedParams {
    fn default() -> Self {
        SeedParams { shape: SeedShape::Contiguous, k: 15, w: 5, max_occ: 2000 }
    }
}

impl SeedParams {
    /// Span on the sequence covered by one seed.
    pub fn span(&self) -> usize {
        match self.shape {
            SeedShape::Contiguous => self.k,
            // k informative positions with mask 110: span = k + k/2 (rounded so the
            // final position is informative)
            SeedShape::CodonSpaced => {
                let full = self.k / 2;
                let rem = self.k % 2;
                full * 3 + rem
            }
        }
    }

    fn mask(&self) -> Vec<bool> {
        let span = self.span();
        match self.shape {
            SeedShape::Contiguous => vec![true; span],
            SeedShape::CodonSpaced => (0..span).map(|i| i % 3 != 2).collect(),
        }
    }
}

#[inline]
fn mix64(mut x: u64) -> u64 {
    // invertible integer hash (splitmix64 finaliser)
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^= x >> 31;
    x
}

/// A seed occurrence: `pos` is the leftmost base of the seed span on the
/// forward strand of the read; `rev` is true when the canonical form came from
/// the reverse complement.
#[derive(Clone, Copy, Debug)]
pub struct Seed {
    pub hash: u64,
    pub pos: u32,
    pub rev: bool,
}

/// Compute (w,k) minimizers of one read.
pub fn minimizers(seq: &[u8], p: &SeedParams) -> Vec<Seed> {
    let span = p.span();
    if seq.len() < span {
        return Vec::new();
    }
    let mask = p.mask();
    let n = seq.len() - span + 1;
    // all canonical seed hashes; None where the window contains N
    let mut all: Vec<Option<(u64, bool)>> = Vec::with_capacity(n);
    for i in 0..n {
        let win = &seq[i..i + span];
        if win.iter().any(|&b| b > 3) {
            all.push(None);
            continue;
        }
        let mut f: u64 = 0;
        let mut r: u64 = 0;
        for (j, &m) in mask.iter().enumerate() {
            if m {
                f = (f << 2) | win[j] as u64;
                // reverse complement read with the same mask
                r = (r << 2) | (3 - win[span - 1 - j]) as u64;
            }
        }
        let (hf, hr) = (mix64(f), mix64(r));
        all.push(Some(if hf <= hr { (hf, false) } else { (hr, true) }));
    }
    let mut out: Vec<Seed> = Vec::new();
    let w = p.w.max(1);
    let mut last: Option<usize> = None;
    if n < w {
        // one window shorter than w: take its minimum
        if let Some((i, &(h, r))) = all.iter().enumerate().filter_map(|(i, x)| x.as_ref().map(|v| (i, v))).min_by_key(|(_, v)| v.0) {
            out.push(Seed { hash: h, pos: i as u32, rev: r });
        }
        return out;
    }
    for s in 0..=(n - w) {
        let mut best: Option<(u64, usize, bool)> = None;
        for i in s..s + w {
            if let Some((h, r)) = all[i] {
                if best.map_or(true, |b| h < b.0) {
                    best = Some((h, i, r));
                }
            }
        }
        if let Some((h, i, r)) = best {
            if last != Some(i) {
                out.push(Seed { hash: h, pos: i as u32, rev: r });
                last = Some(i);
            }
        }
    }
    out
}

/// Bucketed index over a read set: a flat array of (hash, read, pos, rev)
/// sorted by hash, with lookups by binary search.
pub struct Index {
    pub entries: Vec<(u64, u32, u32, bool)>,
    pub params: SeedParams,
    pub per_read: Vec<Vec<Seed>>,
}

impl Index {
    pub fn build(reads: &[Vec<u8>], p: SeedParams) -> Index {
        use rayon::prelude::*;
        let per_read: Vec<Vec<Seed>> = reads.par_iter().map(|r| minimizers(r, &p)).collect();
        let mut entries: Vec<(u64, u32, u32, bool)> = per_read
            .iter()
            .enumerate()
            .flat_map(|(i, ms)| ms.iter().map(move |m| (m.hash, i as u32, m.pos, m.rev)))
            .collect();
        // within a bucket: by canonical strand, then by seed position, so that a
        // query can binary-search for the entries whose diagonal is nearest its own
        entries.par_sort_unstable_by_key(|e| (e.0, e.3, e.2, e.1));
        // drop over-represented buckets
        let mut kept = Vec::with_capacity(entries.len());
        let mut i = 0;
        while i < entries.len() {
            let mut j = i;
            while j < entries.len() && entries[j].0 == entries[i].0 {
                j += 1;
            }
            if j - i <= p.max_occ {
                kept.extend_from_slice(&entries[i..j]);
            }
            i = j;
        }
        Index { entries: kept, params: p, per_read }
    }

    /// Up to `cap` entries of bucket `h` on canonical strand `rev`, chosen as
    /// those whose seed position is nearest `key`.
    ///
    /// This is what keeps candidate generation linear in depth. In a deep
    /// sample every bucket holds hundreds of reads, and pairing a read with
    /// all of them is quadratic in coverage. For fixed-length reads the
    /// diagonal between two reads sharing a seed is a linear function of the
    /// partner's seed position, so the entries nearest in position are the
    /// partners with the largest overlaps -- the ones layout wants -- and the
    /// same partners are nearest in every bucket the pair shares, so their
    /// anchors still accumulate. `cap = 0` disables the cap.
    pub fn nearest(&self, h: u64, rev: bool, key: i64, cap: usize) -> &[(u64, u32, u32, bool)] {
        let b = self.lookup(h);
        let lo = b.partition_point(|e| e.3 < rev);
        let hi = b.partition_point(|e| e.3 <= rev);
        let sub = &b[lo..hi];
        if cap == 0 || sub.len() <= cap {
            return sub;
        }
        let at = sub.partition_point(|e| (e.2 as i64) < key);
        let start = at.saturating_sub(cap / 2).min(sub.len() - cap);
        &sub[start..start + cap]
    }

    pub fn lookup(&self, h: u64) -> &[(u64, u32, u32, bool)] {
        let lo = self.entries.partition_point(|e| e.0 < h);
        let hi = self.entries.partition_point(|e| e.0 <= h);
        &self.entries[lo..hi]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{encode, revcomp};

    #[test]
    fn canonical_across_strands() {
        let s = encode(b"ACGTTGCATGCAAGCTTAGCCGATCGATCGGATCCATGCAAAGTCGTAGCTAGCATCGAT");
        for shape in [SeedShape::Contiguous, SeedShape::CodonSpaced] {
            let p = SeedParams { shape, k: 12, w: 4, max_occ: 100 };
            let f: Vec<u64> = minimizers(&s, &p).iter().map(|m| m.hash).collect();
            let mut r: Vec<u64> = minimizers(&revcomp(&s), &p).iter().map(|m| m.hash).collect();
            let mut f2 = f.clone();
            f2.sort();
            r.sort();
            assert_eq!(f2, r, "shape {:?}", shape);
        }
    }

    #[test]
    fn spaced_seed_ignores_wobble() {
        // two sequences differing only at positions 2,5,8,... (mask phase 0)
        let a = encode(b"ATGGCTAAACTGGGTCGTATTGAAGCTCCGGTTCAGCTG");
        let mut b = a.clone();
        for i in (2..b.len()).step_by(3) {
            b[i] = (b[i] + 1) % 4;
        }
        let p = SeedParams { shape: SeedShape::CodonSpaced, k: 12, w: 1, max_occ: 100 };
        let ha: std::collections::HashSet<u64> = minimizers(&a, &p).iter().map(|m| m.hash).collect();
        let hb: std::collections::HashSet<u64> = minimizers(&b, &p).iter().map(|m| m.hash).collect();
        assert!(ha.intersection(&hb).count() > 3);
        let pc = SeedParams { shape: SeedShape::Contiguous, k: 12, w: 1, max_occ: 100 };
        let ca: std::collections::HashSet<u64> = minimizers(&a, &pc).iter().map(|m| m.hash).collect();
        let cb: std::collections::HashSet<u64> = minimizers(&b, &pc).iter().map(|m| m.hash).collect();
        assert_eq!(ca.intersection(&cb).count(), 0);
    }
}
