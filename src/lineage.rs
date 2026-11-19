//! Splitting a contig into lineages by alleles linked within reads.
//!
//! Why this exists. A component is one population by design (layout.rs), but
//! real samples hold related lineages that are different genomes. In Berlin
//! wastewater (results/ww_berlin_consensus) one hapax contig carried reads of at
//! least five astrovirus lineages 10-14% apart, and the per-column vote then
//! followed whichever lineage was locally commoner: 107 of 151 residual
//! substitutions against the authors' contig were another lineage's base.
//!
//! What was tried before, and why this differs. `phase.rs` fits two
//! populations to all sites at once by EM; with several lineages present over
//! different stretches of the contig a global two-way fit settles on mosaics
//! (on Berlin it took identity from 97.1% to 90.8%). `layout::split_sparse`
//! cuts the overlap graph, and the cheapest cut of a linear graph runs across
//! the genome rather than between lineages. Both fail for want of position:
//! lineages are local objects that have to be followed along the contig.
//!
//! The method, in two stages, both position-anchored:
//!
//! 1. **Agglomerate reads.** Informative sites are columns where a second
//!    allele is carried by at least `min_alt` reads and `min_af` of the depth;
//!    sequencing error almost never reaches that at tens-fold depth. Pairs of
//!    overlapping reads that agree at `min_shared` or more shared sites and
//!    conflict at none are merged strongest-first, but a merge between two
//!    groups is refused if their allele profiles (per-site majority) conflict at
//!    more than `max_div` per base of the span they share. Profiles carry the
//!    whole history of a group, so a pair of reads that happens to agree at
//!    sites where only a third lineage differs cannot chain two lineages.
//! 2. **Stitch groups along the contig.** Groups of at least `min_group` reads
//!    are chained left to right by mutual best compatible neighbour (same
//!    profile test; gaps up to `max_gap` with no conflicting site), and a group
//!    contained in a compatible chain joins it. Mutual-best keeps an
//!    ambiguous junction unjoined rather than guessing.
//!
//! Remaining reads go to the chain whose profile they conflict with least,
//! ties to the chain with more reads. Chains below `min_lineage_reads` are
//! folded back the same way. More than one chain left means a split.
//!
//! Calibration of `max_div` = 0.04 per base: Berlin lineages differ at
//! ~0.10-0.14 per base; within a lineage, reads differ at ~0.016 (polymorphic
//! sites from mpileup) plus error. A within-lineage haplotype split therefore
//! conflicts at well under 0.04, a between-lineage one well over. The
//! simulator's within-population diversity (0.05-0.20 between tips) straddles
//! this, which is why the ladder is the check on over-splitting.

use crate::consensus::Contig;
use crate::layout::Placement;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct LineageParams {
    pub min_alt: u32,
    pub min_af: f32,
    pub min_shared: usize,
    pub max_div: f32,
    /// ...or conflicting at no more than this fraction of shared sites. Between
    /// Berlin lineages the fraction is ~0.6; within one, a short junction can show
    /// 3 of 12 from its edge reads alone.
    pub max_conflict_frac: f32,
    pub min_group: usize,
    pub max_gap: i64,
    /// widest stretch without informative sites that a unique join may cross
    pub stitch_gap: i64,
    pub min_lineage_reads: usize,
    /// lineages smaller than this fraction of the contig's reads are folded back
    pub min_lineage_frac: f32,
    /// rounds of split-and-rebuild (main.rs)
    pub max_depth: usize,
    /// Chains whose consensus differ by less than this per base, where they
    /// overlap, are one lineage and are merged back. This is the knob that
    /// decides what a separate genome is; see `split`.
    pub min_split_div: f32,
}

impl Default for LineageParams {
    fn default() -> Self {
        LineageParams { min_alt: 3, min_af: 0.10, min_shared: 2, max_div: 0.04, max_conflict_frac: 0.3, min_group: 15, max_gap: 100, stitch_gap: 1000, min_lineage_reads: 30, min_lineage_frac: 0.02, max_depth: 1, min_split_div: 0.0 }
    }
}

/// Per-site allele counts plus the coordinate span they were seen over.
#[derive(Clone, Default)]
struct Profile {
    counts: HashMap<u32, [u32; 4]>,
    lo: i64,
    hi: i64,
    n: usize,
}

fn maj(c: &[u32; 4]) -> u8 {
    let mut b = 0;
    for i in 1..4 {
        if c[i] > c[b] {
            b = i;
        }
    }
    b as u8
}

impl Profile {
    fn add_read(&mut self, v: &[(u32, u8)], lo: i64, hi: i64) {
        for &(s, b) in v {
            self.counts.entry(s).or_insert([0; 4])[b as usize] += 1;
        }
        if self.n == 0 {
            self.lo = lo;
            self.hi = hi;
        } else {
            self.lo = self.lo.min(lo);
            self.hi = self.hi.max(hi);
        }
        self.n += 1;
    }
    fn absorb(&mut self, o: Profile) {
        for (s, c) in o.counts {
            let e = self.counts.entry(s).or_insert([0; 4]);
            for i in 0..4 {
                e[i] += c[i];
            }
        }
        if self.n == 0 {
            self.lo = o.lo;
            self.hi = o.hi;
        } else {
            self.lo = self.lo.min(o.lo);
            self.hi = self.hi.max(o.hi);
        }
        self.n += o.n;
    }
    /// (shared sites, conflicting sites) by per-site majority. A site counts only
    /// where both sides have two or more reads, or is a single read on either
    /// side, so one read's error at a group edge cannot veto a join.
    fn compare(&self, o: &Profile) -> (usize, usize) {
        let (a, b) = if self.counts.len() <= o.counts.len() { (self, o) } else { (o, self) };
        let (mut sh, mut cf) = (0, 0);
        let single = a.n == 1 || b.n == 1;
        for (s, c) in &a.counts {
            if let Some(d) = b.counts.get(s) {
                if !single && (c.iter().sum::<u32>() < 2 || d.iter().sum::<u32>() < 2) {
                    continue;
                }
                sh += 1;
                if maj(c) != maj(d) {
                    cf += 1;
                }
            }
        }
        (sh, cf)
    }
    /// (shared sites, conflicting sites) of group `g` against this lineage. A
    /// site counts where both have two or more reads (or `g` is one read); it
    /// conflicts when `g`'s majority allele is carried by under a fifth of the
    /// lineage's reads. Haplotypes already in the lineage keep their alleles, so
    /// a sister haplotype is not a conflict, while another lineage conflicts at
    /// every site that tells it apart.
    fn fit_group(&self, g: &Profile) -> (usize, usize) {
        let (mut sh, mut cf) = (0, 0);
        for (s, c) in &g.counts {
            if let Some(d) = self.counts.get(s) {
                let dn: u32 = d.iter().sum();
                if dn < 2 || (g.n > 1 && c.iter().sum::<u32>() < 2) {
                    continue;
                }
                sh += 1;
                if 5 * d[maj(c) as usize] < dn {
                    cf += 1;
                }
            }
        }
        (sh, cf)
    }
    fn overlap(&self, o: &Profile) -> i64 {
        self.hi.min(o.hi) - self.lo.max(o.lo)
    }
    /// sites a read shares with this profile, and how many it conflicts at. An
    /// allele carried by a fifth or more of the profile's reads is not a
    /// conflict: haplotypes inside one lineage keep their own alleles.
    fn read_fit(&self, v: &[(u32, u8)]) -> (usize, usize) {
        let (mut sh, mut cf) = (0, 0);
        for &(s, b) in v {
            if let Some(c) = self.counts.get(&s) {
                sh += 1;
                if 5 * c[b as usize] < c.iter().sum::<u32>() {
                    cf += 1;
                }
            }
        }
        (sh, cf)
    }
}

struct Dsu(Vec<usize>);
impl Dsu {
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
}

/// Read sets, one per lineage, as placements shifted so each set starts at 0.
/// Returns an empty vector when the contig is one lineage.
    // Stitch pieces end to end across a stretch with no informative site, but
    // only where the join is unique both ways: X has exactly one piece that can
    // continue it to the right, and that piece has X as its only left
    // neighbour. Two candidates means two lineages meet there and no read says
    // which continues which, so the pieces stay apart (short and right).
fn stitch(lprof: &mut Vec<Profile>, lineages: &mut Vec<Vec<usize>>, p: &LineageParams) {
    loop {
        let k = lprof.len();
        let cont = |x: usize, y: usize| {
            let (a, b) = (&lprof[x], &lprof[y]);
            let gap = b.lo - a.hi;
            if b.lo <= a.lo || b.hi <= a.hi || gap < -p.max_gap || gap > p.stitch_gap {
                return false;
            }
            let (sh, cf) = a.fit_group(b);
            sh == 0 || cf as f32 <= p.max_div * a.overlap(b).max(1) as f32 || cf as f32 <= p.max_conflict_frac * sh as f32
        };
        // first, pieces that overlap and fit each other on shared sites (in both
        // directions) are one lineage by evidence: best evidence first
        let mut join = None;
        let mut best_ev = 0i64;
        for x in 0..k {
            for y in x + 1..k {
                let ov = lprof[x].overlap(&lprof[y]);
                if ov <= 0 {
                    continue;
                }
                let fits = |a: &Profile, b: &Profile| {
                    let (sh, cf) = a.fit_group(b);
                    (sh >= p.min_shared && (cf as f32 <= p.max_div * ov as f32 || cf as f32 <= p.max_conflict_frac * sh as f32)).then_some(sh as i64 - 2 * cf as i64)
                };
                if let (Some(e1), Some(e2)) = (fits(&lprof[x], &lprof[y]), fits(&lprof[y], &lprof[x])) {
                    if join.is_none() || e1.min(e2) > best_ev {
                        best_ev = e1.min(e2);
                        join = Some((x, y));
                    }
                }
            }
        }
        // only the nearest continuations count (within max_gap of the nearest):
        // a piece further off is not a rival for this junction
        let right: Vec<Vec<usize>> = (0..k)
            .map(|x| {
                let c: Vec<usize> = (0..k).filter(|&y| y != x && cont(x, y)).collect();
                let near = c.iter().map(|&y| lprof[y].lo).min().unwrap_or(0);
                c.into_iter().filter(|&y| lprof[y].lo <= near + p.max_gap).collect()
            })
            .collect();
        for x in 0..k {
            if join.is_some() {
                break;
            }
            if let [y] = right[x][..] {
                if (0..k).filter(|&z| right[z].contains(&y)).count() == 1 {
                    join = Some((x, y));
                    break;
                }
            }
        }
        match join {
            Some((x, y)) => {
                let (py, my) = (lprof.remove(y), lineages.remove(y));
                let x = if x > y { x - 1 } else { x };
                lprof[x].absorb(py);
                lineages[x].extend(my);
            }
            None => break,
        }
    }
}

pub fn split(c: &Contig, p: &LineageParams) -> Vec<Vec<Placement>> {
    let n = c.alns.len();
    if n < 2 * p.min_lineage_reads {
        return Vec::new();
    }
    // informative sites
    let mut site_of = vec![u32::MAX; c.columns.len()];
    let mut ns = 0u32;
    for (i, col) in c.columns.iter().enumerate() {
        let mut v = [col.counts[0], col.counts[1], col.counts[2], col.counts[3]];
        let d: u32 = v.iter().sum();
        v.sort_unstable();
        if d > 0 && v[2] >= p.min_alt && v[2] as f32 >= p.min_af * d as f32 {
            site_of[i] = ns;
            ns += 1;
        }
    }
    if ns < 2 {
        return Vec::new();
    }
    let vecs: Vec<Vec<(u32, u8)>> = c
        .alns
        .iter()
        .map(|a| {
            a.cols
                .iter()
                .enumerate()
                .filter_map(|(k, &b)| {
                    let s = site_of[a.ref_start + k];
                    (s != u32::MAX && b < 4).then_some((s, b))
                })
                .collect()
        })
        .collect();
    let span = |i: usize| (c.alns[i].ref_start as i64, c.alns[i].ref_end as i64);

    // 1. pairs of overlapping reads, full agreement at >= min_shared sites
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| c.alns[i].ref_start);
    let mut pairs: Vec<(usize, usize, usize)> = Vec::new();
    for (ii, &i) in order.iter().enumerate() {
        let (_, ei) = span(i);
        for &j in &order[ii + 1..] {
            if span(j).0 >= ei {
                break;
            }
            let (va, vb) = (&vecs[i], &vecs[j]);
            let (mut x, mut y, mut agree, mut conflict) = (0, 0, 0, 0);
            while x < va.len() && y < vb.len() {
                if va[x].0 == vb[y].0 {
                    if va[x].1 == vb[y].1 {
                        agree += 1;
                    } else {
                        conflict += 1;
                        break;
                    }
                    x += 1;
                    y += 1;
                } else if va[x].0 < vb[y].0 {
                    x += 1;
                } else {
                    y += 1;
                }
            }
            if conflict == 0 && agree >= p.min_shared {
                pairs.push((agree, i, j));
            }
        }
    }
    pairs.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    let mut dsu = Dsu((0..n).collect());
    let mut prof: Vec<Option<Profile>> = (0..n)
        .map(|i| {
            let mut pr = Profile::default();
            let (lo, hi) = span(i);
            pr.add_read(&vecs[i], lo, hi);
            Some(pr)
        })
        .collect();
    let compatible = |a: &Profile, b: &Profile, p: &LineageParams| {
        let (sh, cf) = a.compare(b);
        cf as f32 <= p.max_div * a.overlap(b).max(1) as f32 || cf as f32 <= p.max_conflict_frac * sh as f32
    };
    for &(_, i, j) in &pairs {
        let (x, y) = (dsu.find(i), dsu.find(j));
        if x == y {
            continue;
        }
        if !compatible(prof[x].as_ref().unwrap(), prof[y].as_ref().unwrap(), p) {
            continue;
        }
        let (big, small) = if prof[x].as_ref().unwrap().n >= prof[y].as_ref().unwrap().n { (x, y) } else { (y, x) };
        let s = prof[small].take().unwrap();
        prof[big].as_mut().unwrap().absorb(s);
        dsu.0[small] = big;
    }

    // 2. stitch groups into chains
    let mut groups: Vec<(Profile, Vec<usize>)> = Vec::new();
    {
        let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..n {
            let r = dsu.find(i);
            members.entry(r).or_default().push(i);
        }
        let mut roots: Vec<usize> = members.keys().copied().collect();
        roots.sort_unstable();
        for r in roots {
            let m = members.remove(&r).unwrap();
            if m.len() >= p.min_group {
                groups.push((prof[r].take().unwrap(), m));
            }
        }
    }
    if groups.len() < 2 {
        return Vec::new();
    }
    // Groups become lineages by colouring: largest first, each group joins the
    // lineage whose overlapping groups it is compatible with all of, preferring
    // the most allele evidence and then the bigger lineage; a group compatible
    // with none starts a new lineage. Two groups that never overlap are never
    // forced apart, so a stretch with no informative site (where lineages are
    // indistinguishable) does not create one: it joins the dominant lineage.
    groups.sort_by(|a, b| b.0.n.cmp(&a.0.n).then(a.0.lo.cmp(&b.0.lo)));
    let mut lineages: Vec<Vec<usize>> = Vec::new();
    let mut lprof: Vec<Profile> = Vec::new();
    // (evidence, lineages it shares sites with but does not fit)
    let judge = |lprof: &[Profile], pg: &Profile| -> (Option<(i64, usize)>, Vec<usize>) {
        let (mut best, mut bad): (Option<(i64, usize)>, Vec<usize>) = (None, Vec::new());
        for (li, lp) in lprof.iter().enumerate() {
            let ov = pg.overlap(lp);
            if ov <= 0 {
                continue;
            }
            let (sh, cf) = lp.fit_group(pg);
            if sh == 0 {
                continue;
            }
            if !(cf as f32 <= p.max_div * ov as f32 || cf as f32 <= p.max_conflict_frac * sh as f32) {
                bad.push(li);
                continue;
            }
            let ev = sh as i64 - 2 * cf as i64;
            if best.map_or(true, |b| ev > b.0) {
                best = Some((ev, li));
            }
        }
        (best, bad)
    };
    let mut silent = Vec::new();
    for gi in 0..groups.len() {
        let (best, bad) = judge(&lprof, &groups[gi].0);
        let g = groups[gi].0.clone();
        match best {
            Some((_, li)) => {
                lineages[li].push(gi);
                lprof[li].absorb(g);
            }
            None if bad.is_empty() && !lprof.is_empty() => silent.push(gi),
            None => {
                lineages.push(vec![gi]);
                lprof.push(g);
            }
        }
    }
    // A group that shares no site with any lineage is no evidence of a new one.
    // Once the lineages are laid out by evidence, it joins the biggest lineage
    // whose span it lies in that it was not shown to differ from; failing that
    // it starts a piece of its own, for the stitching below.
    for gi in silent {
        let (best, bad) = judge(&lprof, &groups[gi].0);
        let pg = &groups[gi].0;
        let target = best.map(|b| b.1).or_else(|| {
            (0..lprof.len())
                .filter(|li| !bad.contains(li))
                .map(|li| ((-pg.overlap(&lprof[li])).max(0), usize::MAX - lprof[li].n, li))
                .filter(|x| x.0 == 0)
                .min()
                .map(|x| x.2)
        });
        let g = pg.clone();
        match target {
            Some(li) => {
                lineages[li].push(gi);
                lprof[li].absorb(g);
            }
            None => {
                lineages.push(vec![gi]);
                lprof.push(g);
            }
        }
    }
    stitch(&mut lprof, &mut lineages, p);
    let mut merged_groups: Vec<(Profile, Vec<usize>)> = Vec::new();
    let mut taken: Vec<Option<(Profile, Vec<usize>)>> = groups.into_iter().map(Some).collect();
    for members in lineages {
        let mut pr = Profile::default();
        let mut m = Vec::new();
        for h in members {
            let (ph, mh) = taken[h].take().unwrap();
            pr.absorb(ph);
            m.extend(mh);
        }
        merged_groups.push((pr, m));
    }
    let mut groups = merged_groups;

    // 2b. chains closer than min_split_div are one lineage
    loop {
        let g = groups.len();
        let mut best: Option<(f32, usize, usize)> = None;
        for a in 0..g {
            for b in a + 1..g {
                let ov = groups[a].0.overlap(&groups[b].0);
                if ov < p.min_group as i64 {
                    continue;
                }
                let (_, cf) = groups[a].0.compare(&groups[b].0);
                let d = cf as f32 / ov as f32;
                if best.map_or(true, |x| d < x.0) {
                    best = Some((d, a, b));
                }
            }
        }
        match best {
            Some((d, a, b)) if d < p.min_split_div => {
                let (pb, mb) = groups.remove(b);
                groups[a].0.absorb(pb);
                groups[a].1.extend(mb);
            }
            _ => break,
        }
    }
    if groups.len() < 2 {
        return Vec::new();
    }

    // 3. every read to its best chain; small chains folded back
    // the nearer lineage first (0 when the read lies inside its span: a lineage
    // the read never meets has no sites to conflict at, which is no evidence),
    // then fewer conflicts, then the lineage the read lies deeper inside (so
    // two lineages meeting at a site-free stretch each grow from their own
    // side rather than the bigger one taking the stretch), more shared sites,
    // the bigger lineage
    let assign = |groups: &[(Profile, Vec<usize>)], i: usize| -> (usize, i64) {
        let (lo, hi) = span(i);
        let mut best: Option<((i64, usize, i64, usize, usize), usize)> = None;
        for (k, (pr, _)) in groups.iter().enumerate() {
            let (sh, cf) = pr.read_fit(&vecs[i]);
            let dist = (pr.lo - hi).max(lo - pr.hi).max(0);
            let inside = pr.hi.min(hi) - pr.lo.max(lo);
            let key = (dist, cf, -inside, usize::MAX - sh, usize::MAX - pr.n);
            if best.map_or(true, |b| key < b.0) {
                best = Some((key, k));
            }
        }
        let b = best.unwrap();
        (b.1, b.0 .0)
    };
    // unowned reads, nearest first, each growing its lineage's span, so a lineage
    // extends outward through a stretch without informative sites read by read
    let place = |groups: &mut Vec<(Profile, Vec<usize>)>, mut todo: Vec<usize>| {
        while !todo.is_empty() {
            let (k, _) = todo.iter().enumerate().map(|(k, &i)| (k, assign(groups, i).1)).min_by_key(|x| x.1).unwrap();
            let i = todo.swap_remove(k);
            let (t, _) = assign(groups, i);
            let (lo, hi) = span(i);
            groups[t].0.add_read(&vecs[i], lo, hi);
            groups[t].1.push(i);
        }
    };
    let min_reads = p.min_lineage_reads.max((p.min_lineage_frac * n as f32) as usize);
    loop {
        let small = (0..groups.len()).filter(|&k| groups[k].1.len() < min_reads).min_by_key(|&k| groups[k].1.len());
        match small {
            Some(k) if groups.len() > 1 => {
                let (_, m) = groups.remove(k);
                place(&mut groups, m);
            }
            _ => break,
        }
    }
    if groups.len() < 2 {
        return Vec::new();
    }
    let mut owner = vec![usize::MAX; n];
    for (k, (_, m)) in groups.iter().enumerate() {
        for &i in m {
            owner[i] = k;
        }
    }
    let todo: Vec<usize> = (0..n).filter(|&i| owner[i] == usize::MAX).collect();
    place(&mut groups, todo);
    // the reads placed through site-free stretches can make a junction unique
    // that was not before (each piece has grown towards its own continuation)
    {
        let (mut pr, mut mm): (Vec<Profile>, Vec<Vec<usize>>) = groups.drain(..).unzip();
        stitch(&mut pr, &mut mm, p);
        groups = pr.into_iter().zip(mm).collect();
    }
    if groups.len() < 2 {
        return Vec::new();
    }
    for (k, (_, m)) in groups.iter().enumerate() {
        for &i in m {
            owner[i] = k;
        }
    }
    let mut out: Vec<Vec<Placement>> = vec![Vec::new(); groups.len()];
    for i in 0..n {
        let a = &c.alns[i];
        out[owner[i]].push(Placement { read: a.read, rev: a.rev, x: a.ref_start as i64 });
    }
    for set in &mut out {
        let m = set.iter().map(|q| q.x).min().unwrap_or(0);
        for q in set.iter_mut() {
            q.x -= m;
        }
    }
    out.retain(|s| s.len() >= 2);
    if out.len() < 2 {
        return Vec::new();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::{build, ConsensusParams};

    struct Rng(u64);
    impl Rng {
        fn f(&mut self) -> f64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn mutate(s: &[u8], rate: f64, rng: &mut Rng) -> Vec<u8> {
        s.iter().map(|&b| if rng.f() < rate { (b + 1 + (rng.f() * 3.0) as u8) % 4 } else { b }).collect()
    }

    /// Lineage A runs the whole contig as two haplotypes 1.5% apart; lineage B,
    /// 12% from A, covers only the second half (as NODE_1338 does on Berlin).
    /// Reads carry 0.3% error. Returns reads, placements and each read's lineage.
    fn world(with_b: bool) -> (Vec<Vec<u8>>, Vec<Placement>, Vec<bool>) {
        let mut rng = Rng(11);
        let n = 3000usize;
        let a1: Vec<u8> = (0..n).map(|_| (rng.f() * 4.0) as u8).collect();
        let a2 = mutate(&a1, 0.015, &mut rng);
        let b = mutate(&a1, 0.12, &mut rng);
        let (mut reads, mut pl, mut is_b) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..4000 {
            let from_b = with_b && i % 3 == 0;
            let (src, lo) = if from_b { (&b, n / 2) } else if i % 2 == 0 { (&a1, 0) } else { (&a2, 0) };
            let start = lo + ((rng.f() * (n - lo - 100) as f64) as usize).min(n - lo - 100);
            reads.push(mutate(&src[start..start + 100], 0.003, &mut rng));
            pl.push(Placement { read: i as u32, rev: false, x: start as i64 });
            is_b.push(from_b);
        }
        (reads, pl, is_b)
    }

    #[test]
    fn splits_a_lineage_off_and_keeps_haplotypes_together() {
        let cp = ConsensusParams::default();
        let p = LineageParams::default();
        let (reads, pl, is_b) = world(true);
        let c = build(&reads, &pl, &cp);
        let sets = split(&c, &p);
        let sizes: Vec<(usize, usize)> = sets.iter().map(|s| (s.len(), s.iter().filter(|q| is_b[q.read as usize]).count())).collect();
        assert_eq!(sets.len(), 2, "lineage B should split off as one set; got (reads, from B) {:?}", sizes);
        for s in &sets {
            let nb = s.iter().filter(|q| is_b[q.read as usize]).count();
            let purity = nb.max(s.len() - nb) as f64 / s.len() as f64;
            assert!(purity > 0.95, "set of {} reads is {:.2} pure", s.len(), purity);
        }

        let (reads1, pl1, _) = world(false);
        let c1 = build(&reads1, &pl1, &cp);
        assert!(split(&c1, &p).is_empty(), "two haplotypes 1.5% apart are one lineage");
    }
}
