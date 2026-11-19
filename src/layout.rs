//! Layout: overlap edges -> read placements on contigs.
//!
//! One viral population should become one contig, so the unit of layout is the
//! connected component of the accepted-overlap graph rather than a unitig. A
//! nucleotide string graph splits at every variant bubble; a component does
//! not, and the variation is carried into the consensus as variant structure
//! instead of as graph fragmentation. The cost is that a single wrong edge
//! merges two entities, which is why every contiguity number is reported next
//! to a chimaera number.
//!
//! Each component is embedded by a maximum spanning tree (Kruskal on edge
//! score): reads get an orientation and a leftmost coordinate by walking the
//! tree from its highest-degree read. Offsets accumulate small indel errors
//! along tree paths; consensus re-aligns every read to the draft, so the
//! embedding only has to be right to within the polishing window.

use crate::overlap::Edge;
use std::collections::VecDeque;

pub struct Dsu {
    p: Vec<u32>,
    r: Vec<u8>,
}

impl Dsu {
    pub fn new(n: usize) -> Dsu {
        Dsu { p: (0..n as u32).collect(), r: vec![0; n] }
    }
    pub fn find(&mut self, mut x: u32) -> u32 {
        while self.p[x as usize] != x {
            let gp = self.p[self.p[x as usize] as usize];
            self.p[x as usize] = gp;
            x = gp;
        }
        x
    }
    pub fn union(&mut self, a: u32, b: u32) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        let (ra, rb) = if self.r[ra as usize] < self.r[rb as usize] { (rb, ra) } else { (ra, rb) };
        self.p[rb as usize] = ra;
        if self.r[ra as usize] == self.r[rb as usize] {
            self.r[ra as usize] += 1;
        }
        true
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Placement {
    pub read: u32,
    /// leftmost coordinate of the oriented read on the component's axis
    pub x: i64,
    pub rev: bool,
}

/// Keep only edges that rank in the top `k` (by per-column score) for at least
/// one endpoint -- or for both, when `mutual`. In a swarm every read has many
/// same-population partners that score better than any cross-population
/// partner, so a cross-population edge survives only where a read genuinely
/// has nothing better. This is a best-overlap graph in the miniasm/Shasta
/// sense, and it is the main defence against one weak edge merging two
/// entities. `k = 0` disables the filter.
pub fn best_k(n: usize, edges: &[Edge], k: usize, mutual: bool) -> Vec<Edge> {
    if k == 0 {
        return edges.to_vec();
    }
    let mut per: Vec<Vec<(f32, usize)>> = vec![Vec::new(); n];
    for (i, e) in edges.iter().enumerate() {
        per[e.a as usize].push((e.per_col, i));
        per[e.b as usize].push((e.per_col, i));
    }
    let mut votes = vec![0u8; edges.len()];
    for v in per.iter_mut() {
        v.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
        for &(_, i) in v.iter().take(k) {
            votes[i] += 1;
        }
    }
    let need = if mutual { 2 } else { 1 };
    edges.iter().zip(votes).filter(|(_, v)| *v >= need).map(|(e, _)| *e).collect()
}

/// Returns placements grouped by component (largest first). Singletons are
/// returned as their own components so that nothing is silently dropped.
pub fn layout(lens: &[usize], edges: &[Edge]) -> Vec<Vec<Placement>> {
    let n = lens.len();
    let mut order: Vec<usize> = (0..edges.len()).collect();
    order.sort_by(|&x, &y| edges[y].score.partial_cmp(&edges[x].score).unwrap());
    let mut dsu = Dsu::new(n);
    let mut adj: Vec<Vec<(u32, bool, i64, bool)>> = vec![Vec::new(); n]; // (nbr, rev, offset, i_is_a)
    for &e in &order {
        let ed = &edges[e];
        if dsu.union(ed.a, ed.b) {
            adj[ed.a as usize].push((ed.b, ed.rev, ed.offset, true));
            adj[ed.b as usize].push((ed.a, ed.rev, ed.offset, false));
        }
    }
    let mut comps: std::collections::HashMap<u32, Vec<u32>> = std::collections::HashMap::new();
    for i in 0..n as u32 {
        let r = dsu.find(i);
        comps.entry(r).or_default().push(i);
    }
    let mut out = Vec::new();
    let mut placed: Vec<Option<(i64, bool)>> = vec![None; n];
    for (_, members) in comps {
        let root = *members.iter().max_by_key(|&&m| adj[m as usize].len()).unwrap();
        placed[root as usize] = Some((0, false));
        let mut q = VecDeque::from([root]);
        let mut comp = Vec::with_capacity(members.len());
        while let Some(u) = q.pop_front() {
            let (xu, ou) = placed[u as usize].unwrap();
            comp.push(Placement { read: u, x: xu, rev: ou });
            for &(v, rel, off, u_is_a) in &adj[u as usize] {
                if placed[v as usize].is_some() {
                    continue;
                }
                let lu = lens[u as usize] as i64;
                let lv = lens[v as usize] as i64;
                // offset of v relative to u in u's forward frame
                let off_uv = if u_is_a { off } else if rel { lu - lv + off } else { -off };
                let ov = ou ^ rel;
                let xv = if !ou { xu + off_uv } else { xu + lu - off_uv - lv };
                placed[v as usize] = Some((xv, ov));
                q.push_back(v);
            }
        }
        let minx = comp.iter().map(|p| p.x).min().unwrap();
        for p in comp.iter_mut() {
            p.x -= minx;
        }
        comp.sort_by_key(|p| p.x);
        out.push(comp);
    }
    out.sort_by(|a, b| b.len().cmp(&a.len()));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::Edge;
    fn e(a: u32, b: u32, rev: bool, offset: i64) -> Edge {
        Edge { a, b, rev, offset, score: 100.0, per_col: 1.0, identity: 1.0, tier: 1 }
    }
    #[test]
    fn embeds_mixed_orientations() {
        // genome positions: r0 at 0 (+), r1 at 50 (-), r2 at 100 (+); all len 150
        // r0->r1: r1 reversed starts at 50 in r0 coords
        // r1->r2: in r1 forward frame, rc(r2) starts at 150-150+... compute via flip rule
        let lens = vec![150, 150, 150];
        // r0 (a) to r2 (b), same orientation, offset 100
        let edges = vec![e(0, 1, true, 50), e(0, 2, false, 100)];
        let comps = layout(&lens, &edges);
        assert_eq!(comps.len(), 1);
        let mut c = comps[0].clone();
        c.sort_by_key(|p| p.read);
        let r0 = c[0];
        // relative geometry is what matters
        let rel = |p: Placement| if r0.rev { -(p.x - r0.x) } else { p.x - r0.x };
        assert_eq!(rel(c[2]).abs(), 100);
        assert_eq!(rel(c[1]).abs(), 50);
        assert_eq!(c[1].rev != r0.rev, true);
        assert_eq!(c[2].rev == r0.rev, true);
    }
}
