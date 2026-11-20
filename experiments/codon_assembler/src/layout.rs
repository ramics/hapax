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

/// Group reads into entities in two stages.
///
/// Stage 1 takes connected components of the best-k graph. That graph is
/// precise but, at high per-variant coverage, it tends to phase a population
/// apart: every read's top partners are reads of its own haplotype, so
/// cross-variant edges never make anyone's top k and one population splits
/// into several haplotype clusters. That is the opposite of what is wanted.
///
/// Stage 2 therefore merges stage-1 clusters back together when the *full*
/// accepted-edge set links them by many edges: at least `min_bridge`, and at
/// least `bridge_frac` of the smaller cluster's reads. Two haplotypes of one
/// swarm share hundreds of accepted overlaps; two sister populations meet in a
/// handful of conserved stretches. Counting bridges is what tells them apart,
/// where no single overlap score can.
///
/// `stage1`, when given, replaces the best-k clusters (the tier-1 components
/// from dispatch are the natural choice: near-identical, indel-free joins).
pub fn group(n: usize, edges: &[Edge], k: usize, min_bridge: usize, bridge_frac: f32, stage1: Option<&[u32]>) -> Vec<u32> {
    let c1: Vec<u32> = match stage1 {
        Some(s) => s.to_vec(),
        None => {
            let strong = best_k(n, edges, k, false);
            let mut d1 = Dsu::new(n);
            for e in &strong {
                d1.union(e.a, e.b);
            }
            (0..n as u32).map(|i| d1.find(i)).collect()
        }
    };
    let mut size = std::collections::HashMap::<u32, usize>::new();
    for &c in &c1 {
        *size.entry(c).or_insert(0) += 1;
    }
    let mut bridges = std::collections::HashMap::<(u32, u32), usize>::new();
    for e in edges {
        let (x, y) = (c1[e.a as usize], c1[e.b as usize]);
        if x != y {
            *bridges.entry((x.min(y), x.max(y))).or_insert(0) += 1;
        }
    }
    let mut b: Vec<((u32, u32), usize)> = bridges.into_iter().collect();
    b.sort_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(&y.0)));
    let mut d2 = Dsu::new(n);
    for ((x, y), cnt) in b {
        let small = size[&x].min(size[&y]);
        if cnt >= min_bridge && cnt as f32 >= bridge_frac * small as f32 {
            d2.union(x, y);
        }
    }
    (0..n).map(|i| d2.find(c1[i])).collect()
}

/// Returns placements grouped by component (largest first). Singletons are
/// returned as their own components so that nothing is silently dropped.
pub fn layout(lens: &[usize], edges: &[Edge], groups: Option<&[u32]>) -> Vec<Vec<Placement>> {
    let n = lens.len();
    let mut order: Vec<usize> = (0..edges.len()).collect();
    order.sort_by(|&x, &y| edges[y].score.partial_cmp(&edges[x].score).unwrap());
    let mut dsu = Dsu::new(n);
    let mut adj: Vec<Vec<(u32, bool, i64, bool)>> = vec![Vec::new(); n]; // (nbr, rev, offset, i_is_a)
    for &e in &order {
        let ed = &edges[e];
        if let Some(g) = groups {
            if g[ed.a as usize] != g[ed.b as usize] {
                continue;
            }
        }
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
        let comps = layout(&lens, &edges, None);
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

/// Split a read group at a sparse cut, when it has one.
///
/// Two swarms of the same virus from different hosts end up in one group
/// because they share overlaps through conserved windows that no per-pair
/// score can refuse (see src/phase.rs for the measurements, and
/// experiments/sister_pairs.py). But those shared overlaps are *few*: on the
/// real HIV *env* pair the overlaps between the two hosts' reads are 1.3% of
/// the overlaps within the sparser host, and on the HBV pair 3.4%. Two
/// populations are therefore two dense clusters joined by a thin cut, while
/// one swarm is a single dense cluster: the haplotypes inside it connect
/// through hundreds of overlaps each. So look for a cut and take it only if it
/// is thin.
///
/// The cut comes from two-label propagation seeded at the ends of a BFS
/// diameter path: each read repeatedly takes the majority label of its
/// neighbours, which settles on the natural boundary when there is one. It is
/// taken only when the crossing edges are below `max_cut` of the smaller
/// side's internal edges and both sides hold at least `min_side_frac` of the
/// reads. Recursion allows up to `2^depth` populations per group.
pub fn split_sparse(n: usize, groups: &[u32], edges: &[Edge], max_cut: f32, min_side_frac: f32, depth: usize, min_reads: usize) -> Vec<u32> {
    let mut out = groups.to_vec();
    let mut members: std::collections::HashMap<u32, Vec<u32>> = std::collections::HashMap::new();
    for i in 0..n as u32 {
        members.entry(groups[i as usize]).or_default().push(i);
    }
    let mut by_group: std::collections::HashMap<u32, Vec<&Edge>> = std::collections::HashMap::new();
    for e in edges {
        if groups[e.a as usize] == groups[e.b as usize] {
            by_group.entry(groups[e.a as usize]).or_default().push(e);
        }
    }
    let mut next_label = groups.iter().copied().max().unwrap_or(0) + 1;
    for (g, mem) in members.iter() {
        if mem.len() < min_reads {
            continue;
        }
        let ge = match by_group.get(g) {
            Some(v) => v.clone(),
            None => continue,
        };
        let mut queue: Vec<(Vec<u32>, Vec<&Edge>, usize)> = vec![(mem.clone(), ge, 0)];
        while let Some((mem, ge, d)) = queue.pop() {
            if d >= depth || mem.len() < min_reads {
                continue;
            }
            match two_cut(&mem, &ge, max_cut, min_side_frac) {
                None => {}
                Some(side) => {
                    let (mut ma, mut mb) = (Vec::new(), Vec::new());
                    let lab = next_label;
                    next_label += 1;
                    for (&r, &s) in mem.iter().zip(side.iter()) {
                        if s {
                            out[r as usize] = lab;
                            mb.push(r);
                        } else {
                            ma.push(r);
                        }
                    }
                    let ea: Vec<&Edge> = ge.iter().copied().filter(|e| out[e.a as usize] == out[e.b as usize] && !side_of(&mem, &side, e.a)).collect();
                    let eb: Vec<&Edge> = ge.iter().copied().filter(|e| out[e.a as usize] == out[e.b as usize] && side_of(&mem, &side, e.a)).collect();
                    queue.push((ma, ea, d + 1));
                    queue.push((mb, eb, d + 1));
                }
            }
        }
    }
    out
}

fn side_of(mem: &[u32], side: &[bool], r: u32) -> bool {
    match mem.binary_search(&r) {
        Ok(i) => side[i],
        Err(_) => false,
    }
}

/// Label propagation from the two ends of a BFS diameter path. Returns the
/// bipartition if its cut is thin enough to be two populations.
fn two_cut(mem: &[u32], edges: &[&Edge], max_cut: f32, min_side_frac: f32) -> Option<Vec<bool>> {
    let idx: std::collections::HashMap<u32, usize> = mem.iter().enumerate().map(|(i, &r)| (r, i)).collect();
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); mem.len()];
    for e in edges {
        if let (Some(&a), Some(&b)) = (idx.get(&e.a), idx.get(&e.b)) {
            adj[a].push(b as u32);
            adj[b].push(a as u32);
        }
    }
    let far = |start: usize, adj: &Vec<Vec<u32>>| -> usize {
        let mut seen = vec![false; adj.len()];
        let mut q = VecDeque::new();
        q.push_back(start);
        seen[start] = true;
        let mut last = start;
        while let Some(u) = q.pop_front() {
            last = u;
            for &v in &adj[u] {
                if !seen[v as usize] {
                    seen[v as usize] = true;
                    q.push_back(v as usize);
                }
            }
        }
        last
    };
    let hub = (0..mem.len()).max_by_key(|&i| adj[i].len())?;
    let a = far(hub, &adj);
    let b = far(a, &adj);
    if a == b {
        return None;
    }
    // 0 = a's side, 1 = b's side, 2 = undecided
    let mut lab = vec![2u8; mem.len()];
    lab[a] = 0;
    lab[b] = 1;
    for _ in 0..30 {
        let mut moved = 0usize;
        for i in 0..mem.len() {
            if i == a || i == b {
                continue;
            }
            let (mut c0, mut c1) = (0u32, 0u32);
            for &v in &adj[i] {
                match lab[v as usize] {
                    0 => c0 += 1,
                    1 => c1 += 1,
                    _ => {}
                }
            }
            let new = if c0 > c1 { 0 } else if c1 > c0 { 1 } else { lab[i] };
            if new != lab[i] {
                lab[i] = new;
                moved += 1;
            }
        }
        if moved == 0 {
            break;
        }
    }
    let side: Vec<bool> = lab.iter().map(|&l| l == 1).collect();
    let n1 = side.iter().filter(|s| **s).count();
    let n0 = side.len() - n1;
    if n0.min(n1) < (min_side_frac * side.len() as f32) as usize {
        return None;
    }
    let (mut w0, mut w1, mut cross) = (0usize, 0usize, 0usize);
    for e in edges {
        if let (Some(&ia), Some(&ib)) = (idx.get(&e.a), idx.get(&e.b)) {
            match (side[ia], side[ib]) {
                (false, false) => w0 += 1,
                (true, true) => w1 += 1,
                _ => cross += 1,
            }
        }
    }
    let inner = w0.min(w1);
    if inner == 0 || (cross as f32) > max_cut * inner as f32 {
        return None;
    }
    Some(side)
}
