//! End-to-end: two synonymous variants of one coding genome, reads from both
//! strands, must assemble into one contig carrying the variant sites.

use hapax::consensus::{self, ConsensusParams};
use hapax::genetic::{aa, base_of, codon_index, is_stop};
use hapax::index::{Index, SeedParams};
use hapax::layout;
use hapax::model::{ModelParams, Scorer};
use hapax::overlap::{self, DispatchStats, OverlapParams};
use hapax::seq::revcomp;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

#[test]
fn swarm_of_two_becomes_one_contig() {
    let mut r = Lcg(42);
    let mut g = Vec::new();
    while g.len() < 3000 {
        let c = (r.next() % 64) as usize;
        if !is_stop(c) {
            (0..3).for_each(|p| g.push(base_of(c, p)));
        }
    }
    // variant: synonymous change at every 7th codon where one exists
    let mut v = g.clone();
    let mut sites = Vec::new();
    for c in (0..g.len() / 3).step_by(7) {
        let i = c * 3 + 2;
        let orig = codon_index(&g[c * 3..c * 3 + 3]).unwrap();
        for nb in 0..4u8 {
            let t = [g[c * 3], g[c * 3 + 1], nb];
            if nb != g[i] && aa(codon_index(&t).unwrap()) == aa(orig) {
                v[i] = nb;
                sites.push(i);
                break;
            }
        }
    }
    let mut reads = Vec::new();
    for (k, s) in (0..g.len() - 150).step_by(7).enumerate() {
        let src = if k % 2 == 0 { &g } else { &v };
        let rd = src[s..s + 150].to_vec();
        reads.push(if k % 3 == 0 { revcomp(&rd) } else { rd });
    }
    let idx = Index::build(&reads, SeedParams::default());
    let op = OverlapParams::default();
    let cands = overlap::candidates(&reads, &idx, &op);
    let st = DispatchStats::default();
    let edges = overlap::dispatch(&reads, &cands, &Scorer::new(ModelParams::default()), &op, &st).edges;
    let lens: Vec<usize> = reads.iter().map(|r| r.len()).collect();
    let groups = layout::group(reads.len(), &edges, 6, 5, 0.1, None);
    let comps = layout::layout(&lens, &edges, Some(&groups));
    let sizes: Vec<usize> = comps.iter().map(|c| c.len()).collect();
    eprintln!("sizes {:?} edges {} t1 {} t2 {} t2acc {}", &sizes[..sizes.len().min(10)], edges.len(), st.get(&st.t1_accept), st.get(&st.t2), st.get(&st.t2_accept));
    assert_eq!(comps.iter().filter(|c| c.len() > 1).count(), 1, "expected one component");
    let c = consensus::build(&reads, &comps[0], &ConsensusParams::default());
    assert!(c.seq.len() > 2900 && c.seq.len() < 3100, "len {}", c.seq.len());
    // most synonymous sites are reported as variants, annotated synonymous
    let called = c.variants.iter().filter(|x| x.effects.iter().any(|e| e.1)).count();
    assert!(called as f64 > 0.7 * sites.len() as f64, "called {} of {}", called, sites.len());
    assert!(st.hmm_fraction() > 0.0);
}
