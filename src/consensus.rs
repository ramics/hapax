//! Variant-aware consensus.
//!
//! A single consensus throws away exactly what this tool exists to keep, so a
//! contig here is a consensus *plus* its variant structure: per-column allele
//! counts, the variant sites they imply, the linkage between neighbouring
//! sites as observed within single reads, and an annotation of which regions
//! look coding (one frame), dual-coding (two frames) or untranslated. That is
//! the output the later closure-scoring work needs: it can ask whether a
//! contig is one entity or several without going back to the reads.
//!
//! Polishing is plain nucleotide alignment of each read to the draft. By the
//! time reads reach here the codon model has already done its job (deciding
//! which reads belong together); placement against the draft is a much easier
//! problem and the nucleotides are what get emitted.

use crate::genetic::{aa, codon_index};
use crate::layout::Placement;
pub use crate::orf::{Orf, OrfCaller};
use crate::seq::revcomp;
use rayon::prelude::*;

const MATCH: i32 = 2;
const MISMATCH: i32 = -3;
const GAP_OPEN: i32 = -5;
const GAP_EXT: i32 = -2;
const PAD: i64 = 40;
/// Cost of leaving a read end unaligned. See `align_to_ref`.
const CLIP: i32 = 10;

#[derive(Clone, Debug)]
pub struct ReadAln {
    pub read: u32,
    pub rev: bool,
    pub ref_start: usize,
    pub ref_end: usize,
    /// read base aligned to each ref column in [ref_start, ref_end); 4 = deletion
    pub cols: Vec<u8>,
    /// insertions: (ref column after which it sits, inserted bases)
    pub ins: Vec<(usize, Vec<u8>)>,
    pub score: i32,
    /// read bases left unaligned before ref_start / after ref_end (see `CLIP`)
    pub clip_l: Vec<u8>,
    pub clip_r: Vec<u8>,
}

/// Affine alignment of a read into a reference window. The window has free ends;
/// each read end may be left unaligned (clipped) for a fixed cost `CLIP`.
///
/// Why clipping. Aligning every read base end to end forces anything at a read end
/// that is not template into the pileup as mismatches and indels. On real libraries
/// that is common: in Berlin wastewater run SRR23973193, 96% of mate-2 reads start
/// with ~12 bp of non-template sequence (random-primer/linker; minimap2 soft-clips
/// it), at full base quality, so quality trimming does not remove it. End to end,
/// those bases made mate-2 reads look ~9% divergent from the draft and left 9-241 bp
/// junk indels in the consensus, the main cause of hapax's 97.0% identity there.
/// CLIP = 10 matches minimap2's short-read end bonus: a tail is clipped only if
/// aligning it scores below -10, which a genuine SNP near a read end (-3) never does.
/// Clipped bases are not voted; the read itself is kept whole.
pub fn align_to_ref(q: &[u8], r: &[u8]) -> Option<(usize, usize, Vec<u8>, Vec<(usize, Vec<u8>)>, i32)> {
    align_clipped(q, r).map(|(s, e, c, i, sc, _, _)| (s, e, c, i, sc))
}

/// `align_to_ref`, also returning how many read bases were clipped at each end.
#[allow(clippy::type_complexity)]
fn align_clipped(q: &[u8], r: &[u8]) -> Option<(usize, usize, Vec<u8>, Vec<(usize, Vec<u8>)>, i32, usize, usize)> {
    let m = q.len();
    let n = r.len();
    if m == 0 || n == 0 {
        return None;
    }
    const NEG: i32 = i32::MIN / 4;
    let w = n + 1;
    let mut h = vec![NEG; (m + 1) * w];
    let mut e = vec![NEG; (m + 1) * w]; // gap in read (consume ref) -> deletion
    let mut f = vec![NEG; (m + 1) * w]; // gap in ref (consume read) -> insertion
    let mut tb = vec![0u8; (m + 1) * w];
    for j in 0..=n {
        h[j] = 0;
    }
    for i in 1..=m {
        f[i * w] = GAP_OPEN + GAP_EXT * (i as i32 - 1);
        h[i * w] = f[i * w];
        tb[i * w] = 2 | (1 << 3);
        if -CLIP > h[i * w] {
            h[i * w] = -CLIP;
            tb[i * w] = 3 | (1 << 3);
        }
        for j in 1..=n {
            let c = i * w + j;
            // E: deletion (move in j)
            let e_open = h[c - 1] + GAP_OPEN;
            let e_ext = e[c - 1] + GAP_EXT;
            let (ev, eb) = if e_ext > e_open { (e_ext, 1u8) } else { (e_open, 0u8) };
            e[c] = ev;
            // F: insertion (move in i)
            let f_open = h[c - w] + GAP_OPEN;
            let f_ext = f[c - w] + GAP_EXT;
            let (fv, fb) = if f_ext > f_open { (f_ext, 1u8) } else { (f_open, 0u8) };
            f[c] = fv;
            let s = if q[i - 1] == r[j - 1] && q[i - 1] < 4 { MATCH } else { MISMATCH };
            let dv = h[c - w - 1] + s;
            let (hv, hb) = if dv >= ev && dv >= fv { (dv, 0u8) } else if ev >= fv { (ev, 1u8) } else { (fv, 2u8) };
            // hb = 3: the alignment starts after this cell, read[..i] clipped
            let (hv, hb) = if -CLIP > hv { (-CLIP, 3u8) } else { (hv, hb) };
            h[c] = hv;
            tb[c] = hb | (eb << 2) | (fb << 3);
        }
    }
    // best end: the full read (row m), or a clipped tail at cost CLIP
    let mut best = NEG;
    let (mut i, mut j) = (m, 0);
    for ii in 1..=m {
        let pen = if ii < m { CLIP } else { 0 };
        for jj in 0..=n {
            let v = h[ii * w + jj] - pen;
            if v > best || (v == best && ii == m && i != m) {
                best = v;
                i = ii;
                j = jj;
            }
        }
    }
    let end = j;
    let rclip = m - i;
    let mut state = 0u8; // 0=H, 1=E, 2=F
    let mut ops: Vec<(u8, u8)> = Vec::new(); // (op, base) op 0=match/mismatch, 1=del, 2=ins
    while i > 0 {
        let c = i * w + j;
        match state {
            0 => {
                let hb = tb[c] & 3;
                if hb == 3 {
                    break; // read[..i] is clipped
                }
                if hb == 0 {
                    ops.push((0, q[i - 1]));
                    i -= 1;
                    j -= 1;
                } else {
                    state = hb;
                }
            }
            1 => {
                ops.push((1, 4));
                let ext = (tb[c] >> 2) & 1 == 1;
                j -= 1;
                if !ext {
                    state = 0;
                }
            }
            _ => {
                ops.push((2, q[i - 1]));
                let ext = (tb[c] >> 3) & 1 == 1;
                i -= 1;
                if !ext {
                    state = 0;
                }
            }
        }
        if j == 0 && i > 0 && state == 0 && tb[i * w] & 3 == 2 {
            // remaining read bases are an insertion before the window
            state = 2;
        }
    }
    ops.reverse();
    let lclip = i;
    let start = j;
    let mut cols = Vec::with_capacity(end - start);
    let mut ins: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut col = start;
    for (op, b) in ops {
        match op {
            0 | 1 => {
                cols.push(b);
                col += 1;
            }
            _ => {
                let after = col.saturating_sub(1);
                if let Some(last) = ins.last_mut() {
                    if last.0 == after && col > start {
                        last.1.push(b);
                        continue;
                    }
                }
                ins.push((after, vec![b]));
            }
        }
    }
    Some((start, end, cols, ins, best, lclip, rclip))
}

#[derive(Clone, Debug, Default)]
pub struct Column {
    pub counts: [u32; 5],
}

impl Column {
    pub fn depth(&self) -> u32 {
        self.counts.iter().sum()
    }
}

#[derive(Clone, Debug)]
pub struct Variant {
    pub pos: usize,
    pub depth: u32,
    pub counts: [u32; 5],
    pub cons: u8,
    pub alt: u8,
    pub alt_freq: f32,
    pub region: Region,
    /// per coding frame covering the site: (frame, synonymous?)
    pub effects: Vec<(u8, bool)>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Region {
    Noncoding,
    Coding,
    Dual,
}

impl Region {
    pub fn name(&self) -> &'static str {
        match self {
            Region::Noncoding => "noncoding",
            Region::Coding => "coding",
            Region::Dual => "dual",
        }
    }
}


#[derive(Clone, Debug)]
pub struct Contig {
    pub seq: Vec<u8>,
    pub columns: Vec<Column>,
    pub alns: Vec<ReadAln>,
    pub unplaced: Vec<u32>,
    pub variants: Vec<Variant>,
    pub orfs: Vec<Orf>,
    pub region: Vec<Region>,
}

#[derive(Clone, Debug)]
pub struct ConsensusParams {
    pub rounds: usize,
    pub min_alt: u32,
    pub min_af: f32,
    pub min_depth: u32,
    pub min_orf: usize,
    pub orf_caller: OrfCaller,
}

impl Default for ConsensusParams {
    fn default() -> Self {
        ConsensusParams { rounds: 3, min_alt: 2, min_af: 0.05, min_depth: 4, min_orf: 300, orf_caller: OrfCaller::Scored }
    }
}

fn argmax(c: &[u32; 5]) -> u8 {
    let mut b = 0;
    for i in 1..5 {
        if c[i] > c[b] {
            b = i;
        }
    }
    b as u8
}

fn realign(reads: &[Vec<u8>], pl: &[Placement], draft: &[u8]) -> (Vec<ReadAln>, Vec<u32>) {
    let res: Vec<Result<ReadAln, u32>> = pl
        .par_iter()
        .map(|p| {
            let r = &reads[p.read as usize];
            let q = if p.rev { revcomp(r) } else { r.clone() };
            let ws = (p.x - PAD).max(0) as usize;
            let we = ((p.x + q.len() as i64 + PAD) as usize).min(draft.len());
            if ws >= we {
                return Err(p.read);
            }
            match align_clipped(&q, &draft[ws..we]) {
                Some((s, e, cols, ins, score, lc, rc)) if score > (q.len() as i32) / 2 => Ok(ReadAln {
                    read: p.read,
                    rev: p.rev,
                    ref_start: ws + s,
                    ref_end: ws + e,
                    cols,
                    ins: ins.into_iter().map(|(c, v)| (c + ws, v)).collect(),
                    score,
                    clip_l: q[..lc].to_vec(),
                    clip_r: q[q.len() - rc..].to_vec(),
                }),
                _ => Err(p.read),
            }
        })
        .collect();
    let mut alns = Vec::new();
    let mut un = Vec::new();
    for r in res {
        match r {
            Ok(a) => alns.push(a),
            Err(x) => un.push(x),
        }
    }
    (alns, un)
}

fn pileup(len: usize, alns: &[ReadAln]) -> (Vec<Column>, Vec<std::collections::HashMap<Vec<u8>, u32>>) {
    let mut cols = vec![Column::default(); len];
    let mut ins: Vec<std::collections::HashMap<Vec<u8>, u32>> = vec![Default::default(); len];
    for a in alns {
        for (k, &b) in a.cols.iter().enumerate() {
            cols[a.ref_start + k].counts[b.min(4) as usize] += 1;
        }
        for (c, s) in &a.ins {
            if *c < len {
                *ins[*c].entry(s.clone()).or_insert(0) += 1;
            }
        }
    }
    (cols, ins)
}

/// Build a new consensus by per-column majority; returns it with a map old col -> new col.
fn call(len: usize, alns: &[ReadAln], draft: &[u8]) -> (Vec<u8>, Vec<usize>) {
    let (cols, ins) = pileup(len, alns);
    let fill = missing_sequence(alns, draft, &cols);
    let mut out = Vec::with_capacity(len);
    let mut map = Vec::with_capacity(len);
    let mut skip_to = 0usize;
    for (c, col) in cols.iter().enumerate() {
        map.push(out.len());
        let d = col.depth();
        if d == 0 || c < skip_to {
            continue;
        }
        let b = argmax(&col.counts);
        if b < 4 {
            out.push(b);
        }
        if let Some((s, &n)) = ins[c].iter().max_by_key(|x| *x.1) {
            let total_ins: u32 = ins[c].values().sum();
            if 2 * total_ins > d && 2 * n > total_ins {
                out.extend_from_slice(s);
            }
        }
        if let Some((drop, f)) = fill.get(&c) {
            out.extend_from_slice(f);
            skip_to = c + 1 + drop;
        }
    }
    (out, map)
}

/// Sequence missing from the draft, recovered from clipped read ends.
///
/// Why. Read-end clipping (`CLIP`) keeps non-template ends out of the vote, but
/// it also means a read spanning a stretch the draft lacks or has wrong clips
/// there instead of carrying a long insertion -- clipping costs 10, a 45 bp gap
/// 93. Long insertions were never callable anyway: reads cover them partially,
/// so no single inserted string wins a majority. On Berlin that left a 45 bp
/// hole in the target-lineage contig with 68 reads spanning it
/// (results/ww_berlin_consensus).
///
/// Rule: a column where at least `FILL_MIN` reads, and `FILL_FRAC` of the
/// depth, are clipped on the right is a left edge of a hole; the same on the
/// left is a right edge. The clipped tails' per-offset majority `t` reads into
/// the hole from the left, the heads' majority `h` from the right. A repair is
/// made only when it is anchored on both sides by agreement with sequence:
/// `t` rejoining the draft after its edge (an insertion), or `t` and `h`
/// overlapping each other by `FILL_ANCHOR` bases (the hole, possibly with wrong
/// draft columns between the two edges, is replaced by `t` + `h`). Clipped
/// non-template ends never pass: they differ between reads, so no per-offset
/// majority forms. Returns, per left-edge column, (columns to drop after it,
/// bases to insert).
fn missing_sequence(alns: &[ReadAln], draft: &[u8], cols: &[Column]) -> std::collections::HashMap<usize, (usize, Vec<u8>)> {
    use std::collections::HashMap;
    const FILL_MIN: usize = 4;
    const FILL_FRAC: f32 = 0.25;
    const FILL_ANCHOR: usize = 12;
    const FILL_SPAN: usize = 200;
    let mut tails: HashMap<usize, Vec<Vec<u8>>> = HashMap::new();
    let mut heads: HashMap<usize, Vec<Vec<u8>>> = HashMap::new();
    for a in alns {
        if a.ref_end <= a.ref_start {
            continue;
        }
        if a.clip_r.len() >= FILL_ANCHOR && a.ref_end < draft.len() {
            tails.entry(a.ref_end - 1).or_default().push(a.clip_r.clone());
        }
        if a.clip_l.len() >= FILL_ANCHOR && a.ref_start > 0 {
            // reversed, so offset 0 is the base next to the edge
            heads.entry(a.ref_start).or_default().push(a.clip_l.iter().rev().copied().collect());
        }
    }
    let enough = |n: usize, c: usize| n >= FILL_MIN && n as f32 >= FILL_FRAC * cols[c].depth() as f32;
    let consensus = |ts: &[Vec<u8>]| {
        let mut out = Vec::new();
        for k in 0.. {
            let mut cnt = [0usize; 4];
            for t in ts {
                if k < t.len() && t[k] < 4 {
                    cnt[t[k] as usize] += 1;
                }
            }
            let n: usize = cnt.iter().sum();
            let (b, top) = cnt.iter().enumerate().max_by_key(|x| *x.1).map(|(b, &v)| (b, v)).unwrap();
            if n < FILL_MIN || 2 * top <= n {
                break;
            }
            out.push(b as u8);
        }
        out
    };
    let close = |x: &[u8], y: &[u8]| x.len() == y.len() && x.iter().zip(y).filter(|(a, b)| a != b).count() <= 1.max(x.len() / 20);
    let head_cons: HashMap<usize, Vec<u8>> = heads
        .iter()
        .filter(|(c, h)| enough(h.len(), **c))
        .map(|(c, h)| (*c, consensus(h).into_iter().rev().collect::<Vec<u8>>()))
        .collect();
    let mut out = HashMap::new();
    let mut lefts: Vec<usize> = tails.keys().copied().filter(|c| enough(tails[c].len(), *c)).collect();
    lefts.sort_unstable();
    let mut done_to = 0usize;
    for c in lefts {
        if c < done_to {
            continue;
        }
        let t = consensus(&tails[&c]);
        if t.len() < FILL_ANCHOR {
            continue;
        }
        // (a) t rejoins the draft: insertion of t[..l]; l = 0 means nothing is missing
        let mut fixed = false;
        for l in 0..=t.len() - FILL_ANCHOR {
            if c + 1 + FILL_ANCHOR > draft.len() {
                break;
            }
            if close(&t[l..l + FILL_ANCHOR], &draft[c + 1..c + 1 + FILL_ANCHOR]) {
                if l > 0 {
                    out.insert(c, (0, t[..l].to_vec()));
                }
                fixed = true;
                break;
            }
        }
        if fixed {
            continue;
        }
        // (b) bridge to the right edge with the most support within FILL_SPAN
        let mut edges: Vec<(&usize, &Vec<u8>)> = head_cons.iter().filter(|(e, _)| **e > c && **e <= c + FILL_SPAN).collect();
        edges.sort_by_key(|(e, _)| std::cmp::Reverse(heads[*e].len()));
        for (&e, h) in edges {
            // h ends at the edge; find the shift s with t[s..] overlapping h's start
            let mut bridge = None;
            for sft in 0..t.len() {
                let ov = (t.len() - sft).min(h.len());
                if ov < FILL_ANCHOR {
                    break;
                }
                if close(&t[sft..sft + ov], &h[..ov]) {
                    let mut seq = t[..sft].to_vec();
                    seq.extend_from_slice(h);
                    bridge = Some(seq);
                    break;
                }
            }
            // or h reaches back past the left edge: t's start lies inside h
            if bridge.is_none() {
                for sft in 1..h.len() {
                    let ov = (h.len() - sft).min(t.len());
                    if ov < FILL_ANCHOR {
                        break;
                    }
                    if close(&h[sft..sft + ov], &t[..ov]) {
                        bridge = Some(h[sft..].to_vec());
                        break;
                    }
                }
            }
            if let Some(seq) = bridge {
                out.insert(c, (e - c - 1, seq));
                done_to = e;
                break;
            }
        }
    }
    out
}

pub fn annotate_regions(len: usize, orfs: &[Orf]) -> Vec<Region> {
    let mut n = vec![0u8; len];
    for o in orfs {
        // clamped: a split read set (src/phase.rs) can leave a consensus
        // shorter than the sequence its ORFs were called on
        let (s, e) = (o.start.min(len), o.end.min(len));
        for x in &mut n[s..e] {
            *x += 1;
        }
    }
    n.into_iter()
        .map(|k| match k {
            0 => Region::Noncoding,
            1 => Region::Coding,
            _ => Region::Dual,
        })
        .collect()
}

fn effect(seq: &[u8], pos: usize, alt: u8, o: &Orf) -> Option<bool> {
    if o.frame < 3 {
        let cs = o.start + (pos - o.start) / 3 * 3;
        if cs + 3 > seq.len() {
            return None;
        }
        let ref_c = codon_index(&seq[cs..cs + 3])?;
        let mut c = [seq[cs], seq[cs + 1], seq[cs + 2]];
        c[pos - cs] = alt;
        Some(aa(codon_index(&c)?) == aa(ref_c))
    } else {
        // reverse frame: work on the reverse complement coordinates
        let rc = revcomp(&seq[o.start..o.end]);
        let rpos = o.end - 1 - pos;
        let cs = rpos / 3 * 3;
        if cs + 3 > rc.len() {
            return None;
        }
        let ref_c = codon_index(&rc[cs..cs + 3])?;
        let mut c = [rc[cs], rc[cs + 1], rc[cs + 2]];
        c[rpos - cs] = 3 - alt;
        Some(aa(codon_index(&c)?) == aa(ref_c))
    }
}

pub fn build(reads: &[Vec<u8>], comp: &[Placement], p: &ConsensusParams) -> Contig {
    // round 0 draft: ungapped column votes from the embedding
    let len = comp.iter().map(|pl| pl.x as usize + reads[pl.read as usize].len()).max().unwrap_or(0);
    let mut cols = vec![Column::default(); len];
    for pl in comp {
        let r = &reads[pl.read as usize];
        let q = if pl.rev { revcomp(r) } else { r.clone() };
        for (k, &b) in q.iter().enumerate() {
            cols[pl.x as usize + k].counts[b.min(4) as usize] += 1;
        }
    }
    let mut draft: Vec<u8> = cols.iter().filter(|c| c.depth() > 0).map(|c| argmax(&c.counts).min(3)).collect();
    let mut pl: Vec<Placement> = comp.to_vec();
    let mut alns = Vec::new();
    let mut unplaced = Vec::new();
    for round in 0..=p.rounds {
        let (a, un) = realign(reads, &pl, &draft);
        alns = a;
        unplaced = un;
        if round == p.rounds {
            break;
        }
        let (new, map) = call(draft.len(), &alns, &draft);
        // move placements into the new coordinates. Reads that failed to realign
        // this round keep their previous placement (mapped across) and are retried
        // against the next draft, rather than being dropped for good after one bad
        // round: early drafts are the noisiest, so that is when misses happen.
        let mv = |x: i64| map.get(x.max(0) as usize).copied().unwrap_or(map.last().copied().unwrap_or(0)) as i64;
        let old: std::collections::HashMap<u32, Placement> = pl.iter().map(|q| (q.read, *q)).collect();
        pl = alns.iter().map(|a| Placement { read: a.read, rev: a.rev, x: mv(a.ref_start as i64) }).collect();
        pl.extend(unplaced.iter().filter_map(|r| old.get(r)).map(|q| Placement { x: mv(q.x), ..*q }));
        draft = new;
    }
    let (columns, _) = pileup(draft.len(), &alns);
    let orfs = crate::orf::call_orfs(&draft, p.min_orf, p.orf_caller);
    let region = annotate_regions(draft.len(), &orfs);
    let mut variants = Vec::new();
    for (pos, col) in columns.iter().enumerate() {
        let d = col.depth();
        if d < p.min_depth {
            continue;
        }
        let cons = draft[pos];
        let mut alt = 5u8;
        let mut best = 0;
        for b in 0..5u8 {
            if b != cons && col.counts[b as usize] > best {
                best = col.counts[b as usize];
                alt = b;
            }
        }
        if alt < 5 && best >= p.min_alt && best as f32 / d as f32 >= p.min_af {
            let effects = if alt < 4 {
                orfs.iter().filter(|o| pos >= o.start && pos < o.end).filter_map(|o| effect(&draft, pos, alt, o).map(|s| (o.frame, s))).collect()
            } else {
                Vec::new()
            };
            variants.push(Variant { pos, depth: d, counts: col.counts, cons, alt, alt_freq: best as f32 / d as f32, region: region[pos], effects });
        }
    }
    Contig { seq: draft, columns, alns, unplaced, variants, orfs, region }
}

/// Within-read linkage between consecutive variant sites: haplotype counts and r^2.
pub fn linkage(c: &Contig) -> Vec<(usize, usize, Vec<((u8, u8), u32)>, f32)> {
    use std::collections::HashMap;
    let sites: Vec<usize> = c.variants.iter().map(|v| v.pos).collect();
    let mut tab: HashMap<(usize, usize), HashMap<(u8, u8), u32>> = HashMap::new();
    for a in &c.alns {
        let lo = sites.partition_point(|&s| s < a.ref_start);
        let hi = sites.partition_point(|&s| s < a.ref_end);
        for k in lo..hi.saturating_sub(1) {
            let (s1, s2) = (sites[k], sites[k + 1]);
            let b1 = a.cols[s1 - a.ref_start];
            let b2 = a.cols[s2 - a.ref_start];
            *tab.entry((s1, s2)).or_default().entry((b1, b2)).or_insert(0) += 1;
        }
    }
    let mut out: Vec<_> = tab
        .into_iter()
        .map(|((s1, s2), h)| {
            let mut v: Vec<((u8, u8), u32)> = h.into_iter().collect();
            v.sort_by(|x, y| y.1.cmp(&x.1));
            let r2 = r_squared(&v, c.seq[s1], c.seq[s2]);
            (s1, s2, v, r2)
        })
        .collect();
    out.sort_by_key(|x| x.0);
    out
}

/// r^2 between "consensus allele or not" at the two sites.
fn r_squared(h: &[((u8, u8), u32)], c1: u8, c2: u8) -> f32 {
    let (mut n11, mut n10, mut n01, mut n00) = (0f32, 0f32, 0f32, 0f32);
    for &((a, b), n) in h {
        let n = n as f32;
        match (a == c1, b == c2) {
            (true, true) => n11 += n,
            (true, false) => n10 += n,
            (false, true) => n01 += n,
            (false, false) => n00 += n,
        }
    }
    let t = n11 + n10 + n01 + n00;
    if t == 0.0 {
        return 0.0;
    }
    let (p1, q1) = ((n11 + n10) / t, (n11 + n01) / t);
    let d = n11 / t - p1 * q1;
    let den = p1 * (1.0 - p1) * q1 * (1.0 - q1);
    if den <= 0.0 {
        0.0
    } else {
        d * d / den
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode;
    #[test]
    fn semi_global_places_read() {
        let r = encode(b"TTTTTTTTTTACGTACGGATCCATGCAGGGTTTTTTTTTT");
        let q = encode(b"ACGTACGGATCAATGCAGGG");
        let (s, e, cols, ins, _) = align_to_ref(&q, &r).unwrap();
        assert_eq!(s, 10);
        assert_eq!(e, 30);
        assert_eq!(cols.len(), 20);
        assert!(ins.is_empty());
        let q2 = encode(b"ACGTACGGATCCCATGCAGGG");
        let (_, _, cols2, ins2, _) = align_to_ref(&q2, &r).unwrap();
        assert_eq!(cols2.len(), 20);
        assert_eq!(ins2.len(), 1);
    }

    #[test]
    fn non_template_read_ends_are_clipped_not_voted() {
        // 12 bp of non-template sequence before the true read, as on mate 2 of the
        // Berlin wastewater library: it must not appear as an insertion or mismatches.
        let r = encode(b"TTGCATTGCAACGTACGGATCCATGCAGGGTACCTTAGGCATTGACCGTA");
        let q = encode(b"GGGGCCCCGGGGACGTACGGATCCATGCAGGGTACCTTAGG");
        let (s, e, cols, ins, _) = align_to_ref(&q, &r).unwrap();
        assert_eq!((s, e), (10, 39));
        assert!(ins.is_empty());
        assert!(cols.iter().zip(&r[s..e]).all(|(a, b)| a == b));
        // a genuine SNP two bases from the end is kept, not clipped
        let q3 = encode(b"ACGTACGGATCCATGCAGGGTACCTTAGGCATTGACCGAA");
        let (s3, e3, cols3, _, _) = align_to_ref(&q3, &r).unwrap();
        assert_eq!((s3, e3), (10, 50));
        assert_eq!(cols3[38], 0);
    }
}
