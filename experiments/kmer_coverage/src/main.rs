//! Exact whole-library 31-mer coverage for the reads in a few blocks of a paired run.
//!
//! WHY. The residue-structure experiment asks whether residue reads behave like pieces of
//! a genome (sampled from a molecule present at some depth, so their k-mers recur in the
//! rest of the library at a consistent level along the read and between mates) or like
//! error (k-mers unique to the read). A count-min sketch only ever overstates counts, which
//! is safe for "rare" but blurs exactly the flatness and mate-agreement readouts this
//! needs, so counts here are exact for every k-mer that occurs in a query block, over the
//! whole library.
//!
//! Each pair's own occurrences are subtracted, so a count is "occurrences elsewhere in the
//! library". This matters because 53% of this library's pairs overlap their mate: without
//! it every k-mer in the overlap would look like it had a second witness.
//!
//! A planted positive control can be passed as a TSV of `name\tr1\tr2`: its pairs are
//! added to the library (so they witness each other exactly as a real low-abundance genome's
//! reads would) and every one is scored, as block 0.
//!
//! Usage: kmer_coverage R1.gz R2.gz TRIM_R2 PAIRS_PER_BLOCK SKIP1,SKIP2,... [SPIKE.tsv|-] [QUERY_ONLY.tsv] > out.tsv
//!
//! QUERY_ONLY pairs (block 9) are scored against the library but never added to it: that is
//! how the shuffled null is measured, so a null read cannot witness its own siblings.
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader, Write};

const K: usize = 31;
const MASK: u64 = (1u64 << (2 * K)) - 1;

fn enc(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Canonical 31-mers in read order; None where the window holds an N.
fn kmers(s: &[u8], out: &mut Vec<Option<u64>>) {
    out.clear();
    if s.len() < K {
        return;
    }
    let (mut f, mut r, mut valid) = (0u64, 0u64, 0usize);
    for (i, &b) in s.iter().enumerate() {
        match enc(b) {
            Some(c) => {
                f = ((f << 2) | c) & MASK;
                r = (r >> 2) | ((3 - c) << (2 * (K - 1)));
                valid += 1;
            }
            None => valid = 0,
        }
        if i + 1 >= K {
            out.push(if valid >= K { Some(f.min(r)) } else { None });
        }
    }
}

struct Pairs {
    h1: Box<dyn BufRead>,
    h2: Box<dyn BufRead>,
    trim: usize,
}
impl Pairs {
    /// Opens either gzip or plain FASTQ, decided by the file's own magic bytes.
    ///
    /// WHY SNIFF RATHER THAN TAKE A FLAG. The cross-sample recurrence pass scans one subject
    /// library per panel member, and every one of them arrives as plain FASTQ from
    /// `fasterq-dump`. Re-compressing 130 GB only so this function can decompress it again
    /// is a gzip pass per library for nothing. A flag would do, but the magic bytes cannot
    /// be passed wrong: a plain file opened as gzip fails at the first read, and a gzip file
    /// read as text yields no parsable record, so getting this from the file removes a way
    /// for a 116-library run to die an hour in.
    fn open(a: &str, b: &str, trim: usize) -> Self {
        let o = |p: &str| -> Box<dyn BufRead> {
            let mut probe = std::fs::File::open(p).unwrap();
            let mut magic = [0u8; 2];
            use std::io::Read;
            let gz = probe.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b];
            let file = std::fs::File::open(p).unwrap();
            if gz {
                Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
            } else {
                Box::new(BufReader::with_capacity(1 << 20, file))
            }
        };
        Pairs { h1: o(a), h2: o(b), trim }
    }
    fn next(&mut self) -> Option<(String, Vec<u8>, Vec<u8>)> {
        let mut l = [String::new(), String::new(), String::new(), String::new()];
        let mut m = [String::new(), String::new(), String::new(), String::new()];
        for i in 0..4 {
            if self.h1.read_line(&mut l[i]).ok()? == 0 {
                return None;
            }
            self.h2.read_line(&mut m[i]).ok()?;
        }
        let name = l[0][1..].split_whitespace().next()?.to_string();
        let s1 = l[1].trim_end().as_bytes().to_vec();
        let s2 = m[1].trim_end().as_bytes();
        let s2 = if s2.len() > self.trim { s2[self.trim..].to_vec() } else { Vec::new() };
        Some((name, s1, s2))
    }
}

fn q(v: &[u32], p: f64) -> u32 {
    if v.is_empty() {
        return 0;
    }
    let mut s = v.to_vec();
    s.sort_unstable();
    s[((s.len() - 1) as f64 * p).round() as usize]
}

/// Largest |log2| ratio of left-median to right-median over breakpoints with >=20 k-mers a
/// side: a step in coverage along one read is what a chimaera or a junction looks like.
fn step(v: &[u32]) -> f64 {
    let mut best = 0f64;
    if v.len() < 40 {
        return 0.0;
    }
    let mut i = 20;
    while i + 20 <= v.len() {
        let a = q(&v[..i], 0.5) as f64;
        let b = q(&v[i..], 0.5) as f64;
        let d = ((a + 1.0) / (b + 1.0)).log2().abs();
        if d > best {
            best = d;
        }
        i += 5;
    }
    best
}

/// Longest run of consecutive k-mers with no witness anywhere else.
///
/// One substitution error destroys the (up to) 31 k-mers spanning it, so a zero-run near 31
/// in an otherwise covered read is the signature of an error rather than of a rare molecule,
/// which is unwitnessed end to end. This is what separates the two explanations of a
/// coverage step.
fn zero_run(v: &[u32]) -> usize {
    let (mut best, mut cur) = (0usize, 0usize);
    for &c in v {
        cur = if c == 0 { cur + 1 } else { 0 };
        if cur > best {
            best = cur;
        }
    }
    best
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (r1, r2) = (&a[1], &a[2]);
    let trim: usize = a[3].parse().unwrap();
    let per: usize = a[4].parse().unwrap();
    let mut skips: Vec<usize> = a[5].split(',').map(|x| x.parse().unwrap()).collect();
    skips.sort();

    // Pass 1: collect query pairs and their k-mers.
    let mut q_pairs: Vec<(usize, String, Vec<Option<u64>>, Vec<Option<u64>>)> = Vec::new();
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    let mut p = Pairs::open(r1, r2, trim);
    let mut idx = 0usize;
    let last = *skips.last().unwrap() + per;
    let mut buf = Vec::new();
    while let Some((name, s1, s2)) = p.next() {
        if idx >= last {
            break;
        }
        if let Some(bi) = skips.iter().position(|&s| idx >= s && idx < s + per) {
            kmers(&s1, &mut buf);
            let k1 = buf.clone();
            kmers(&s2, &mut buf);
            let k2 = buf.clone();
            for k in k1.iter().chain(k2.iter()).flatten() {
                counts.entry(*k).or_insert(0);
            }
            q_pairs.push((bi + 1, name, k1, k2));
        }
        idx += 1;
    }
    let mut spike: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    if a.len() > 6 && a[6] != "-" {
        for line in std::fs::read_to_string(&a[6]).unwrap().lines() {
            let f: Vec<&str> = line.split('\t').collect();
            let (s1, s2) = (f[1].as_bytes().to_vec(), f[2].as_bytes().to_vec());
            kmers(&s1, &mut buf);
            let k1 = buf.clone();
            kmers(&s2, &mut buf);
            let k2 = buf.clone();
            for k in k1.iter().chain(k2.iter()).flatten() {
                counts.entry(*k).or_insert(0);
            }
            q_pairs.push((0, f[0].to_string(), k1, k2));
            spike.push((s1, s2));
        }
    }
    if a.len() > 7 {
        for line in std::fs::read_to_string(&a[7]).unwrap().lines() {
            let f: Vec<&str> = line.split('\t').collect();
            kmers(f[1].as_bytes(), &mut buf);
            let k1 = buf.clone();
            kmers(f[2].as_bytes(), &mut buf);
            let k2 = buf.clone();
            for k in k1.iter().chain(k2.iter()).flatten() {
                counts.entry(*k).or_insert(0);
            }
            q_pairs.push((9, f[0].to_string(), k1, k2));
        }
    }
    eprintln!("query pairs {} distinct kmers {} spike pairs {}", q_pairs.len(), counts.len(), spike.len());

    // Pass 2: whole library.
    let mut p = Pairs::open(r1, r2, trim);
    let mut n = 0usize;
    while let Some((_, s1, s2)) = p.next() {
        for s in [&s1, &s2] {
            kmers(s, &mut buf);
            for k in buf.iter().flatten() {
                if let Some(c) = counts.get_mut(k) {
                    *c = c.saturating_add(1);
                }
            }
        }
        n += 1;
    }
    for (s1, s2) in &spike {
        for s in [s1, s2] {
            kmers(s, &mut buf);
            for k in buf.iter().flatten() {
                if let Some(c) = counts.get_mut(k) {
                    *c = c.saturating_add(1);
                }
            }
        }
    }
    eprintln!("library pairs {} (+{} planted)", n, spike.len());

    let out = std::io::stdout();
    let mut o = std::io::BufWriter::new(out.lock());
    writeln!(o, "block\tname\tn1\tn2\tmed1\tmed2\tq10_1\tq90_1\tq10_2\tq90_2\tzero1\tzero2\tstep1\tstep2\tzrun1\tzrun2\tstepnz1\tstepnz2\tlibrary_pairs").unwrap();
    for (b, name, k1, k2) in &q_pairs {
        let mut own: FxHashMap<u64, u32> = FxHashMap::default();
        for k in k1.iter().chain(k2.iter()).flatten() {
            *own.entry(*k).or_insert(0) += 1;
        }
        let ext = |ks: &Vec<Option<u64>>| -> Vec<u32> {
            ks.iter().flatten().map(|k| if *b == 9 { counts[k] } else { counts[k].saturating_sub(own[k]) }).collect()
        };
        let (e1, e2) = (ext(k1), ext(k2));
        let z = |e: &Vec<u32>| if e.is_empty() { f64::NAN } else { e.iter().filter(|&&c| c == 0).count() as f64 / e.len() as f64 };
        // step on the non-zero k-mers only: if the step survives dropping them it is a
        // real change of depth, and if it collapses it was unwitnessed sequence.
        let nz1: Vec<u32> = e1.iter().copied().filter(|&c| c > 0).collect();
        let nz2: Vec<u32> = e2.iter().copied().filter(|&c| c > 0).collect();
        writeln!(
            o,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.3}\t{:.3}\t{}\t{}\t{:.3}\t{:.3}\t{}",
            b, name, e1.len(), e2.len(), q(&e1, 0.5), q(&e2, 0.5), q(&e1, 0.1), q(&e1, 0.9),
            q(&e2, 0.1), q(&e2, 0.9), z(&e1), z(&e2), step(&e1), step(&e2),
            zero_run(&e1), zero_run(&e2), step(&nz1), step(&nz2), n
        )
        .unwrap();
    }
}
