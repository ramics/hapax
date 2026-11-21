//! The linkage census: for a query read, who else in the library covers the same
//! sequence, and what are THEIR mates?
//!
//! WHY THIS AND NOT COVERAGE. `kmer_coverage` answers how many times a query's k-mers
//! recur in the library. That is a depth, and a depth cannot say what a read sits next to.
//! This asks the adjacency question instead: find the library pairs that share sequence
//! with the query, then look at the OTHER mate of each of those pairs. That other mate is
//! a piece of sequence physically linked - same molecule, one insert away - to something
//! the query also covers. Label those mates and the residue splits into two populations
//! that matter: reads that sit beside sequence we can name (rRNA precursor spacers, UTRs,
//! intergenic sequence, the untranslated flanks of things we identified) and reads whose
//! neighbourhood is as unknown as they are, which are whole unknown molecules.
//!
//! WHAT IS EMITTED. One row per (query, witnessing library pair): which mate of the
//! witness matched, how many 31-mers it shared, and how many the witness's OTHER mate
//! shared. The labelling and every aggregate is left to Python, because the labels come
//! from a table built elsewhere and this binary must not need rebuilding when that table's
//! classes change.
//!
//! THE CONTROL IS THE CALLER'S JOB AND MUST NOT BE SKIPPED. Run the same census over
//! already-identified reads. A residue read's neighbourhood looking unknown means nothing
//! until it is set beside what an identified read's neighbourhood looks like, measured the
//! same way in the same pass - in this project a separator that reads the same on both
//! piles has turned out, every time, to be describing the database.
//!
//! BIASES, BOTH DIRECTIONS NAMED.
//!  - A witness found through a SHARED k-mer is biased towards abundant sequence: a
//!    singleton molecule has no witnesses at all and lands in "no witnesses", which is a
//!    different statement from "whole unknown molecule" and is reported separately.
//!  - About a third of this library's pairs overlap their own mate, so for those the
//!    "other mate" is partly the same sequence. Rows carry both mates' shared counts so a
//!    caller can drop or split on that rather than averaging over it.
//!
//! Usage: census R1.gz R2.gz TRIM_R2 QUERIES.tsv MIN_SHARED MAX_WITNESSES > rows.tsv
//!   QUERIES.tsv is name\tseq1\tseq2, the same format kmer_coverage's query-only file uses.
use flate2::read::MultiGzDecoder;
use rustc_hash::{FxHashMap, FxHashSet};
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

fn kmers(s: &[u8], out: &mut Vec<u64>) {
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
        if i + 1 >= K && valid >= K {
            out.push(f.min(r));
        }
    }
}

struct Pairs {
    h1: Box<dyn BufRead>,
    h2: Box<dyn BufRead>,
    trim: usize,
}
impl Pairs {
    fn open(a: &str, b: &str, trim: usize) -> Self {
        let o = |p: &str| -> Box<dyn BufRead> {
            Box::new(BufReader::with_capacity(
                1 << 20,
                MultiGzDecoder::new(std::fs::File::open(p).unwrap()),
            ))
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

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 7 {
        eprintln!("usage: census R1.gz R2.gz TRIM_R2 QUERIES.tsv MIN_SHARED MAX_WITNESSES");
        std::process::exit(2);
    }
    let (r1, r2) = (&a[1], &a[2]);
    let trim: usize = a[3].parse().unwrap();
    let min_shared: usize = a[5].parse().unwrap();
    let max_wit: usize = a[6].parse().unwrap();

    // Pass 1: index the queries' k-mers. One entry per distinct query k-mer, pointing at
    // the queries that carry it; a k-mer in many queries is kept whole rather than
    // deduplicated, because dropping it would silently remove exactly the conserved
    // sequence this census needs to see.
    let mut names: Vec<String> = Vec::new();
    let mut index: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
    let mut buf: Vec<u64> = Vec::new();
    for line in std::fs::read_to_string(&a[4]).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let qi = names.len() as u32;
        names.push(f[0].to_string());
        let mut seen: FxHashSet<u64> = FxHashSet::default();
        for s in [f[1].as_bytes(), f[2].as_bytes()] {
            kmers(s, &mut buf);
            for k in buf.iter() {
                if seen.insert(*k) {
                    index.entry(*k).or_default().push(qi);
                }
            }
        }
    }
    eprintln!("queries {} distinct query kmers {}", names.len(), index.len());

    // Pass 2: stream the library. For each pair, work out which queries each of its mates
    // touches, then emit a row per query it witnesses. A pair that IS the query (same
    // name) is skipped: a read witnessing itself is not adjacency.
    // Witnesses are RESERVOIR SAMPLED per query rather than taken in file order. The
    // library is in flow-cell order, so the first N witnesses of a query are one region of
    // one lane, and a cap applied in stream order would quietly make the census a
    // statement about that region. The reservoir is seeded, so the sample is reproducible
    // from the seed like everything else here.
    let mut seen_count: Vec<u64> = vec![0; names.len()];
    let mut keep: Vec<Vec<(String, u32, u32, Vec<u8>)>> = vec![Vec::new(); names.len()];
    let mut rng: u64 = a.get(7).and_then(|x| x.parse().ok()).unwrap_or(1);
    let mut next_rand = move || {
        // xorshift64*, enough for choosing which witness to keep and nothing else.
        rng ^= rng >> 12;
        rng ^= rng << 25;
        rng ^= rng >> 27;
        rng.wrapping_mul(0x2545F4914F6CDD1D)
    };
    let qset: FxHashSet<&str> = names.iter().map(|s| s.as_str()).collect();
    let mut p = Pairs::open(r1, r2, trim);
    let mut lib = 0usize;
    let (mut h1, mut h2): (FxHashMap<u32, u32>, FxHashMap<u32, u32>) =
        (FxHashMap::default(), FxHashMap::default());
    while let Some((name, s1, s2)) = p.next() {
        lib += 1;
        if qset.contains(name.as_str()) {
            continue;
        }
        h1.clear();
        h2.clear();
        for (h, s) in [(&mut h1, &s1), (&mut h2, &s2)] {
            kmers(s, &mut buf);
            for k in buf.iter() {
                if let Some(qs) = index.get(k) {
                    for q in qs {
                        *h.entry(*q).or_insert(0) += 1;
                    }
                }
            }
        }
        // One row per (query, witnessing PAIR), carrying both mates' shared counts. A pair
        // whose two mates both match the query is a witness with no new neighbourhood to
        // offer, and the caller can see that from the two counts rather than from two rows.
        for q in h1.keys().chain(h2.keys()).copied().collect::<FxHashSet<u32>>() {
            let (c1, c2) = (h1.get(&q).copied().unwrap_or(0), h2.get(&q).copied().unwrap_or(0));
            if (c1.max(c2) as usize) < min_shared {
                continue;
            }
            let qi = q as usize;
            seen_count[qi] += 1;
            // The sequence kept is the OTHER mate's - the one that shares less with the
            // query - because that is the neighbourhood this census exists to look at, and
            // carrying it here means the labelling downstream needs no second pass over
            // the library to find out what the neighbour was.
            let other: &Vec<u8> = if c1 >= c2 { &s2 } else { &s1 };
            if keep[qi].len() < max_wit {
                keep[qi].push((name.clone(), c1, c2, other.clone()));
            } else {
                let j = (next_rand() % seen_count[qi]) as usize;
                if j < max_wit {
                    keep[qi][j] = (name.clone(), c1, c2, other.clone());
                }
            }
        }
    }

    let out = std::io::stdout();
    let mut o = std::io::BufWriter::new(out.lock());
    writeln!(o, "query\twitness\tshared_mate1\tshared_mate2\twitnesses_total\tother_mate\tother_seq")
        .unwrap();
    let mut rows = 0usize;
    for (qi, ws) in keep.iter().enumerate() {
        for (w, c1, c2, seq) in ws {
            rows += 1;
            writeln!(
                o,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                names[qi],
                w,
                c1,
                c2,
                seen_count[qi],
                if c1 >= c2 { 2 } else { 1 },
                String::from_utf8_lossy(seq)
            )
            .unwrap();
        }
    }
    let with = seen_count.iter().filter(|&&c| c > 0).count();
    eprintln!(
        "library pairs {} rows {} queries with at least one witness {} of {} (cap {} per query, reservoir sampled)",
        lib, rows, with, names.len(), max_wit
    );
}
