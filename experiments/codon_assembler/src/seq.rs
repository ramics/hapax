//! Sequence encoding and FASTA/FASTQ input.
//!
//! Bases are stored 2-bit style in a `u8` (A=0, C=1, G=2, T=3, N=4). Every
//! downstream stage works on this encoding; the original nucleotide sequence is
//! never replaced by a translation. That is the first non-negotiable of the
//! design: codon space is where overlaps are scored, not where sequence lives.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

pub const N: u8 = 4;

#[inline]
pub fn encode_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => N,
    }
}

#[inline]
pub fn decode_base(b: u8) -> u8 {
    b"ACGTN"[b.min(4) as usize]
}

pub fn encode(s: &[u8]) -> Vec<u8> {
    s.iter().map(|&b| encode_base(b)).collect()
}

pub fn decode(s: &[u8]) -> String {
    s.iter().map(|&b| decode_base(b) as char).collect()
}

/// Reverse complement in the 0..4 encoding (N stays N).
pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| if b < 4 { 3 - b } else { N }).collect()
}

#[derive(Clone, Debug)]
pub struct Record {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Read FASTA or FASTQ (auto-detected from the first character).
pub fn read_seqs(path: &str) -> io::Result<Vec<Record>> {
    let f = File::open(path)?;
    let mut rd = BufReader::new(f);
    let mut out = Vec::new();
    let mut line = String::new();
    // peek first non-empty line
    loop {
        line.clear();
        if rd.read_line(&mut line)? == 0 {
            return Ok(out);
        }
        if !line.trim().is_empty() {
            break;
        }
    }
    if line.starts_with('>') {
        let mut name = header_name(&line[1..]);
        let mut seq = Vec::new();
        loop {
            line.clear();
            if rd.read_line(&mut line)? == 0 {
                break;
            }
            if let Some(h) = line.strip_prefix('>') {
                out.push(Record { name, seq: encode(&seq) });
                name = header_name(h);
                seq = Vec::new();
            } else {
                seq.extend(line.trim_end().as_bytes());
            }
        }
        out.push(Record { name, seq: encode(&seq) });
    } else if line.starts_with('@') {
        let mut header = line.clone();
        loop {
            let mut s = String::new();
            let mut plus = String::new();
            let mut q = String::new();
            rd.read_line(&mut s)?;
            rd.read_line(&mut plus)?;
            rd.read_line(&mut q)?;
            out.push(Record { name: header_name(&header[1..]), seq: encode(s.trim_end().as_bytes()) });
            header.clear();
            if rd.read_line(&mut header)? == 0 || header.trim().is_empty() {
                break;
            }
        }
    } else {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "input is neither FASTA nor FASTQ"));
    }
    Ok(out)
}

fn header_name(h: &str) -> String {
    h.split_whitespace().next().unwrap_or("").to_string()
}

pub fn write_fasta<W: Write>(w: &mut W, name: &str, seq: &[u8]) -> io::Result<()> {
    writeln!(w, ">{}", name)?;
    for chunk in seq.chunks(80) {
        w.write_all(&chunk.iter().map(|&b| decode_base(b)).collect::<Vec<u8>>())?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn revcomp_roundtrip() {
        let s = encode(b"ACGTNAAC");
        assert_eq!(decode(&revcomp(&s)), "GTTNACGT");
        assert_eq!(revcomp(&revcomp(&s)), s);
    }
}
