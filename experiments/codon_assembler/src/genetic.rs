//! Standard genetic code over the 0..4 base encoding.
//!
//! A codon index is `16*b0 + 4*b1 + b2` (0..64); any codon containing N has no
//! index. Amino acids are single-letter ASCII with `*` for stop.

const AA: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

#[inline]
pub fn codon_index(c: &[u8]) -> Option<usize> {
    if c[0] > 3 || c[1] > 3 || c[2] > 3 {
        return None;
    }
    Some(16 * c[0] as usize + 4 * c[1] as usize + c[2] as usize)
}

#[inline]
pub fn aa(idx: usize) -> u8 {
    AA[idx]
}

#[inline]
pub fn is_stop(idx: usize) -> bool {
    AA[idx] == b'*'
}

/// Indices of the 61 sense codons, in index order.
pub fn sense_codons() -> Vec<usize> {
    (0..64).filter(|&c| !is_stop(c)).collect()
}

#[inline]
pub fn base_of(idx: usize, pos: usize) -> u8 {
    ((idx >> (2 * (2 - pos))) & 3) as u8
}

#[inline]
pub fn is_transition(a: u8, b: u8) -> bool {
    // A<->G (0,2) and C<->T (1,3)
    a != b && (a & 1) == (b & 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode;
    #[test]
    fn code_spot_checks() {
        let ix = |s: &[u8]| codon_index(&encode(s)).unwrap();
        assert_eq!(aa(ix(b"ATG")), b'M');
        assert_eq!(aa(ix(b"TGG")), b'W');
        assert!(is_stop(ix(b"TAA")) && is_stop(ix(b"TAG")) && is_stop(ix(b"TGA")));
        assert_eq!(aa(ix(b"GCT")), b'A');
        assert_eq!(aa(ix(b"TTA")), b'L');
        assert_eq!(aa(ix(b"AGA")), b'R');
        assert_eq!(sense_codons().len(), 61);
        assert_eq!(base_of(ix(b"ACG"), 0), 0);
        assert_eq!(base_of(ix(b"ACG"), 2), 2);
    }
}
