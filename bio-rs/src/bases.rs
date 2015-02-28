use std::ops::Add;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Base {
    A,
    T,
    G,
    C,
    N,
}

impl Base {
    pub fn complement(self) -> Base {
        use self::Base::*;

        match self {
            A => T,
            T => A,
            G => C,
            C => G,
            N => N,
        }
    }

    pub fn to_char(self) -> char {
        use self::Base::*;

        match self {
            A => 'A',
            T => 'T',
            G => 'G',
            C => 'C',
            N => 'N',
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Bases {
    pub bases: Vec<Base>,
}

impl Bases {
    /// Builds a new Bases sequence from a &str
    /// Panics if bases contains anything other than AaTtGgCcNn
    pub fn from_str(bases: &str) -> Bases {
        use self::Base::*;

        let bases = bases.chars().map(
            |b| {
                match b.to_uppercase() {
                    'A' => A,
                    'T' => T,
                    'G' => G,
                    'C' => C,
                    'N' => N,
                    _ => unreachable!(),
                }
            }
        ).collect();

        Bases { bases: bases }
    }

    /// Gets the bases sequence as a String
    pub fn as_string(&self) -> String {
        self.bases.iter().map(|b| b.to_char()).collect()
    }

    /// Returns the length of the bases sequence
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    // TODO look into deprecating head and tail by supporting indexing and slicing on Bases

    /// Creates new Bases from the first n bases
    pub fn head(&self, n: usize) -> Bases {
        Bases { bases: self.bases.iter().take(n).map(|b| b.clone()).collect() }
    }

    /// Creates new Bases from the last n bases
    pub fn tail(&self, n: usize) -> Bases {
        let mut bases: Vec<Base> = self.bases.iter().rev().take(n).map(|b| b.clone()).collect();
        bases.reverse();
        Bases { bases: bases }
    }

    /// Reverse complements the sequence
    pub fn reverse_complement(&mut self) {
        self.bases.reverse();
        for base in self.bases.iter_mut() {
            *base = base.complement();
        }
    }

    /// Tries to debarcode the sequence
    /// On success returns debarcoded sequence
    /// On failure returns original sequence
    pub fn debarcode(&mut self, forward_barcode: &Bases, reverse_barcode: &Bases, diffs_allowed: u16) -> bool
    {
        let mut diff_count = 0;

        // Make sure the sequence is long enough to debarcode with these barcodes
        if self.bases.len() < forward_barcode.bases.len() + reverse_barcode.bases.len() {
            return false;
        }

        // Compare forward
        for (base, barcode_base) in self.bases.iter().zip(forward_barcode.bases.iter()) {
            if *base != *barcode_base {
                diff_count += 1;
            }
            if diff_count > diffs_allowed {
                return false;
            }
        }

        // Reset diff count
        diff_count = 0;

        // Compare reverse
        for (base, barcode_base) in self.bases.iter().rev().zip(reverse_barcode.bases.iter()) {
            if base.complement() != *barcode_base {
                diff_count += 1;
            }
            if diff_count > diffs_allowed {
                return false;
            }
        }

        let start = forward_barcode.bases.len();
        let end = self.bases.len()-reverse_barcode.bases.len();
        self.bases = self.bases[start..end].to_vec();

        true
    }
}

impl<'a> Add<&'a Bases> for Bases {
    type Output = Bases;

    fn add(self, _rhs: &'a Bases) -> Bases {
        let Bases { bases: mut bases } = self;
        bases.push_all(_rhs.bases.as_slice());
        Bases { bases: bases }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Unit tests

#[test]
fn bases_from_str() {
    use self::Base::*;

    let bases = Bases::from_str("AaTtGgCcNn");
    let expected = Bases { bases: vec![A, A, T, T, G, G, C, C, N, N] };
    assert_eq!(bases, expected);
}

#[test]
fn bases_as_string() {
    use self::Base::*;

    let bases = Bases { bases: vec![A, A, T, T, G, G, C, C, N, N] };
    assert_eq!(bases.as_string().as_slice(), "AATTGGCCNN");
}

#[test]
fn head_bases() {
    let bases = Bases::from_str("GATACT");
    
    assert_eq!(bases.head(3), Bases::from_str("GAT"));
}

#[test]
fn tail_bases() {
    let bases = Bases::from_str("GATACT");
    
    assert_eq!(bases.tail(3), Bases::from_str("ACT"));
}

#[test]
fn reverse_complement() {
    let mut bases = Bases::from_str("GATACA");
    
    bases.reverse_complement();

    assert_eq!(bases, Bases::from_str("TGTATC"));
}

#[test]
fn debarcode_works_properly() {
    let mut bases = Bases::from_str("ATAGGATACACTAT");

    let ref forward = Bases::from_str("ATAG");
    let ref reverse = Bases::from_str("ATAG");
    let debarcoded = bases.debarcode(forward, reverse, 0);

    assert!(debarcoded);
    assert_eq!(bases, Bases::from_str("GATACA"));
}

#[test]
fn debarcode_works_properly_with_diff() {
    let mut bases = Bases::from_str("ATTGGATACATTAT");

    let ref forward = Bases::from_str("ATAG");
    let ref reverse = Bases::from_str("ATAG");
    let debarcoded = bases.debarcode(forward, reverse, 1);

    assert!(debarcoded);
    assert_eq!(bases, Bases::from_str("GATACA"));
}

#[test]
fn debarcode_fails_properly() {
    let mut bases = Bases::from_str("ATTGGATACACTAT");

    let ref forward = Bases::from_str("ATAG");
    let ref reverse = Bases::from_str("ATAG");
    let debarcoded = bases.debarcode(forward, reverse, 0);

    assert!(!debarcoded);
    assert_eq!(bases, Bases::from_str("ATTGGATACACTAT"));
}

#[test]
fn add_bases_ref() {
    let a = Bases::from_str("ATG");
    let b = Bases::from_str("TAG");

    assert_eq!(a + &b, Bases::from_str("ATGTAG"));
}
