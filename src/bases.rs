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
    pub fn debarcode(mut self, forward_barcode: &Bases, reverse_barcode: &Bases, diffs_allowed: u16)
        -> Result<Bases, Bases>
    {
        let mut diff_count = 0;

        // Make sure the sequence is long enough to debarcode with these barcodes
        if self.bases.len() < forward_barcode.bases.len() + reverse_barcode.bases.len() {
            return Err(self);
        }

        // Compare forward
        for (base, barcode_base) in self.bases.iter().zip(forward_barcode.bases.iter()) {
            if *base != *barcode_base {
                diff_count += 1;
            }
            if diff_count > diffs_allowed {
                break;
            }
        }

        if diff_count > diffs_allowed {
            return Err(self);
        }

        // Reset diff count
        diff_count = 0;

        // Compare reverse
        for (base, barcode_base) in self.bases.iter().rev().zip(reverse_barcode.bases.iter()) {
            if base.complement() != *barcode_base {
                diff_count += 1;
            }
            if diff_count > diffs_allowed {
                break;
            }
        }

        if diff_count > diffs_allowed {
            return Err(self);
        }

        let start = forward_barcode.bases.len();
        let end = self.bases.len()-reverse_barcode.bases.len();
        self.bases = self.bases[start..end].to_vec();
        Ok(self)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Unit tests

#[test]
fn bases_from_str() {
    Bases::from_str("AaTtGgCcNn");
}

#[test]
fn reverse_complement() {
    let mut bases = Bases::from_str("GATACA");
    
    bases.reverse_complement();

    assert_eq!(bases, Bases::from_str("TGTATC"));
}

#[test]
fn debarcode_works_properly() {
    let bases = Bases::from_str("ATAGGATACACTAT");

    let ref forward = Bases::from_str("ATAG");
    let ref reverse = Bases::from_str("ATAG");
    let debarcoded = bases.debarcode(forward, reverse, 0);

    assert_eq!(debarcoded, Ok(Bases::from_str("GATACA")));
}

#[test]
fn debarcode_works_properly_with_diff() {
    let bases = Bases::from_str("ATTGGATACATTAT");

    let ref forward = Bases::from_str("ATAG");
    let ref reverse = Bases::from_str("ATAG");
    let debarcoded = bases.debarcode(forward, reverse, 1);

    assert_eq!(debarcoded, Ok(Bases::from_str("GATACA")));
}

#[test]
fn debarcode_fails_properly() {
    let bases = Bases::from_str("ATTGGATACACTAT");

    let ref forward = Bases::from_str("ATAG");
    let ref reverse = Bases::from_str("ATAG");
    let debarcoded = bases.debarcode(forward, reverse, 0);

    assert_eq!(debarcoded, Err(Bases::from_str("ATTGGATACACTAT")));
}
