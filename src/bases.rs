use std::string::String; 

#[derive(Debug, PartialEq)]
pub struct Bases(String);

impl Bases {
    pub fn from_str(bases: &str) -> Bases {
        Bases(bases.to_string())
    }

    /// Reverse complements the sequence
    pub fn reverse_complement(&mut self) {
        let &mut Bases(ref mut bases) = self;
        *bases = bases.as_slice().chars().rev().map(
            |b| {
                match b {
                    'A' => 'T',
                    'a' => 't',
                    'T' => 'A',
                    't' => 'a',
                    'G' => 'C',
                    'g' => 'c',
                    'C' => 'G',
                    'c' => 'g',
                    _ => unreachable!(),
                }
            }
        ).collect();
    }

    /// Tries to debarcode the sequence
    /// On success returns debarcoded sequence
    /// On failure returns original sequence
    pub fn debarcode(self, diffs_allowed: u16) -> Result<Bases, Bases> {
        Err(self)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Unit tests

#[test]
fn reverse_complement() {
    let mut bases = Bases::from_str("GATACA");
    
    bases.reverse_complement();

    assert_eq!(bases, Bases::from_str("TGTATC"));
}

#[test]
fn debarcode_works_properly() {
    let bases = Bases::from_str("ATAGGATACACTAT");

    let debarcoded = bases.debarcode(0);

    assert_eq!(debarcoded, Ok(Bases::from_str("GATACA")));
}

#[test]
fn debarcode_works_properly_with_diff() {
    let bases = Bases::from_str("ATTGGATACATTAT");

    let debarcoded = bases.debarcode(1);

    assert_eq!(debarcoded, Ok(Bases::from_str("GATACA")));
}

#[test]
fn debarcode_fails_properly() {
    let bases = Bases::from_str("ATTGGATACACTAT");

    let debarcoded = bases.debarcode(0);

    assert_eq!(debarcoded, Err(Bases::from_str("ATTGGATACACTAT")));
}
