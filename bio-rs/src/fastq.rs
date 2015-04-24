use std::io;
use bases::Bases;

#[derive(Clone, Debug, PartialEq)]
pub struct Sequence {
    pub header: String,
    pub bases: Bases,
    pub qual: String,
}

impl Sequence {
    pub fn head(&self, n: usize) -> Sequence {
        Sequence {
            header: self.header.clone(),
            bases: self.bases.head(n),
            qual: self.qual[0..n].to_string(),
        }
    }

    pub fn tail(&self, n: usize) -> Sequence {
        Sequence {
            header: self.header.clone(),
            bases: self.bases.tail(n),
            qual: self.qual[self.qual.len()-n..].to_string(),
        }
    }

    pub fn debarcode(&mut self, forward_barcode: &Bases, reverse_barcode: &Bases, diffs_allowed: u16) -> bool {
        let start = forward_barcode.bases.len();
        let end = self.bases.bases.len() - reverse_barcode.bases.len();

        let debarcoded = self.bases.debarcode(forward_barcode, reverse_barcode, diffs_allowed);
        
        if debarcoded {
            self.qual = self.qual[start..end].to_string();
        }

        debarcoded
    }
}

pub fn read_fastq<R: io::Read>(fastq: &mut io::BufReader<R>) -> Vec<Sequence> {
    use std::io::BufRead;

    let mut seqs = vec!();
    loop {
        // Read the four lines making up a fastq sequence
        let mut header_line = String::new();
        let mut bases_line = String::new();
        let mut plus_line = String::new();
        let mut qual_line = String::new();
        
        if fastq.read_line(&mut header_line).is_err() { break; }
        if fastq.read_line(&mut bases_line).is_err() { break; }
        if fastq.read_line(&mut plus_line).is_err() { break; }
        if fastq.read_line(&mut qual_line).is_err() { break; }

        if header_line.len() == 0 {
            break;
        }

        // Trim the '@' off the header line
        let header = header_line[1..header_line.len()].trim_right().to_string();

        // Add the sequence
        seqs.push(Sequence {
            header: header,
            bases: Bases::from_str(bases_line.as_slice().trim_right()),
            qual: qual_line.as_slice().trim_right().to_string(),
        });
    }

    seqs
}

pub fn write_fastq<'a, W, I>(fastq: &mut io::BufWriter<W>, seqs: I) -> io::Result<()>
    where
        W: io::Write,
        I: Iterator<Item=&'a Sequence>,
{
    for seq in seqs {
        try!(write_fastq_seq(fastq, seq));
    }
    Ok(())
}

// Owned version of write_fastq function
// You should use write_fastq if you can, as this consumes the Sequences in the iterator
pub fn write_fastq_owned<W, I>(fastq: &mut io::BufWriter<W>, seqs: I) -> io::Result<()>
    where
        W: io::Write,
        I: Iterator<Item=Sequence>,
{
    for seq in seqs {
        try!(write_fastq_seq(fastq, &seq));
    }
    Ok(())
}

fn write_fastq_seq<W: io::Write>(fastq: &mut io::BufWriter<W>, seq: &Sequence) -> io::Result<()> {
    use std::io::Write;

    try!(fastq.write(format!("@{}\n", seq.header).as_slice().as_bytes()));
    try!(fastq.write(format!("{}\n", seq.bases.as_string()).as_slice().as_bytes()));
    try!(fastq.write(b"+\n"));
    try!(fastq.write(format!("{}\n", seq.qual).as_slice().as_bytes()));
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#[test]
fn test_sequence_head() {
    let seq = Sequence {
        header: "foo".to_string(),
        bases: Bases::from_str("ATGAAAAACAT"),
        qual: "ABCDEFGHIJK".to_string(),
    };

    let seq = seq.head(3);
    assert_eq!(seq.header, "foo");
    assert_eq!(seq.bases, Bases::from_str("ATG"));
    assert_eq!(seq.qual, "ABC".to_string());
}

#[test]
fn test_sequence_tail() {
    let seq = Sequence {
        header: "foo".to_string(),
        bases: Bases::from_str("ATGAAAAACAT"),
        qual: "ABCDEFGHIJK".to_string(),
    };

    let seq = seq.tail(3);
    assert_eq!(seq.header, "foo");
    assert_eq!(seq.bases, Bases::from_str("CAT"));
    assert_eq!(seq.qual, "IJK".to_string());
}

#[test]
fn test_sequence_debarcodes_properly() {
    let mut seq = Sequence {
        header: "foo".to_string(),
        bases: Bases::from_str("ATGAAAAACAT"),
        qual: "ABCDEFGHIJK".to_string(),
    };

    let debarcoded = seq.debarcode(&Bases::from_str("ATG"), &Bases::from_str("ATG"), 0);
    assert!(debarcoded);
    assert_eq!(seq.header, "foo");
    assert_eq!(seq.bases, Bases::from_str("AAAAA"));
    assert_eq!(seq.qual, "DEFGH".to_string());
}

#[test]
fn test_sequence_fails_debarcode_properly() {
    let mut seq = Sequence {
        header: "foo".to_string(),
        bases: Bases::from_str("ATGAAAAACCT"),
        qual: "ABCDEFGHIJK".to_string(),
    };

    let debarcoded = seq.debarcode(&Bases::from_str("ATG"), &Bases::from_str("ATG"), 0);
    assert!(!debarcoded);
    assert_eq!(seq.header, "foo");
    assert_eq!(seq.bases, Bases::from_str("ATGAAAAACCT"));
    assert_eq!(seq.qual, "ABCDEFGHIJK".to_string());
}

#[test]
fn test_read_fastq_qual_header() {
    use std::old_io::MemReader;

    let mut fastq = BufferedReader::new(MemReader::new(b"@header\nGATACA\n+header\nAAAAAA".to_vec()));

    let expected_seq = Sequence {
        header: "header".to_string(),
        bases: Bases::from_str("GATACA"),
        qual: "AAAAAA".to_string(),
    };

    assert_eq!(vec![expected_seq], read_fastq(&mut fastq));
}

#[test]
fn test_read_fastq_no_qual_header() {
    use std::old_io::MemReader;

    let mut fastq = BufferedReader::new(MemReader::new(b"@header\nGATACA\n+\nAAAAAA".to_vec()));

    let expected_seq = Sequence {
        header: "header".to_string(),
        bases: Bases::from_str("GATACA"),
        qual: "AAAAAA".to_string(),
    };

    assert_eq!(vec![expected_seq], read_fastq(&mut fastq));
}

#[test]
fn test_write_fastq() {
    use std::old_io::MemWriter;
    use std::str::from_utf8;

    let mut fastq = BufferedWriter::new(MemWriter::new());

    let seq1 = Sequence {
        header: "foo1".to_string(),
        bases: Bases::from_str("ATG"),
        qual: "AAA".to_string(),
    };

    let seq2 = Sequence {
        header: "foo2".to_string(),
        bases: Bases::from_str("GATACA"),
        qual: "AAAAAA".to_string(),
    };

    let expected_fastq =
        "@foo1\n\
        ATG\n\
        +\n\
        AAA\n\
        @foo2\n\
        GATACA\n\
        +\n\
        AAAAAA\n";

    let result = write_fastq(&mut fastq, vec![seq1, seq2].iter());

    assert!(result.is_ok());

    assert_eq!(from_utf8(fastq.into_inner().into_inner().as_slice()).unwrap(), expected_fastq);
}
