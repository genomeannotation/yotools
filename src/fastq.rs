use std::old_io::{
    BufferedReader,
    Reader,
};
use std::str::StrExt;
use std::string::String;

use seq::Bases;

#[derive(Debug, PartialEq)]
pub struct Sequence {
    pub header: String,
    pub bases: Bases,
    pub qual: String,
}

pub fn read_fastq<R: Reader>(fastq: &mut BufferedReader<R>) -> Vec<Sequence> {
    let mut seqs = vec!();
    loop {
        // Read the header line
        let header_line = if let Ok(line) = fastq.read_line() { line } else { break; };
        // Read the bases line
        let bases_line = if let Ok(line) = fastq.read_line() { line } else { break; };
        // Read the quality header line
        let _ = if let Ok(line) = fastq.read_line() { line } else { break; };
        // Read the quality line
        let qual_line = if let Ok(line) = fastq.read_line() { line } else { break; };

        // Trim the '+' off the header line
        let header = header_line[1..header_line.len()].trim_right().to_string();

        // Add the sequence
        seqs.push(Sequence {
            header: header,
            bases: Bases::new(bases_line.as_slice().trim_right().to_string()),
            qual: qual_line.as_slice().trim_right().to_string(),
        });
    }

    seqs
}

#[test]
fn test_read_fastq_qual_header() {
    use std::old_io::MemReader;

    let mut fastq = BufferedReader::new(MemReader::new(b"@header\nGATACA\n+header\nAAAAAA".to_vec()));

    let expected_seq = Sequence {
        header: "header".to_string(),
        bases: Bases::new("GATACA".to_string()),
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
        bases: Bases::new("GATACA".to_string()),
        qual: "AAAAAA".to_string(),
    };

    assert_eq!(vec![expected_seq], read_fastq(&mut fastq));
}
