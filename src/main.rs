#![feature(io)]
#![feature(core)]
#![feature(path)]

use std::old_io::{
    BufferedReader,
    File,
    Reader,
};
use std::str::StrExt;
use std::string::String;

#[derive(Debug)]
pub struct Sequence {
    pub header: String,
    pub bases: String,
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
            bases: bases_line.as_slice().trim_right().to_string(),
            qual: qual_line.as_slice().trim_right().to_string(),
        });
    }

    seqs
}

fn main() {
    let mut sample_fastq = BufferedReader::new(File::open(&Path::new("sample_data/foo.fastq")));
    
    let seqs = read_fastq(&mut sample_fastq);

    println!("{:?}", seqs[0]);
}
