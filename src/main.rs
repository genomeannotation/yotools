#![feature(io)]
#![feature(core)]
#![feature(path)]
#![feature(unicode)]

use std::old_io::{
    BufferedReader,
    File,
};

mod bases;
mod fastq;

#[allow(dead_code)] // So there aren't warnings on unit tests
fn main() {
    let mut fastq = BufferedReader::new(File::open(&Path::new("sample_data/deoligo_test/test1.fastq")));
    
    let seqs = fastq::read_fastq(&mut fastq);

    // Read the oligos file

    println!("{:?}", seqs[0]);
}
