#![feature(io)]
#![feature(core)]
#![feature(path)]

use std::old_io::{
    BufferedReader,
    File,
};

mod fastq;

#[allow(dead_code)] // So there aren't warnings on unit tests
fn main() {
    let mut sample_fastq = BufferedReader::new(File::open(&Path::new("sample_data/foo.fastq")));
    
    let seqs = fastq::read_fastq(&mut sample_fastq);

    println!("{:?}", seqs[0]);
}
