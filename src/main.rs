#![feature(io)]
#![feature(core)]
#![feature(path)]

use std::old_io::{
    BufferedReader,
    File,
};

use fastq::read_fastq;

mod fastq;

fn main() {
    let mut sample_fastq = BufferedReader::new(File::open(&Path::new("sample_data/foo.fastq")));
    
    let seqs = read_fastq(&mut sample_fastq);

    println!("{:?}", seqs[0]);
}
