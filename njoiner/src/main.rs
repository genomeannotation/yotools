extern crate bio;

use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::old_io::{
    BufferedReader,
    File,
};

use bio::fastq;

fn main() {
    let mut forward_fastq = BufferedReader::new(File::open(&Path::new("forward.fastq")));
    let mut reverse_fastq = BufferedReader::new(File::open(&Path::new("reverse.fastq")));
    
    let forward_seqs = fastq::read_fastq(&mut forward_fastq);
    let reverse_seqs = fastq::read_fastq(&mut reverse_fastq);
}
