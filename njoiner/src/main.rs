extern crate bio;

use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::old_io::{
    BufferedReader,
    BufferedWriter,
    File,
};

use bio::bases;
use bio::fastq;

fn main() {
    let mut forward_fastq = BufferedReader::new(File::open(&Path::new("forward.fastq")));
    let mut reverse_fastq = BufferedReader::new(File::open(&Path::new("reverse.fastq")));
    
    let forward_seqs = fastq::read_fastq(&mut forward_fastq);
    let reverse_seqs = fastq::read_fastq(&mut reverse_fastq);

    let joined_seqs_iter = forward_seqs.into_iter().zip(reverse_seqs.into_iter()).map(|(forward_seq, reverse_seq)| {
        fastq::Sequence {
            header: forward_seq.header,
            bases: forward_seq.bases + &bases::Bases::from_str("NNNNNNNNNN") + &reverse_seq.bases,
            qual: forward_seq.qual + "NNNNNNNNNN" + &reverse_seq.qual,
        }
    });

    let mut sorted_fastq = BufferedWriter::new(File::create(&Path::new("joined.fastq"))); 

    fastq::write_fastq_owned(&mut sorted_fastq, joined_seqs_iter);
}
