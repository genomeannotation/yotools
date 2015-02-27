#![feature(io)]
#![feature(path)]
#![feature(core)]

extern crate bio;

use std::old_io::{
    BufferedReader,
    BufferedWriter,
    File,
};

use bio::fastq;

use find::Find;

mod find;

fn main() {
    let oligo_size = 10;

    let mut forward_fastq = BufferedReader::new(File::open(&Path::new("forward.fastq")));
    let mut reverse_fastq = BufferedReader::new(File::open(&Path::new("reverse.fastq")));
    
    let forward_seqs = fastq::read_fastq(&mut forward_fastq);
    let reverse_seqs = fastq::read_fastq(&mut reverse_fastq);

    let joined_seqs_iter = forward_seqs.into_iter().zip(reverse_seqs.into_iter()).map(|(forward_seq, reverse_seq)| {
        // Get forward oligo from forward sequence
        let forward_oligo = forward_seq.bases.head(oligo_size);

        // Get reverse oligo from the reverse sequence
        let mut reverse_oligo = reverse_seq.bases.head(oligo_size);
        reverse_oligo.reverse_complement();

        if let Some(index) = forward_seq.bases.bases.find(reverse_oligo.bases.as_slice()) {
            println!("Found oligo at {:?}", index);
        }

        fastq::Sequence {
            header: forward_seq.header,
            bases: forward_oligo,
            qual: forward_seq.qual,
        }
    });

    let mut sorted_fastq = BufferedWriter::new(File::create(&Path::new("joined.fastq"))); 

    fastq::write_fastq_owned(&mut sorted_fastq, joined_seqs_iter).ok().expect("Failed to write fastq file");
}
