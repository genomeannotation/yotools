#![feature(io)]
#![feature(path)]

extern crate bio;

use std::collections::HashMap;
use std::old_io::{
    BufferedReader,
    BufferedWriter,
    File,
};

use bio::bases;
use bio::fastq;

fn main() {
    let mut fastq_file = BufferedReader::new(File::open(&Path::new("sorted.fastq")));

    let mut samples_file = BufferedReader::new(File::open(&Path::new("samples")));
    let mut loci_file = BufferedReader::new(File::open(&Path::new("loci")));
    
    // Read fastq file
    let seqs = fastq::read_fastq(&mut fastq_file);

    // Map samples to maps of loci to piles of reads
    let mut seq_matrix: HashMap<String, HashMap<String, Vec<fastq::Sequence>>> = HashMap::new();

    for sample in samples_file.lines() {
        if let Ok(sample) = sample {
            let mut sample_map = HashMap::new();
            for loci in loci_file.lines() {
                if let Ok(loci) = loci {
                    // Initialize the sample/loci pile to an empty Vec
                    sample_map.insert(loci, vec!());
                }
            }
            // Add the sample row to the matrix
            seq_matrix.insert(sample, sample_map);
        }
    }

    // Sort all of the sequences into a matrix
    for seq in seqs {
        use std::iter::IteratorExt;

        let sample = seq.header.as_slice().split(' ').nth(1).unwrap().split(':').nth(3).unwrap().to_string();
        let loci = seq.header.as_slice().split(' ').nth(3).unwrap().to_string();

        seq_matrix[sample][loci].push(seq);
    }

    //let mut matrix_file = BufferedWriter::new(File::create(&Path::new("matrix.tsv"))); 

    //fastq::write_fastq_owned(&mut sorted_fastq, joined_seqs_iter);
}
