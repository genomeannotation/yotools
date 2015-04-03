#![feature(io)]
#![feature(path)]

extern crate bio;

use std::cmp;
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
    let mut fastq_file = BufferedReader::new(File::open(&Path::new("sorted.fastq")));

    let mut samples_file = BufferedReader::new(File::open(&Path::new("samples")));
    let mut loci_file = BufferedReader::new(File::open(&Path::new("loci")));
    
    // Read fastq file
    let seqs = fastq::read_fastq(&mut fastq_file);

    // Store samples and loci
    let loci: Vec<String> = loci_file.lines().map(|x| x.ok().unwrap().as_slice().trim_right().to_string()).collect();
    let samples: Vec<String> = samples_file.lines().map(|x| x.ok().unwrap().as_slice().trim_right().to_string()).collect();

    // Map samples to maps of loci to piles of reads
    let mut seq_matrix: HashMap<String, HashMap<String, Vec<fastq::Sequence>>> = HashMap::new();

    println!("Building matrix...");

    for loci in &loci {
        let mut sample_map = HashMap::new();
        for sample in &samples {
            // Initialize the loci/sample pile to an empty Vec
            sample_map.insert(sample.clone(), vec!());
        }
        // Add the sample row to the matrix
        seq_matrix.insert(loci.clone(), sample_map);
    }

    // Sort all of the sequences into a matrix
    println!("Filling matrix...");
    for seq in seqs {
        use std::iter::IteratorExt;

        let sample = seq.header.as_slice().split(' ').nth(1).unwrap().split(':').nth(3).unwrap().to_string();
        let loci = seq.header.as_slice().split(' ').nth(3).unwrap().trim_right().to_string();

        seq_matrix[loci][sample].push(seq);
    }

    // Unique-ify sequence piles
    println!("Calling consensus...");
    let mut consensus_file = BufferedWriter::new(File::create(&Path::new("consensus.tsv"))); 
    let mut consensus_matrix = BufferedWriter::new(File::create(&Path::new("yay_nay_matrix.tsv"))); 
    let mut count_matrix = BufferedWriter::new(File::create(&Path::new("counts_matrix.tsv"))); 
    
    let mut fasta1_file = BufferedWriter::new(File::create(&Path::new("1.fasta"))); 
    let mut fasta2_file = BufferedWriter::new(File::create(&Path::new("2.fasta"))); 
    let mut fasta3_file = BufferedWriter::new(File::create(&Path::new("3.fasta"))); 
    let mut fasta4_file = BufferedWriter::new(File::create(&Path::new("4.fasta"))); 

    // Write column labels to matrices
    let samples_row = "\t".to_string() + &samples.connect("\t") + "\n";
    consensus_matrix.write_str(samples_row.as_slice());
    count_matrix.write_str(samples_row.as_slice());

    for (loci, sample_map) in &seq_matrix {
        // Begin rows of matrices
        consensus_matrix.write_str(loci.as_slice());
        count_matrix.write_str(loci.as_slice());

        for (sample, seqs) in sample_map {
            // Output count matrix entry
            count_matrix.write_str(format!("\t{}", seqs.len()).as_slice());

            if seqs.len() < 42 {
                consensus_matrix.write_str("\t0");
                continue;
            } 

            let mut seqs_map: HashMap<bases::Bases, u32> = HashMap::new();
            for seq in seqs {
                match seqs_map.entry(seq.bases.clone()) {
                    Entry::Occupied(mut entry) => { *entry.get_mut() += 1; },
                    Entry::Vacant(entry) => { entry.insert(1); },
                }
            }
            if let Some(most_abundant) = seqs_map.iter().map(|(_, count)| count).max().cloned() {
                let seqs_map: HashMap<bases::Bases, u32> =
                    seqs_map.into_iter()
                        .filter(|&(_, count)| (count as f64)/(most_abundant as f64) > 0.1)
                        .collect();

                let mut seq_counts: Vec<(bases::Bases, u32)> = seqs_map.into_iter().collect();
                seq_counts.sort_by(|&(_, a), &(_, b)| a.cmp(&b));
                seq_counts.reverse();

                consensus_matrix.write_str(format!("\t{}", cmp::min(4, seq_counts.len())).as_slice());

                // Write fastas
                if seq_counts.len() > 0 {
                    let &(ref bases1, count1) = &seq_counts[0];
                    fasta1_file.write_str(format!(">{}_{} {}\n{}\n", loci, sample, count1, bases1.as_string()).as_slice());
                }

                if seq_counts.len() > 1 {
                    let &(ref bases2, count2) = &seq_counts[1];
                    fasta2_file.write_str(format!(">{}_{} {}\n{}\n", loci, sample, count2, bases2.as_string()).as_slice());
                }

                if seq_counts.len() > 2 {
                    let &(ref bases3, count3) = &seq_counts[2];
                    fasta3_file.write_str(format!(">{}_{} {}\n{}\n", loci, sample, count3, bases3.as_string()).as_slice());
                }

                if seq_counts.len() > 3 {
                    let &(ref bases4, count4) = &seq_counts[3];
                    fasta4_file.write_str(format!(">{}_{} {}\n{}\n", loci, sample, count4, bases4.as_string()).as_slice());
                }

                // Write consensus file stuff
                consensus_file.write_str(format!("{}\t{}", loci, sample).as_slice());
                for &(_, count) in &seq_counts {
                    consensus_file.write_str(format!("\t{}", count).as_slice());
                }
                consensus_file.write_str("\n");
                for &(ref bases, _) in &seq_counts {
                    consensus_file.write_str(format!("\n{}", bases.as_string()).as_slice());
                }
                consensus_file.write_str("\n\n");
            } else {
                consensus_matrix.write_str("\t0");
            }
        }
        consensus_matrix.write_str("\n");
        count_matrix.write_str("\n");
    }
}
