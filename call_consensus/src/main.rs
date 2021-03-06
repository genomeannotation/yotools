#![feature(io)]
#![feature(path)]
#![feature(core)]
#![feature(convert)]

extern crate bio;

use std::cmp;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::fs::{create_dir, File, OpenOptions};
use std::io::{
    BufRead,
    BufReader,
    BufWriter,
    Write,
};
use std::path::Path;

use bio::bases;
use bio::fastq;

fn main() {
    let mut fastq_file = BufReader::new(File::open(&Path::new("sorted.fastq")).unwrap());

    let mut samples_file = BufReader::new(File::open(&Path::new("samples")).unwrap());
    let mut loci_file = BufReader::new(File::open(&Path::new("loci")).unwrap());

    // Store samples and loci
    println!("Reading loci and samples...");
    let loci: Vec<String> = loci_file.lines().map(|x| x.ok().unwrap().as_str().trim_right().to_string()).collect();
    let samples: Vec<String> = samples_file.lines().map(|x| x.ok().unwrap().as_str().trim_right().to_string()).collect();
    
    // Read fastq file
    println!("Reading fastq...");
    let seqs = fastq::read_fastq(&mut fastq_file);

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
        let sample = seq.header.as_str().split(' ').nth(1).unwrap().split(':').nth(3).unwrap().to_string();
        let loci = seq.header.as_str().split(' ').nth(3).unwrap().trim_right().to_string();

        seq_matrix.get_mut(&loci).unwrap().get_mut(&sample).unwrap().push(seq);
    }

    // Call consensus
    println!("Calling consensus...");

    let mut consensus_file = BufWriter::new(File::create(&Path::new("consensus.tsv")).unwrap()); 
    let mut consensus_matrix = BufWriter::new(File::create(&Path::new("yay_nay_matrix.tsv")).unwrap()); 
    let mut count_matrix = BufWriter::new(File::create(&Path::new("counts_matrix.tsv")).unwrap()); 
    

    // Write column labels to matrices
    let samples_row = "\t".to_string() + &samples.connect("\t") + "\n";
    consensus_matrix.write_all(samples_row.as_bytes());
    count_matrix.write_all(samples_row.as_bytes());

    for (loci, sample_map) in &seq_matrix {
        // Create the loci's fasta folder
        create_dir(loci);

        // Begin rows of matrices
        consensus_matrix.write_all(loci.as_bytes());
        count_matrix.write_all(loci.as_bytes());

        for sample in &samples {
            let seqs = &sample_map[sample];

            // Output count matrix entry
            count_matrix.write_all(format!("\t{}", seqs.len()).as_bytes());

            if seqs.len() < 42 {
                consensus_matrix.write_all(b"\t0");
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

                consensus_matrix.write_all(format!("\t{}", cmp::min(4, seq_counts.len())).as_bytes());

                // Write fasta
                let mut fasta_file = BufWriter::new(File::create(format!("{}/{}.fasta", loci, sample)).unwrap()); 

                for &(ref bases, count) in &seq_counts {
                    fasta_file.write_all(format!(">{}_{} {}\n{}\n", loci, sample, count, bases.as_string()).as_bytes());
                }

                // Write consensus file stuff
                consensus_file.write_all(format!("{}\t{}", loci, sample).as_bytes());
                for &(_, count) in &seq_counts {
                    consensus_file.write_all(format!("\t{}", count).as_bytes());
                }
                consensus_file.write_all(b"\n");
                for &(ref bases, _) in &seq_counts {
                    consensus_file.write_all(format!("\n{}", bases.as_string()).as_bytes());
                }
                consensus_file.write_all(b"\n\n");
            } else {
                consensus_matrix.write_all(b"\t0");
            }
        }
        consensus_matrix.write_all(b"\n");
        count_matrix.write_all(b"\n");
    }
}
