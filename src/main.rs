#![feature(io)]
#![feature(core)]
#![feature(path)]
#![feature(unicode)]
#![feature(std_misc)]

use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::old_io::{
    BufferedReader,
    File,
};

use bases::Bases;
use fastq::Sequence;

mod bases;
mod fastq;

#[allow(dead_code)] // So there aren't warnings on unit tests
fn main() {
    let mut fastq = BufferedReader::new(File::open(&Path::new("sample_data/deoligo_test/test1.fastq")));
    
    let seqs = fastq::read_fastq(&mut fastq);

    println!("Sequences: {:?}", seqs);

    // Read the oligos file

    // Maps oligo name to forward and reverse barcode sequences
    let mut oligos: HashMap<String, (String, String)> = HashMap::new();

    let oligos_file = File::open(&Path::new("sample_data/deoligo_test/test.oligos"));
    let mut oligos_file = BufferedReader::new(oligos_file);

    let mut line_number = 0u32;
    while let Ok(line) = oligos_file.read_line() {
        let line_split: Vec<String> = line.as_slice().split('\t').map(|s| s.to_string()).collect();
        
        // Verify that the line is properly formatted
        if line_split.len() != 3 {
            panic!(
                "Expected 3 columns in oligos file, found {} columns at line {}",
                line_split.len(),
                line_number
            );
        }

        // Add the current oligo
        oligos.insert(line_split[2].clone(), (line_split[0].clone(), line_split[1].clone()));

        // Increment line number
        line_number += 1;
    }

    println!("Oligos: {:?}", oligos);

    // Sort the sequences by oligo

    // Maps oligo name to list of debarcoded seqs
    let mut seqs_sorted: HashMap<String, Vec<Sequence>> = HashMap::new();

    for seq in seqs.iter() {
        let mut bases = seq.bases.clone();
        for (oligo_name, &(ref forward, ref reverse)) in oligos.iter() {
            // Attempt to debarcode the sequence
            let debarcoded =
                bases.debarcode(
                    &Bases::from_str(forward.as_slice()),
                    &Bases::from_str(reverse.as_slice()),
                    0,
                );

            // Check if debarcoding succeeded
            if debarcoded {
                // Build the new sorted sequence
                let sorted_seq = Sequence {
                    header: seq.header.clone(),
                    bases: bases,
                    qual: seq.qual.clone(),
                };

                // Update our sorted sequences map
                match seqs_sorted.entry(oligo_name.clone()) {
                    Entry::Occupied(mut entry) => { entry.get_mut().push(sorted_seq); },
                    Entry::Vacant(entry) => { entry.insert(vec![sorted_seq]); },
                }

                // Done with this sequence, move on
                break;
            }
        }
    }

    println!("Sorted seqs {:?}", seqs_sorted);
}
