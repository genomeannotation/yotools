use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() {
    // Read the fasta file
    let input_file = BufReader::new(File::open("file.fa").unwrap());
    
    let mut seqs: Vec<(String, String)> = Vec::new();

    let mut header = "".to_string();
    for line in input_file.lines() {
        let line = line.unwrap();

        if line.as_bytes()[0] == b'>' {
            header = line[1..line.len()-1].to_string(); // Remember to skip newline
        } else {
            if header == "" {
                panic!("Empty header or invalid fasta");
            }

            let seq = line[1..line.len()-1].to_string(); // Remember to skip newline
            seqs.push((header.clone(), seq));
        }
    }

    // Sort the sequences
    let mut sorted: HashMap<String, Vec<String>> = HashMap::new();

    for &(ref header, ref bases) in &seqs {
        if header.as_bytes()[header.len() - 2] == b'_' {
            let seqs = sorted.entry(header[0..header.len() - 2].to_string()).or_insert(Vec::new());
            seqs.push(bases.clone());
        }
    }
}
