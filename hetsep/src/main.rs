use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{
    BufRead,
    BufReader,
    BufWriter,
    Write,
};

fn main() {
    // Read the fasta file
    let input_file = BufReader::new(File::open("file.fa").unwrap());

    println!("Reading fasta...");
    
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
    println!("Sorting sequences...");

    let mut sorted: HashMap<String, Vec<String>> = HashMap::new();

    for &(ref header, ref bases) in &seqs {
        if header.as_bytes()[header.len() - 2] == b'_' {
            let seqs = sorted.entry(header[0..header.len() - 2].to_string()).or_insert(Vec::new());
            // Make sure any new sequences are the right length
            if let Some(seq) = seqs.get(0) {
                if seq.len() != bases.len() {
                    println!("WARNING: trashing sequence with different length: {}", header);
                    break;
                }
            }
            seqs.push(bases.clone());
        }
    }

    // Output mixed base
    println!("Writing output...");

    let mut out_fasta = BufWriter::new(File::create("output.fa").unwrap());

    for (header, het_seqs) in &sorted {
        write!(&mut out_fasta, ">{}\n", header).unwrap();
        for i in (0..het_seqs[0].len()) {
            let mut base_variants = HashSet::new();
            for seq in het_seqs {
                base_variants.insert(seq.as_bytes()[i]);
            }
            if base_variants.len() == 1 {
                let base = base_variants.iter().next().unwrap();
                out_fasta.write_all(&[*base]).unwrap();
            } else {
                out_fasta.write_all(b"[").unwrap();
                for (i, base) in base_variants.iter().enumerate() {
                    if i != 0 {
                        out_fasta.write_all(b", ").unwrap();
                    }
                    out_fasta.write_all(&[*base]).unwrap();
                }
                out_fasta.write_all(b"]").unwrap();
            }
        }
        out_fasta.write_all(b"\n").unwrap();
    }
}
