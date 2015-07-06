#![feature(convert)]

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{
    BufRead,
    BufReader,
    BufWriter,
    Write,
};
use std::str;

fn main() {
    let args: Vec<String> = std::env::args().take(3).collect();
    if args.len() < 3 {
        println!("usage: hetsep <input_fasta> <output_fasta>");
        return;
    }

    // Read the fasta file
    let input_file = BufReader::new(File::open(args[1].as_str()).unwrap());

    println!("Reading fasta...");
    
    let mut seqs: Vec<(String, String)> = Vec::new();

    let mut header = "".to_string();
    for line in input_file.lines() {
        let line = line.unwrap();

        if line.len() == 0 {
            continue;
        } else if line.as_bytes()[0] == b'>' {
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
                    continue;
                }
            }
            seqs.push(bases.clone());
        } else {
            sorted.insert(header.clone(), vec![bases.clone()]);
        }
    }

    // Build IUPAC map
    let mut iupac_map: HashMap<&[u8], &[u8]> = HashMap::new();
    iupac_map.insert(b"-A", b"A");
    iupac_map.insert(b"-T", b"T");
    iupac_map.insert(b"-G", b"G");
    iupac_map.insert(b"-C", b"C");
    iupac_map.insert(b"AG", b"R");
    iupac_map.insert(b"CT", b"Y");
    iupac_map.insert(b"CG", b"S");
    iupac_map.insert(b"AT", b"W");
    iupac_map.insert(b"GT", b"K");
    iupac_map.insert(b"AC", b"M");
    iupac_map.insert(b"CGT", b"B");
    iupac_map.insert(b"AGT", b"D");
    iupac_map.insert(b"ACT", b"H");
    iupac_map.insert(b"ACG", b"V");

    // Output mixed base
    println!("Writing output...");

    let mut out_fasta = BufWriter::new(File::create(args[2].as_str()).unwrap());

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
            } else if base_variants.contains(&b'N') {
                out_fasta.write_all(b"N").unwrap();
            } else {
                let mut sorted_variants: Vec<u8> = base_variants.iter().cloned().collect();
                sorted_variants.sort();
                match iupac_map.get(sorted_variants.as_slice()) {
                    Some(iupac_code) => {
                        out_fasta.write_all(iupac_code).unwrap();
                    },
                    None => {
                        println!("WARNING: Unrecognized IUPAC combination: [{}]", str::from_utf8(sorted_variants.as_slice()).unwrap());
                        out_fasta.write_all(b"[").unwrap();
                        out_fasta.write_all(sorted_variants.as_slice()).unwrap();
                        out_fasta.write_all(b"]").unwrap();
                    },
                }
            }
        }
        out_fasta.write_all(b"\n").unwrap();
    }
}
