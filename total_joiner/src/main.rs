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
    use ProcessedSequence::*;

    let oligo_size = 10;

    let mut forward_fastq = BufferedReader::new(File::open(&Path::new("forward.fastq")));
    let mut reverse_fastq = BufferedReader::new(File::open(&Path::new("reverse.fastq")));
    
    let forward_seqs = fastq::read_fastq(&mut forward_fastq);
    let reverse_seqs = fastq::read_fastq(&mut reverse_fastq);

    let maybe_joined_seqs =
        forward_seqs.into_iter().zip(reverse_seqs.into_iter()).map(|(forward_seq, reverse_seq)| {
            // Get forward oligo from forward sequence
            let forward_oligo = forward_seq.bases.head(oligo_size);

            // Get reverse oligo from the reverse sequence
            let mut reverse_oligo = reverse_seq.bases.head(oligo_size);
            reverse_oligo.reverse_complement();

            match forward_seq.bases.bases.find(reverse_oligo.bases.as_slice()) { 
                Some(index) => {
                    let merged_len = index + reverse_oligo.len();
                    Joined(fastq::Sequence {
                        header: forward_seq.header,
                        bases: forward_seq.bases.head(merged_len),
                        qual: forward_seq.qual[0..merged_len].to_string(),
                    })
                },
                None => Unjoined(forward_seq, reverse_seq),
            }
        });

    let maybe_joined_seqs: Vec<ProcessedSequence> = maybe_joined_seqs.collect();

    {
        let joined_seqs = maybe_joined_seqs.iter().filter_map(|s| s.clone().map_joined());
        
        // Write joined fastq
        let mut joined_fastq = BufferedWriter::new(File::create(&Path::new("joined.fastq"))); 
        fastq::write_fastq_owned(&mut joined_fastq, joined_seqs).ok()
            .expect("Failed to write joined fastq file");
    }

    // Prepare unjoined iterators
    let unjoined_seqs: Vec<(fastq::Sequence, fastq::Sequence)> = maybe_joined_seqs.into_iter().filter_map(|s| s.map_unjoined()).collect();

    let unjoined_forward_seqs = unjoined_seqs.iter().map(|&(ref f, _)| f);
    let unjoined_reverse_seqs = unjoined_seqs.iter().map(|&(_, ref r)| r);

    // Write unjoined forward fastq
    let mut unjoined_forward_fastq = BufferedWriter::new(File::create(&Path::new("unjoined_forward.fastq"))); 
    fastq::write_fastq(&mut unjoined_forward_fastq, unjoined_forward_seqs).ok()
        .expect("Failed to write unjoined forward fastq file");

    // Write unjoined reverse fastq
    let mut unjoined_reverse_fastq = BufferedWriter::new(File::create(&Path::new("unjoined_reverse.fastq"))); 
    fastq::write_fastq(&mut unjoined_reverse_fastq, unjoined_reverse_seqs).ok()
        .expect("Failed to write unjoined reverse fastq file");
}

#[derive(Clone)]
enum ProcessedSequence {
    Joined(fastq::Sequence),
    Unjoined(fastq::Sequence, fastq::Sequence),
}

impl ProcessedSequence {
    fn map_joined(self) -> Option<fastq::Sequence> {
        use ProcessedSequence::*;

        if let Joined(sequence) = self {
            Some(sequence)
        } else {
            None
        }
    }

    fn map_unjoined(self) -> Option<(fastq::Sequence, fastq::Sequence)> {
        use ProcessedSequence::*;
        
        if let Unjoined(forward, reverse) = self {
            Some((forward, reverse))
        } else {
            None
        }
    }
}
