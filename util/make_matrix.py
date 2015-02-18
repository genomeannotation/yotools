#!/usr/bin/env python

import sys
import glob
import re

# Verify command line arguments
if len(sys.argv) != 2:
    sys.stderr.write("usage: python make_matrix.py <path/to/target/folder>\n")
    sys.stderr.write("Target folder should contain a file with the extension '.fastq',\n")
    sys.stderr.write("a file called 'samples' and a file called 'loci'.\n")
    sys.stderr.write("The latter two files should contain one sample/locus name per line.\n\n")
    sys.exit()

path = sys.argv[1]

# Read samples
samples = []
with open(path + "/samples", 'r') as samplefile:
    for line in samplefile:
        samples.append(line.strip())
if not samples:
    sys.stderr.write("Something went wrong reading " + path
            + "/samples, exiting...\n")
    sys.exit()

# Read loci
loci = []
with open(path + "/loci", 'r') as locusfile:
    for line in locusfile:
        loci.append(line.strip())
if not loci:
    sys.stderr.write("Something went wrong reading " + path
            + "/loci, exiting...\n")
    sys.exit()

# Find .fastq file
fastq_files = glob.glob(path + "/*.fastq")
if not fastq_files:
    sys.stderr.write("Couldn't find any fastq files in " + path
            + ", exiting...\n")
    sys.exit()
fastq_file = fastq_files[0]

# Create a big dictionary of dictionaries to hold counts
#   per sample/locus
sample_dict = {}
for sample in samples:
    sample_dict[sample] = {}
    for locus in loci:
        sample_dict[sample][locus] = 0


# Read .fastq file, count occurrences of each locus and sample
with open(fastq_file, 'r') as fastq:
    line_number = 0
    current_sample = ""
    current_locus = ""
    for line in fastq:
        line_number += 1
        if line_number % 4 == 1:
            # Header line
            for sample in samples:
                if sample + " " in line or sample + "\n" in line:
                    current_sample = sample
                    break
            for locus in loci:
                if locus + " " in line or locus + "\n" in line:
                    current_locus = locus
                    break
            if not current_locus or not current_sample:
                sys.stderr.write("Error: failed to find locus and sample ")
                sys.stderr.write("for line " + str(line_number) + ": ")
                sys.stderr.write(line)
                sys.stderr.write("Skipping...\n")
                continue
            sample_dict[current_sample][current_locus] += 1
            current_sample = ""
            current_locus = ""

## Print matrix
# Open output file
with open(path + "/sample_locus_matrix.tsv", 'w') as outfile:
    # Write header line
    outfile.write("\t" + "\t".join(samples) + "\n")
    for locus in loci:
        counts = []
        for sample in samples:
            counts.append(str(sample_dict[sample][locus]))
        outfile.write(locus + "\t")
        outfile.write("\t".join(counts) + "\n")

