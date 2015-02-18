#!/bin/bash

. test/helper_functions.sh

# Clear any previous test results
rm -r temp/

# Run deoligo on the test data
# Usage is "yotools deoligo <in_dir> <out_dir>"
./yotools deoligo sample_data/deoligo_test/ temp/

# Make sure command ran successfully
check_exit_code "run yotools deoligo command"

# Verify temp/ directory created by command
check_directory temp "yotools deoligo sample_data/deoligo_test/ temp/"

# Verify files created
# The input files are called: 
#   sample_data/deoligo_test/test1.fastq
#   sample_data/deoligo_test/test.oligos
# So the program should create these files:
#   temp/test1.deoligoed.fastq
#   temp/deoligo_report.txt
# Check for existence of each one

for x in "test1.deoligoed.fastq" "deoligo_report.txt"
do
	check_file "temp/$x" "yotools deoligo sample_data/deoligo_test/ temp/"
done

# Expected output files are located in sample_data/deoligo_test.expected/
# Check each one against the actual output files
cd sample_data/deoligo_test.expected/
for x in *
do
	check_files_match $x ../../temp/$x
done
cd ../..

