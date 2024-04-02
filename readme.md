Sequence Alignment Pipeline
Author: Tanvir Saini

This pipeline automates the analysis of pooled sequencing data; including trimming reads, 
aligning them to a reference genome, variant sequence analysis, and generating a report 
summarizing the findings.The workflow is orchestrated by the pipeline.py script, 
which takes in three input files: 
a FASTQ file containing pooled sequencing data, a tabulated text file with clinical information, 
and a FASTA file representing the reference sequence.

## Workflow Overview

1.Data Preprocessing: Initial FASTQ file (hawkins_pooled_sequences.fastq) is processed to
    have the leading barcode trimmed from the 5' end. While on the 3' end poor quality reads and
    No calls (N) are removed when the QC string has two occuring instances of FF, FD, or DD.
    The trimmed results are saved to a new FASTQ in the directory fastqs.
2.Sequence Alignment: Aligns the trimmed FASTQ reads to the provided reference genome 
    using BWA (Burrows-Wheeler Aligner). Successful alignments are stored in .SAM files
    in the directory sam_files.
3.SAM TO BAM File Convertion: Aligned SAM files are converted to BAM format, sorts them, 
    indexes them, and removes intermediate SAM files to conserve disk space. Successful
    outputs are stored in the directory bams.
4.Variant Discovery: Performs variant discovery by analyzing the pileup of reads at each position 
    in the aligned BAM files. It identifies mutations, calculates their frequencies, and determines 
    their positions.
5.Report Generation: A report summarizing variant discovery results with the default name "report.txt".

## Files Overview:

1.parseFastq.py
Description: Python script for parsing FASTQ files, acts similarly to a generator function.

2.trim_reads.py
Description: Python script for trimming DNA barcodes from the FASTQ sequences, and trims 
    the FASTQ based on the QC string. Outputs are stored in new FASTQ files using the following
    convention "<patient name>.fastq" within the fastqs directory.

3.align_fastq.py
Description: Python script that is indexed using the BASH variant of BWA. After which, 
    the trimmed FASTQs are aligned using the BASH varaint of SAMtools and outputs a SAM file
    with the following convention "<patient name>.sam" within the sam_files directory.

4.sam_to_bam.py
Description: Python script for converting SAM files to sorted.bam files. Intermediate files (.sam and .bam)
    are deleted at the end of the convertion along with the sam_files directory. The final files will be named with the following convention
    "<patient name>.sorted.bam" within the bam_files directory.

5.variant_discovery.py
Description: Python script that uses Pysam to analyze pileup data from a sorted and indexed BAM file.
    The results are then compared against the reference sequence provided, and the nucleotide frequcies
    are calculated. The output is a list containig the mutation reads, frequency of occurance, the position
    of the mutation, and the mutated nucleotide.

6.generate_report.py
Description: Python script that constructs a text based report from a dictionary input that contains

7.pipeline.py
Description: Python script to orechestrate the execution of the other Python scripts and generates the final report.

## Requirements

This project requires the following Python packages:

- pandas
- logging
- pysam
- shutil
- sys 
- os 
- typing 
- subprocess

This project requires the following Python scripts:

- parseFastq from the directory necessary_scripts

## How to Run

Execute the following command in the terminal when inside the week5 directory using the provided files and the directory
necessary_scripts must be present. Inputs such as pooled sequence FASTQs and tabulated clinical data text filesshould ideally be in the input directory. 
The reference fasta for the genome of interest should be placed within the ref_genome directory. This is the default file structure.

Example with provided files in the suggested directories:

python3 pipeline.py "inputs/hawkins_pooled_sequences.fastq" "inputs/harrington_clinical_data.txt" "ref_genome/dgorgon_reference.fa"

If your fastq file (`hawkins_pooled_sequenecs.fastq`), clincal data text file (`harrington_clincal_data.txt`), or 
reference fasta (`dgorgon_reference.fa`) are located in a different set of directories, provide the full path to the files 
when executing the command within the week5 directory. 

For example:

python3 pipeline.py "/path/to/your/hawkins_pooled_sequenecs.fastq" "/path/to/your/harrington_clincal_data.txt" "/path/to/your/dgorgon_reference.fa"

## Input
Pooled FASTQ(hawkins_pooled_sequenecs.fastq): FASTQ file with pooled sequences, with barcode sequences
    that should coinside with the barcodes found in the associated Clincal Data File.

Clinical Data File(harrington_clincal_data.txt): Must be a tab-delimited text file, the following
    column names are required:
    Name, Color, Barcode

Reference Genome(dgorgon_reference.fa): A fasta file containing the genome of the refrence organism that the FASTQ
    files will be aligned to.


## Logging

The pipeline script logs events and errors at different stages of the workflow to provide information about the progress 
and to capture any issues encountered during execution. Log messages are written to the console and can be helpful for 
troubleshooting and debugging purposes.

Example:

1. Incorrect key file extension:
    [WARNING]- YYYY-MM-DD H:M:S ::: "Checking file extensions for the inputs given."
    [WARNING]- YYYY-MM-DD H:M:S ::: "<path_to_file> does not have extension <required_extension>"

2. Incorrect barcode provided in clinical data:
    [WARNING]- YYYY-MM-DD H:M:S ::: "Could not find barcode! Index was : <barcode>"
    [WARNING]- YYYY-MM-DD H:M:S ::: "Exiting with error code 1."

3. subprocess.CalledProcessError when running samtools via subprocess:
    [ERROR]- YYYY-MM-DD H:M:S ::: Failed to create an alignment with <sam_file> and <reference_genome>
    [WARNING]- YYYY-MM-DD H:M:S ::: "<subprocess.CallProcessError>"