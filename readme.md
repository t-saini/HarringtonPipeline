# HarringtonPipeline
# Author: Tanvir Saini

This pipeline automates the analysis of pooled sequencing data to discover variants and generate a comprehensive report. It takes as input a FASTQ file containing pooled sequencing data, a tabulated text file providing clinical information associated with each sample, and a FASTA file serving as the reference sequence. The workflow is orchestrated by the `pipeline.py` script, which efficiently processes the input data and generates a comprehensive report summarizing the findings.

## Workflow Overview:

1. **Trim Reads**: Initial processing of the input FASTQ file to remove adapter sequences, low-quality bases, and other artifacts. This step enhances the quality of the sequencing data and improves downstream analysis accuracy.

2. **Align to Reference Genome**: The trimmed reads are aligned to a reference genome using the Burrows-Wheeler Aligner (BWA). This step ensures that the sequencing reads are correctly mapped to the reference sequence, enabling further analysis.

3. **SAM to BAM Conversion**: The aligned reads are initially stored in the SAM format, which is then converted into the more efficient and indexed BAM format. This conversion simplifies data storage and facilitates faster data retrieval for subsequent analysis steps.

4. **Variant Discovery**: Analysis of the aligned reads to identify genetic variants, such as single nucleotide polymorphisms (SNPs). This step involves comparing the aligned reads to the reference genome to detect variations and assess their frequency within the sample population.

5. **Variant Sequence Analysis**: Further analysis of detected variants to characterize their sequence context, assess their potential impact on gene function, and prioritize variants of clinical significance. This step aids in understanding the biological relevance of the identified variants and their implications for disease risk or treatment response.

6. **Generate Report**: Compilation of the variant analysis results, along with associated clinical information, into a structured report. This report provides a comprehensive summary of the detected variants, their frequency, genomic locations, and any relevant clinical annotations. It serves as a valuable resource for researchers, clinicians, and other stakeholders involved in genomic data interpretation and decision-making.

## Input Files:

- **FASTQ File**: Contains the raw sequencing data obtained from the sample pool.
  
- **Tabulated Text File**: Provides clinical metadata associated with each sample, such as patient identifiers, clinical characteristics, and experimental conditions.
  
- **FASTA File**: Represents the reference genome or target genomic region against which the sequencing reads are aligned.

## Usage:

To run the pipeline, execute the `pipeline.py` script within the week5 directory.

```
python3 pipeline.py "<input_fastq>" "<clinical_data_file>" "<reference_genome_fasta>"
```
If your input files are outside of the week5 directory, provide the full path.

```
python3 pipeline.py "<path>/<to>/<your>/<input_fastq>" "<path>/<to>/<your>/<clinical_data_file>" "<path>/<to>/<your>/<reference_genome_fasta>"
```

# Logging

The pipeline script logs events and errors at different stages of the workflow to provide information about the progress and to capture any issues encountered during execution. Log messages are written to the terminal and can be helpful for troubleshooting and debugging purposes.
