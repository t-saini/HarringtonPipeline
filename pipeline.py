from necessary_scripts import parseFastq
import logging
import pandas as pd
import sys
import os
from trim_reads import (identify_sequence, trim_leading_barcode, trim_trailing_read,generate_fastq)
from align_fastq import (index_ref_genome, check_index_success, pull_fastq_files, execute_alignment)
from sam_to_bam import (fetch_sam, convert_to_bam, sort_bam_files, index_bam_files, rm_sam_bam)
from variant_discovery import (pileup, prep_reference, detect_mut)
from generate_report import (gen_report_contents, gen_report)


def _run_trim_reads(parse_fastq,clinical_df:pd.DataFrame)->None:
    #this helper function orchestrates the workflow for trimming
    #the fastq reads
    #a dictionary is created where the key is the patient name
    #and it points to an empty list by default
    fastq_dict = {key:[] for key in clinical_df['Name'].unique()}
    #as we itterate through the parse_fastq generator
    for parsed_fastq in parse_fastq:
        #we assign key details into the variables
        #seq_header, sequence, and read_quality
        seq_header = parsed_fastq[0]
        sequence = parsed_fastq[1]
        read_quality = parsed_fastq[3]
        #if the sequence is an empty we continue the for loop
        #in order to ensure we don't process an empty string
        #and hit an error.
        if sequence == '':
            continue
        #we fetch the barcode and other dataframe
        #elments associated with a given barcode
        barcode, result = identify_sequence(sequence, clinical_df)
        #the patient name is extracted from the 0th index
        #within dataframe result 
        patient_name = result['Name'].iloc[0]
        logging.info(f"Looking for barcode: {barcode} within {seq_header}.")
        #sequencing trimming is executed when the barcode, sequence, and read quality
        #strings are passed.
        trimmed_seq, trimmed_qc = trim_leading_barcode(barcode, sequence, read_quality)
        logging.info(f"Leading barcode and qc string removed.")
        logging.info(f"Removing poor trailing reads from {seq_header}.")
        #the edited strings from the previous function are passed into
        #the trim trailing reads in order to achieve the final sequence
        #and final qc string
        final_seq, final_qc = trim_trailing_read(trimmed_seq, trimmed_qc)
        logging.info(f"Removed poor trailing reads.")
        #the data is packaged into an array
        patient_data = [seq_header, final_seq, '+', final_qc]
        #and then saved into our dictionary from the start
        #using the patient name as the key.
        fastq_dict[patient_name].append(patient_data)
    logging.info("Generating new FASTQs")
    #the dictionary is passed into the function
    #generate function, it returns no data types
    #but creates files into the default location for
    #trimmed fastq files.
    generate_fastq(fastq_dict)

def _run_align_fastq(reference_genome:str)->None:
    #a helper function that is used to orchestrate
    #the functions found in align_fastq
    #the input is a string path to the reference genome
    #index_ref_genome will use bwa to index the reference genome
    index_ref_genome(reference_genome)
    #the next function checks for certain extensions that
    #should have been created if index ref genome
    #was successful
    check_index_success(reference_genome)
    #no inputs are used for pull fastq files because
    #the function has a default parameter
    #and outputs an array of string file paths
    fastq_files = pull_fastq_files()
    #the array and the reference genome are passed
    #as arguments together, and a for loop will
    #align each fastq read to the reference
    execute_alignment(reference_genome, fastq_files)

def _run_sam_to_bam()->list:
    #fetch_sam does not need an input as it
    #will by default look at the default
    #outputs an array with strings that are
    #the paths to the sam files
    sam_files = fetch_sam()
    #our array from the previouis function is passed into 
    #convert to bam, the second input is optional, by default
    #the second argument is set to the bam directory path
    convert_to_bam(sam_files)
    #no inputs are required for sort_bam_files
    #as the default arguments point it to the bam directory
    #the sorting is completed using sam tools
    sorted_files = sort_bam_files()
    #the array of sorted files is passed into
    #the index bam files function.
    #the bam files are indexed using samtools
    index_bam_files(sorted_files)
    #using the default parameters, the function
    #rm_sam_bam will remove sam and bam files
    #from the default locations
    rm_sam_bam()
    return sorted_files

def _run_variant_discovery(bam_files:list, reference_genome:str, 
                           clinical_df:pd.DataFrame)->dict:
    #this helper function manages all of the functions required for 
    #variant discovery. The inputs are a list of string paths to sorted bam files,
    #a string path to the reference genome, and a dataframe containing clinical data
    #the final output is a dictionary with ints as keys and lists as values.
    #the dictionary variants discovered is initialized and will contain
    #the final patient details that will be transcribed into a report.
    variants_discovered = {}
    #we iterate through the array bam files and take note of the index
    #associated with the file that is selected
    for index, bam_file in enumerate(bam_files):
        #a dictionary of basebair frequencies per positoin
        #is generated when we pass our input bam file
        bp_frequency = pileup(bam_file)
        #the refrence sequence as a string is extracted from
        #the reference fasta
        reference = prep_reference(reference_genome)
        #by bassing the dictionary bp_frequency and the string reference
        #the the mutation positoin, mutation nucleotide, mutation reads, and mutation
        #frequency is returned as a list.
        variant_details = detect_mut(bp_frequency,reference)
        #the basename of the file is extracted
        #in order to generate the patient name that is
        #associated with the given bam split 
        basename = os.path.basename(bam_file)
        patient_name = basename.split('.')[0]
        #the patient name is then used to exctract the mold color from the
        #clinical data frame, these values are then saved into the array
        #additional_details
        mold_color = clinical_df['Color'][clinical_df['Name']==patient_name].iloc[0]
        #lower() is used to ensure that the molde color is in lowercase letters
        additional_details = [patient_name, mold_color.lower()]
        #by using the extend method additional details will have the
        #elements from variant details tackt on to it.
        additional_details.extend(variant_details)
        #the index that was presented via enumerate is used here to associate
        #the patient details into our dictionary variants discovered
        variants_discovered[index] = additional_details
    return variants_discovered

def _run_generate_report(report_details:dict)->None:
    #this function handles the function required to fully
    #execute the functions needed to generate report.txt
    #the input report details should be a dictionary structured
    #in the following way: 
    #{0:[patient name, mold color, mutant reads, 
    #frequency as an int, position, and mutant bp]}
    report_contents = gen_report_contents(report_details)
    #gen_report_contents will take the dictionary provided and generate a
    #large string using the default template provided within the function.
    #this function can take in a different template if the user desires.
    #gen_report will create the output file, report.txt, this name is also a default
    #value and the second argument of the function gen_report can be provided if the user
    #wants to change the name of the output file.
    gen_report(report_contents)



def main():
    #system arguments are captured and the files are
    #opened using ParseFastQ and Pandas
    arguments = sys.argv[1:]
    parse_fastq = parseFastq.ParseFastQ(arguments[0])
    clinical_df = pd.read_csv(arguments[1], sep='\t')
    ref_genome = arguments[2]
    #main then orchestrates the order of opperations
    #using the helper functions
    logging.info("Trimming FASTQ Files...")
    #past the dataframe and parsed fastq file into the helper
    #function run trim reads
    _run_trim_reads(parse_fastq, clinical_df)
    logging.info("Aligning FASTQ to Reference Genome")
    #the output from run trim reads are 
    #files on the system, the function _run_align_fastq
    #has the default locations for the files
    #so only the reference genome needs to be passed
    _run_align_fastq(ref_genome)
    logging.info("Converting SAM files to BAM files")
    #similarily to _run_align_fastq sam_to_bam
    #has default arguments so there is no need to pass
    #anything into the function, unless the workflow requires
    #something different
    sorted_bam_files = _run_sam_to_bam()
    logging.info("Determining Variant Information")
    #_run_variant_discovery requires the file paths to 
    #the sorted bam files, the reference genome, and clinical data
    #as a dataframe
    variant_details = _run_variant_discovery(sorted_bam_files,
                                             ref_genome,clinical_df)
    logging.info("Writing report.txt")
    #for run_generate_report to work correctly
    #variant_details must a be a dictionary of lists.
    #the output by default will be report.txt
    _run_generate_report(variant_details)



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, 
    format='[%(levelname)s]-%(asctime)s:::%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')
    main()