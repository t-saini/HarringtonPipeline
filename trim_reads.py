import sys
import pandas as pd
import logging
import os
from typing import Tuple

def identify_sequence(sequence:str,ref_df:pd.DataFrame) -> Tuple[str, pd.DataFrame]:
    #looking for the barcode sequence and according to the
    #assignment all of the barcodes should be at the start of
    #the sequence
    logging.info("Looking for barcode..")
    barcode = sequence[0:5]
    #we verify if the barcode is correct by looking it up
    #in the dataframe and if it is empty, the script exits
    result = ref_df[ref_df['Barcode']==barcode]
    if result.empty is True:
        logging.warning(f"Did not find barcode: {barcode} in reference dataframe!!!")
        sys.exit(1)
    #the function finishes with returning our proposed barcode
    #whch has been verified and the associated results from the dataframe
    logging.info(f"Found Barcode {barcode} and associated results")
    return barcode, result

def trim_leading_barcode(barcode:str,sequence:str,qc_string:str)->Tuple[str,str]:
    #use the proposed barcode with the find method
    #to determine the starting index of the barcode
    barcode_index = sequence.find(barcode)
    #if the barcode is not found and a negative 1 is returned,
    #exit the script, if the barcode index is found but it
    #does not start at 0 we also exit out of the script
    #the barcode should always start at 0, if it does not
    #something has gone wrong
    if barcode_index == -1 or barcode_index != 0:
        logging.warning(f"Could not find barcode! Index was : {barcode_index}")
        sys.exit(1)
    #the ending index of the barcode is calculated by
    #taking the length of the barcode and adding the starting
    #index to that value
    barcode_end = barcode_index + len(barcode)
    logging.info("Trimming sequence and QC string...")
    #using string manipulation, we only save the portion of
    #the sequence and qc string that starts with barcode end
    trimmed_seq = sequence[barcode_end:]
    trimmed_qc = qc_string[barcode_end:]
    logging.info("Successfully trimmed sequence and QC string!")
    return trimmed_seq, trimmed_qc

def trim_trailing_read(sequence:str,qc_string:str)->Tuple[str,str]:
    #a dictionary with all possible combinations of poor quality scoring is created
    remove_conditions = {'FF':None,'FD':None,'DD':None,'DF':None}
    #for each of these remove conditions we itterate through
    #if nothing is found and the find method returns -1, we continue
    #in the search.
    for remove_condition in remove_conditions:
        bad_index = qc_string.find(remove_condition)
        if bad_index == -1:
            continue
        remove_conditions[remove_condition] = bad_index
    #once all possible combinations have been found, the minimal value will
    #be stored in remove index, the value will always be an int, as our condition
    #ignores None types
    #min is being used because the lowest value will be the most upstream
    #and according to the assignment if any dual combination is present
    #all parts of the sequence downstream must be removed.
    remove_index = min(value for value in remove_conditions.values() if value is not None)
    #using the index for both the sequence and qc string indicies
    #from the 0th index up to the removal index will be stored and returned
    trimmed_sequence = sequence[:remove_index]
    trimmed_qc = qc_string[:remove_index]
    return trimmed_sequence, trimmed_qc

def generate_fastq(fastq_dict:dict)->None:
    #check to see if the direcotry trimmed_fastq exists
    #if it does, do nothing, if does not create the dir
    logging.info("Looking for directory trimmed_fastq")
    if os.path.exists('trimmed_fastq') is False:
        logging.info('Did not find trimmed_fastq directory')
        logging.info('Creating directory trimmed_fastq')
        os.mkdir('trimmed_fastq')
    #the names of the patients are pulled from our fastq dict
    #and itterated, for ever name, a new file is created and
    #saved to the directory trimmed_fastq
    names = fastq_dict.keys()
    for name in names:
        fastq_data = fastq_dict[name]
        fname = f"{name}_trimmed.fastq"
        logging.info(f"Writing {fname}")
        with open(f'trimmed_fastq/{fname}','w') as file:
            #because my dictionary is contains list of lists
            #the join method is able to turn the list into a string
            #and every , in the list is replaced with a new line
            #a newline is also manually added to ensure that the 
            #FASTQ is propperly readable.
            for data in fastq_data:
                file.write('\n'.join(data)+'\n')