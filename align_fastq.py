import subprocess
import logging
import sys
import os

#the default FASTQ directory is trimmed_fastq
#this enhances readability and maintainability by providing 
#a central reference point.
FASTQ_DIR_DEFAULT = 'trimmed_fastq'
SAM_DIR_DEFAULT = 'sam_files'

def index_ref_genome(reference_genome: str)->None:
    #takes in a string that is the path to the reference genome
    logging.info(f'Indexing refence genome: {reference_genome}')
    try:
        #subprocess is used to kick off bwa genome indexing
        #try and except statements are used in order to capture
        #any failures that may occur
        subprocess.run(['bwa', 'index', f'{reference_genome}'], check=True)
    except subprocess.CalledProcessError as error:
        #notify the user that bwa failed via logger
        #capture the error and display it in logging.warning
        #if the exception is the CalledProcessError
        logging.warning("An error caused BWA to fail!!!")
        logging.warning(error)
        sys.exit(1)

def check_index_success(reference_genome: str)->None:
    #the function takes in a string that is the path to the
    #reference genome that was used with BWA
    expected_extensions = ['.amb','.ann','.bwt','.pac','.sa']
    logging.info("Ensuring the correct files were created from bwa")
    #The expected file extensions post bwa index are sought out
    for expected_extension in expected_extensions:
        file = reference_genome+expected_extension
        #os,path.exists determines if the file exists
        if os.path.exists(file) is False:
            #if the file does not exist
            #a warning is given for the user aong with highlighting
            #what the missing file is. The script then exits.
            logging.warning(f"{file} is missing!")
            logging.warnig("Exiting workflow!!!")
            sys.exit(1)
    #ootherwise it prompts the user everything ran well.
    logging.info("All extensions found, moving forward!")

def pull_fastq_files(fastq_dir:str=FASTQ_DIR_DEFAULT)->list:
    #function takes in a path to trimmed fastqs
    #the default is trimmed_fastq
    logging.info(f"Pulling all trimmed FASTQ files from {fastq_dir}")
    #using os.listdir all file names are saved to a list
    #if the file ends in .fastq and will be outputed
    fastq_files = [fastq_dir+'/'+fastq for fastq in os.listdir(fastq_dir) if fastq.endswith('.fastq')]
    #the number of files are reported along with the name of the directory
    logging.info(f"Found {len(fastq_files)} files within {fastq_dir}")
    return fastq_files

def execute_alignment(reference_genome:str, fastq_files:list, sam_dir:str=SAM_DIR_DEFAULT)->None:
    #function takes a path to the refrence genome and a list of file paths
    #before starting, the sam_files directory is sought out
    logging.info(f'Looking for {sam_dir} directory...')
    if os.path.exists(f'{sam_dir}') is False:
        #if the directory does not exist it is created
        logging.info(f'Did not find {sam_dir} directory')
        logging.info(f'Creating directory {sam_dir}')
        os.mkdir(sam_dir)
    logging.info(f'{sam_dir} directory found...')
    #with the sam_files directory present
    #each fastq file in the list is itterated through
    for fastq in fastq_files:
        #the patient name is extracted by using the split method
        #and then the 0th index is selected
        file_basename = os.path.basename(fastq)
        name = file_basename.split('_')[0]
        logging.info(f"Running BWA on {fastq}")
        #a try except block is used once again since subprocess
        #is being used, and this seems to be the best form
        #of error hanedling.
        try:
            #bwa mem is ran on every file in fastq_files
            #with the results stored in a sam file within the directory
            #the ouput file is created and is ready to be written
            output_file = open(f'{sam_dir}/{name}.sam','w')
            #stdout is used here because if > output_name is used within the list
            #subprocess.run fails to execute the command.
            subprocess.run(['bwa', 'mem' , reference_genome, 
                            fastq], stdout=output_file, check=True)
            output_file.close()
        #if the error that occurs is a Called Process error
        except subprocess.CalledProcessError as error:
            #the error will be displayed as part of the logger warning
            #and the script will exit
            logging.warning(error)
            sys.exit(1)
    #otherwise the user is notified once all of the flies have been aligned.
    logging.info("BWA alignment complete!!!")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, 
    format='[%(levelname)s]-%(asctime)s:::%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')
    ref_genome = 'ref_genome/dgorgon_reference.fa'
    index_ref_genome(ref_genome)
    check_index_success(ref_genome)
    f_files = pull_fastq_files()
    execute_alignment(ref_genome,f_files)
