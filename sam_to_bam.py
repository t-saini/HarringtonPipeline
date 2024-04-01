import os
import logging
import sys
import subprocess
import shutil

#the default directories for the SAM and BAM
#files are made as global variables for easy
#access if they ever need to be modified in the future
#and the information here is refrenced frequently
SAM_FILE_DEFAULT = 'sam_files'
BAM_FILE_DEFAULT = 'bam_files'

def fetch_sam(sam_dir:str=SAM_FILE_DEFAULT)->list:
    #function generates a list of file paths
    #from a given directory path, the default path
    #leads to sam_files
    logging.info(f"Fetching SAM files from directory {sam_dir}")
    samfiles = [sam_dir+'/'+sam for sam in os.listdir(sam_dir) if sam.endswith('.sam')]
    num_files = len(samfiles)
    if num_files == 0:
        logging.warning(f"No SAM files found in directory {sam_dir}")
        logging.warning("Exiting!")
        sys.exit(1)
    logging.info(f"Found {num_files} files within directory {sam_dir}")
    return samfiles

def convert_to_bam(sam_files:str, bam_dir:str=BAM_FILE_DEFAULT)->None:
    #takes in a list of sam file paths, and a string with the bam file
    #output directory
    #before starting, the sam_files directory is sought out
    logging.info(f'Looking for {bam_dir} directory...')
    if os.path.exists(f'{bam_dir}') is False:
        #if the directory does not exist it is created
        logging.info(f'Did not find {bam_dir} directory')
        logging.info(f'Creating directory {bam_dir}')
        os.mkdir(bam_dir)
    logging.info(f'{bam_dir} directory found...')
    #itterate through the list of sam file paths
    for sam_file in sam_files:
        #extract the patient name from within the sam file
        file_basename = os.path.basename(sam_file)
        name = file_basename.split('.')[0]
        #generate the bam file name 
        bam_file = f'{name}.bam'
        #try and except are used once again here to capture
        #any errors.
        logging.info(f"Converting {sam_file} to {bam_file}")
        try:
            with open(f'{bam_dir}/{bam_file}','w') as output:
                #stdout is used here because if > output_name is used within the list
                #subprocess.run fails to execute the command.
                subprocess.run(['samtools', 'view', '-bS', sam_file]
                               , stdout=output, check=True)
        #if the error that occurs is a Called Process error
        except subprocess.CalledProcessError as error:
            #the error will be displayed as part of the logger warning
            #and the script will exit
            logging.warning(error)
            logging.warning("Exiting...")
            sys.exit(1)
    #otherwise the user is notified once all of the flies have been aligned.

def sort_bam_files(bam_dir:str=BAM_FILE_DEFAULT)->list:
    #this function takes a string path to the bam file directory
    #by default it points to the default directory bam_files
    sorted_bam_files = []
    bam_files = [bam_dir+'/'+bam for bam in os.listdir(bam_dir) if bam.endswith('.bam')]
    #pull files from directory bam_files
    for bam_file in bam_files:
        #extract the patient name from within the sam file
        file_basename = os.path.basename(bam_file)
        name = file_basename.split('.')[0]
        #generate the sorted bam file name 
        sorted_bam = f'{bam_dir}/{name}.sorted.bam'
        #try and except are used once again here to capture any errors.
        logging.info(f"Converting {bam_file} to {sorted_bam}")
        #expecting a potential error using -o but we will see.
        try:
            subprocess.run(["samtools", "sort", "-m", "100M", "-o", sorted_bam, bam_file])
            sorted_bam_files.append(sorted_bam)
        except subprocess.CalledProcessError as error:
            #the error will be displayed as part of the logger warning
            #and the script will exit
            logging.warning(error)
            logging.warning("Exiting...")
            sys.exit(1)
    return sorted_bam_files


def index_bam_files(sorted_files:list)->None:
    #for each file path in the array
    logging.info("Starting BAM file indexing...")
    #itterate through and pass the file path to subprocess
    for sorted_file in sorted_files:
        logging.info(f"Indexing {sorted_file}")
        try:
            #use try and and except to capture any errors
            #that may occur while running samtools
            subprocess.run(["samtools", "index", sorted_file])
        except subprocess.CalledProcessError as error:
            #if an error does occur it will be captured
            #and the error will be displayed.
            logging.warning(error)
            logging.warning("Exiting...")
            sys.exit(1)

def rm_sam_bam(sam_dir:str=SAM_FILE_DEFAULT,bam_dir:str=BAM_FILE_DEFAULT)->None:
    #this function is used to remove the initial SAM files and sam_files directory
    logging.info(f"Removing SAM files and {sam_dir}")
    try:
        #try and except are used on the off chance an error occurs
        shutil.rmtree(sam_dir)
    except Exception as e:
        #if an error occur output the error through logging
        logging.error(f"An error occured while removing {sam_dir}::\n{e}")
        logging.warning("Exiting!!!")
        sys.exit(1)
    logging.info(f"Successfully removed directory {sam_dir}")
    #once the sam files have been removed inform the user
    #then being removing the bam files
    logging.info(f"Removing BAM files within {bam_dir}")
    bam_files = os.listdir(bam_dir)
    for bam_file in bam_files:
        if bam_file.endswidth('.bam'):
            #the extension .bam is explictly used
            #in order to avoid deleting the sorted bam files
            bam_path = f'{bam_dir}/{bam_file}'
            try:
                #try and except are used once agai for error checking
                logging.info(f'Attempting to remove {bam_path}')
                os.remove(bam_path)
                logging.info(f'Successfully removed {bam_path}')
                #the user will be prompted on the start and the success of
                #every deletion attempt
            except Exception as e:
                #if any error occurs, the bam path will be listed along
                #with the error message
                logging.error(f"An error occured while removing {bam_path}::\n{e}")
                logging.warning("Exiting!!!")
                sys.exit(1)

if __name__ == "__main__":
    sam_files = fetch_sam()
    convert_to_bam(sam_files)
    sorted_files = sort_bam_files()
    index_bam_files(sorted_files)