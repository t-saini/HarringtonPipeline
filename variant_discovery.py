import pysam
import logging

def pileup(bamfile:str)->dict:
    #this function takes in a string that is a path
    #to a sorted and indexed bam file.
    #the bam file is passed to pysam.Alignmentfile
    samfile = pysam.AlignmentFile(bamfile, "rb")
    #ntdict is a dictionary that will contian the the positions as keys
    #the values of the dictionaries will be additional dictionaries where
    #the key will be the nucleotide found for that given frequency and 
    #the value will be the frequency of how often that nucleotide appears
    #{0:{'A':1,'C':5,'G':1,'T':2}}
    ntdict = {}
    #by applying the pileup method to our bamfile that has been passed to
    #pysam.AlignmentFile, we can itterate through as if it was a generator function
    for pileupcolumn in samfile.pileup():
        #inform the user what the base is and what position we are currently looking at
        logging.info("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        #pileupcolumn.pos is the position that is currently being looked at
        ntdict[pileupcolumn.pos] = {}
        for pileupread in pileupcolumn.pileups:
            #itterate through the the pileups element that is created from
            #the generator pileup
            #check to see if the selected element has been identified as
            #a deletion or if it is a reference skip
            if not pileupread.is_del and not pileupread.is_refskip:
                #historic print statement this is left for testing/debugging
                #print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                #if the nucleotide we are currenlty on is not seen as a deletion or reference skip, 
                #we add it to the dictionary.
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if base not in ntdict[pileupcolumn.pos]:
                    #if the nucleotide does not exist in ntdict, add it
                    #make the default value 1
                    ntdict[pileupcolumn.pos][base] = 1
                    #once complete we continue with the for loop, and move on to 
                    #the next itteration
                    continue
                ntdict[pileupcolumn.pos][base] +=1
                #if the base pair already exists withing the dicitonary
                #add one to the current count
    #close the samfile to prevent memeory leaks and any issues with the code.
    #and output ntdict.
    samfile.close()
    return ntdict

def prep_reference(file_path:str)->str:
    #this function takes the file path to the reference genome
    #the file is passed to pysam.FastaFile
    ref_genome = pysam.FastaFile(file_path)
    #we are able to extract the name using the refrence
    #method which returns a list so we must explictly pull from the 0th index
    reference_name = ref_genome.references[0]
    #with the name acquired, the referenec sequence can be extracted
    #and returned
    reference_sequence = ref_genome.fetch(reference_name)
    return reference_sequence

def detect_mut(ntdict:dict, reference:str)->list:
    #this function takes in dictionary of dictionaries containing positions (primary keys)
    #of nucleotides (secondary keys) with the frequency of each nucletodie (value)
    #example: {0:{'C':5, 'A':1}}
    pos_keys = ntdict.keys()
    #extract the primary keys and store them as
    #an array so we may itterate through our dict of dicts
    for pos in pos_keys:
        position = ntdict[pos]
        #if the dict only contains 1 nucleotide
        #we ignore it, and move on
        if len(position) <= 1:
            continue
        #the moment we find two nucleotides within
        #one a key:value pair we take the the sum
        #of the frequency values
        nucleotides = position.keys()
        total_reads = sum(position.values())
        #using the reference and position determine the wildtype
        #basepair associated with the position we found in our dict with
        #two basepairs
        wildtype = reference[pos]
        #using list comprehension we itterate through the nucleotides we foudn in our dict
        #but we only identify the mutant basepair if it is not the wild type.
        #since we are working with SNPs the 0th index of this list will always be the mutant
        mutant = [nucleotide for nucleotide in nucleotides if nucleotide != wildtype][0]
        #frequency of occurance is calculated by calling the frequency of the mutant from 
        #our dictionary and dividing it by the total reads we calcualted earlier.
        mutant_reads = position[mutant]
        frequency_occurance = int(round((mutant_reads/total_reads) * 100, 2))
        #python counts starting at zero we have to add one 
        true_position = pos + 1
    #return the number of mutant reads, frequency as an int, position, and mutant bp
    return [mutant_reads, frequency_occurance, true_position, mutant]
