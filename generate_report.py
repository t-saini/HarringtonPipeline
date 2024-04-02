#TEMPLATE contains the scaffold for the report readout. It is a
#global variable because the contents are only ever referenced
#and never edited.
TEMPLATE = "Sample {} had a {} mold, {} reads, and had {}% of the reads at position {} had the mutation {}.\n"

def gen_report_contents(patient_info:dict, report_template:str = TEMPLATE)->str:
    #this functin takes in a dictionary that has ints as keys, and the values
    #are lists containing patient information in a particular order
    #the second argument, report_template defaults to the global variable
    #TEMPLATE if no other template is provided.
    #an empty string contents is instantiated prior to the for loop
    contents = ""
    #we extract the keys from our input dictionary 
    #and use that to itterate through the for loop
    keys = patient_info.keys()
    for key in keys:
        #as a key is selected
        #the patient information is stored
        #and each element in patient is
        #saved to a variable descirbing the information
        #used in the final patient report.
        patient = patient_info[key]
        name = patient[0]
        color = patient[1]
        read = patient[2]
        read_percent = patient[3]
        position = patient[4]
        mutation = patient[5]
        #thes variables are used with the format method and use each variable
        #passed and fill out the template
        report_line = report_template.format(name, color, read, read_percent, position, mutation)
        #report line is then concatinated to our earlier instantiated contents variable
        contents += report_line
    #once all of the itterations are complete we output contents which
    #should be one large string
    return contents

def gen_report(contents: str,filename:str = 'report.txt')->None:
    #this function takes a string for file contents and for the name
    #of the report. By default the name of the report is report.txt
    with open(filename, "w") as report:
        #with is used as a context manager to easily
        #handle file opening and closing.
        #all items within in contents is written
        #to the text file report.txt
        report.write(contents)
