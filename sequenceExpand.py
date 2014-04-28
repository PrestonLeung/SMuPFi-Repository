#!/usr/local/bin/python

import re, os, argparse
import sys, math
from Bio import SeqIO
from Bio.Seq import Seq

#========================"Self defined Errors and Exceptions"=======================#
class Error(Exception):
    """Base class for exceptions."""
    pass
    
#Error classes for handling input errors
class InputError (Error):
    """Exception for errors in the input.
        Attributes:
        expr -- the input that caused the error
        msg  -- message giving details of the error
    """
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg    
    def __str__(self):
        return repr(str(self.msg)+str(self.expr))

#=========================="openFastaFile"==========================#
#A subfunction to try and open the given file. If the file cannot be opened, throws an error.

def openFastaFile(filename):
    try:
        OpenFastaFile = open(filename)    
    except IOError:
            print "Error: cannot find {}".format(filename)            
    else:
        record = SeqIO.parse(OpenFastaFile,  "fasta")
        OpenFastaFile.close        
        return record

#==========================="main function"===========================#
#Main function does the parsing of the arguments given by the user
#It will then open the files and store it as a BioSeq record using biopython and then call the 
#sub functions to process.

if __name__ == '__main__':
    
    readCount = 0
    program_function = """
    *****
    Expand ShoRAH popl fasta files to their frequencies.
    
    Ver. 0.01
    

    *****
    """
    
    parser = argparse.ArgumentParser()
    #required arguments
    parser.add_argument('-f', '--file', help='Filename or directory of sequence(s) in FASTA format.',  required = True)
    parser.add_argument('-b',  '--base', help='Expanded to b*frequency number of reads. E.g. 0.10 freq at b = 100 will become 10 reads',  required = True,  type = int)
    
    if len(sys.argv) <= 1:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()
    
    if(os.path.isfile(args.file)):        
        fastaHandle = openFastaFile(args.file)        #Reminder this is a generator. Can only iterate once.
    else:
        raise InputError(args.file, "Invalid file: ")
    if(args.base < 1):
        raise InputError(args.base,  "Needs to be greater than 0. Given: ")
    
    for record in fastaHandle:
        freq = float(re.search(r'.*_(.*)',  record.description,  re.M).group(1))
        expand_limit = int(freq * args.base)
        
        for i in range(0, expand_limit):
            print ">READ_{}".format(readCount)
            print record.seq
            readCount += 1
