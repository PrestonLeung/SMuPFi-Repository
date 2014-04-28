#!/usr/local/bin/python

# Observed Minus Expected Square tool (OMES)
# ver. 001

# Pre.L


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

#=========================="omes"==========================#

def omes(fastaHash,  length, start_pos):
    
    doubleHash = {}
    singleHash = fillSingleCount(fastaHash, length)
    readSize = len(fastaHash)
    #print readSize
    cov_hash = {}
    for i in range (0, length):        
        for j in range (i+1, length):
            residueHash = {}
            symbolHash_i, symbolHash_j = ({}, {})
            
            omesList = []
            m_infoList = []
            for key in fastaHash.iterkeys():
                doubleKey = (fastaHash[key][i], fastaHash[key][j])
                if(doubleKey in residueHash):                
                    residueHash[doubleKey] += 1
                else:
                    residueHash[doubleKey] = 1
                
                if(fastaHash[key][i] in symbolHash_i):
                    symbolHash_i[fastaHash[key][i]] += 1
                else:
                    symbolHash_i[fastaHash[key][i]] = 1
                
                if(fastaHash[key][j] in symbolHash_j):
                    symbolHash_j[fastaHash[key][j]] += 1 #That is a "j" you are reading. Don't get confused.
                else:
                    symbolHash_j[fastaHash[key][j]] = 1 
                
            for doubleKey in residueHash.iterkeys():
                singleKey1 = doubleKey[0] #key for  position i
                singleKey2 = doubleKey[1] #key for position j
                expected = (singleHash[i][singleKey1] * singleHash[j][singleKey2]) / float(readSize)
                observed = residueHash[doubleKey]                
                omesList.append(math.pow((observed - expected),2))
                
                #jointprob = residueHash[doubleKey]/float(readSize)
                #i_prob = symbolHash_i[singleKey1]/float(readSize)
                #j_prob = symbolHash_j[singleKey2]/float(readSize)                
                #m_infoList.append(jointprob * math.log(jointprob / (j_prob * i_prob)))

            cov = math.fsum(omesList) / readSize
            #mInfo = math.fsum(m_infoList)
            if(cov > 0):
 #               print (i + start_pos), "\t", (j + start_pos), "\t", cov,  "\t", math.fsum(omesList) / len(omesList)
                print (i + start_pos), "\t", (j + start_pos), "\t", cov #,  "\t"#, mInfo
    

#=========================="singleEntropy"==========================#

def singleEntropy(sampleHash, numReads):
    
    listToSum = []
    
    for symbol in sampleHash.iterkeys():
        probabilty = sampleHash[symbol] / float(numReads)
        listToSum.append(probability * math.log(probability))
        
    return -(math.fsum(listToSum))

#=========================="mutualInformation"==========================#

def mutualInformation():
    pass

#=========================="fillSingleCount"==========================#
def fillSingleCount(fastaHash, length):
    
    positionHash = {}
    
    
    for i in range (0, length):
        residueHash = {}
        for key in fastaHash.iterkeys():
            if(fastaHash[key][i] in residueHash):                
                residueHash[fastaHash[key][i]] += 1
            else:
                residueHash[fastaHash[key][i]] = 1

        positionHash[i] = residueHash
    return positionHash
    
#==========================="main function"===========================#
#Main function does the parsing of the arguments given by the user
#It will then open the files and store it as a BioSeq record using biopython and then call the 
#sub functions to process.

if __name__ == '__main__':
    
    bigHash = {}
    seqLength = None
    start_pos = 1
    program_function = """
    *****
    
    Observed Minus Expected Square (OMES) Tool 
    
    Ver. 0.01
    
    (~w~)v
    
    Extra Help:
       - the '-st' option is for user to tell omes.py what's the starting position
         of the given sequence input. Since it only takes in fasta sequences, omes.py
         has no idea what the actual positions of each residue is. Hence defaults the
         first residue to be labelled position 1. If, for example, '-st 200' is entered,
         then the fasta sequence will start at position 200.
    
    *****
    """
    parser = argparse.ArgumentParser()
    #required arguments
    parser.add_argument('-f', '--file', help='Filename or directory of sequence(s) in FASTA format.',  required = True)    
    
    #optional arguments
    parser.add_argument('-st',  '--start_pos',  help = 'Starting position of the fasta file. If left blank defaults at 1.',  type = int)    
    
    if len(sys.argv) <= 1:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()
    
    if(os.path.isfile(args.file)):        
        fastaHandle = openFastaFile(args.file)        #Reminder this is a generator. Can only iterate once.
    else:
        raise InputError(args.file, "Invalid file: ")
    if(args.start_pos):
        start_pos = args.start_pos
    
    for record in fastaHandle:
        bigHash[record.description] =list(record.seq)
        if(not seqLength):
            seqLength = len(record.seq)
    
    omes(bigHash, seqLength, start_pos)
    
