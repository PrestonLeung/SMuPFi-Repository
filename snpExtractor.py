#!/usr/bin/python2.7
#Ver. 0.01
#Pre.L

import argparse,  re
import os,  sys
from Bio import SeqIO
from Bio.Seq import Seq


#========================"Global Variables"=======================#

OutputFileName = ''
IgnoreFirstLine = False
Reparse = False
WantOutput = False
LogFrequencies = False
IncludeSyn = False

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


#==========================="okFile"===========================#        
#Check if the given file is a file. 1 means yes, 0 means no.

def okFile(file):
    if(os.path.isfile(file)):
        return 1
    else:
        return 0

#==========================="extractLine"===========================#        


def extractLine(file,  column,  delimiter):
    
    alreadyIgnored = False
    lineList = []    
    readFile = open(file)
    
    for line in readFile: 
        line.strip()
        if(len(line)<3):
            continue
        if(IgnoreFirstLine):
            
            if(not alreadyIgnored):
                alreadyIgnored = True
                continue
        splitLine = line.split(delimiter)
        columnList = []        
        for position in column:            
            columnList.append(splitLine[int(position)])            
        lineList.append(columnList)
    readFile.close()    
    return lineList    

#==========================="printData"===========================#
#Writes data to an easy-to-read (for humans) format.

def printData(aList):
    
    if(WantOutput):
        openFile = open(OutputFileName, 'w')        
    else:
        openFile = sys.stdout
    openFile.write( "nt-Position".ljust(20)), openFile.write("aa-Position".ljust(20))
    openFile.write("wtCodon".ljust(10)), openFile.write("wtAA".ljust(10) )
    openFile.write("mtCodon".ljust(10)), openFile.write("mtAA".ljust(10)), openFile.write("Frequency".ljust(20)+ '\n')     
    for item in aList:
        openFile.write(str(item[0]).ljust(20)), openFile.write(str(item[1]).ljust(20))
        openFile.write(str(item[2]).ljust(10)), openFile.write(str(item[3]).ljust(10))
        openFile.write(str(item[4]).ljust(10)), openFile.write(str(item[5]).ljust(10))
        openFile.write(str(item[6]).ljust(20) + '\n')
    openFile.close()    
    
#==========================="printEasy"===========================#
#Prints findings in an easy-to-parse format.

def printEasy(aList):
    
    if(OutputFileName):
        short = re.sub('\..*$',  '',  OutputFileName)
        openFile = open(short + '_EasyOutput.txt', 'w')        
    else:
        openFile = open('EasyOutput.txt',  'w')    
    for item in aList:
        openFile.write(str(item[0]) + "|"), openFile.write(str(item[1])+ "|")
        openFile.write(' '.join(str(item[2]))+ "|"), openFile.write(str(item[3])+ "|")
        openFile.write(' '.join(str(item[4]))+ "|"), openFile.write(str(item[5])+ "|")
        openFile.write(str(item[6])+ '\n')
    openFile.close()   


#==========================="changeTypeThenSort"===========================#
#Changes the fields to their respective types. E.g. Position will change from str to int and
#Frequency values will change from str to float.

def changeTypeThenSort(aList,  ntPosition, freqPosition,  mutPosition,  sortFreq):
    
    if(sortFreq):
        sortPosition = freqPosition
    else:
        sortPosition = ntPosition
    
    for item in aList:
        item[ntPosition] = int(item[ntPosition])
        item[freqPosition] = float(item[freqPosition])
        
        if(re.search('>', item[mutPosition])):
            item[mutPosition] = item[mutPosition].split('>')[1]
        
    
    aList.sort(key=lambda x: x[sortPosition])
    
    return aList
    
#==========================="getFile"===========================#
#Opens a file and parses it as a seqRecord file and returns the record.

def getFile(filename):
    try:
        OpenFastaFile = open(filename)    
    except IOError:
            print "Error: cannot find {}".format(filename)            
    else:
        record = SeqIO.parse(OpenFastaFile,  "fasta")
        OpenFastaFile.close        
        return record    

#==========================="translation"===========================#
#Function does the translation of mutated nucleotides in the respective positions in relation to the 
#user-specified starting codon. E.g. Mutation X at position N is in the codon of amino acid position A
#given the starting codon is S.

def translation(aList,  ref,  startCodon, ntPosition,  freqPosition,  mutPosition,  filter):
    
    refSeq = ''
    megaRecordList = []
    
    for record in ref:
        refSeq = str(record.seq)
    
    codingSeq = refSeq[startCodon-1:len(refSeq)]    
    triplet = re.findall('.?.?.?',  codingSeq) 
    
    for item in aList:
        
        codonRecordList  = []    
        
        if(item[ntPosition] < startCodon):
            continue
        if(filter):
            if(item[freqPosition] >= filter):
                continue
        if(item[freqPosition] == 0.0):
            continue
        if(item[mutPosition] == '-'):
            continue

        aaposition = frameShift(item[ntPosition],  0,  startCodon)
        wtCodon = triplet[aaposition-1]
        r = (item[ntPosition] - (startCodon - 1)) % 3        
        wtAA = str(Seq(wtCodon).translate())
        mtCodon = list(wtCodon)
        
        if(r == 0):
            mtCodon[2] = item[mutPosition]
        elif(r == 2):
            mtCodon[1] = item[mutPosition]
        else:
            mtCodon[0] = item[mutPosition]
        mtCodon = ''.join(mtCodon)        
        mtAA = str(Seq(mtCodon).translate())        
        
        if(IncludeSyn):
            codonRecordList.append(item[ntPosition]) # nt-position
            codonRecordList.append(aaposition) # aa-position
            codonRecordList.append(wtCodon) # reference codon
            codonRecordList.append(wtAA) # reference aa
            codonRecordList.append(mtCodon) # mutant codon
            codonRecordList.append(mtAA) # mutant aa
            if(LogFrequencies):
                codonRecordList.append(math.log10(item[freqPosition] * 100))
            else:
                codonRecordList.append(item[freqPosition]) # frequency of occurrence
            megaRecordList.append(codonRecordList)
        else:
            if(mtAA != wtAA):
                codonRecordList.append(item[ntPosition]) # nt-position
                codonRecordList.append(aaposition) # aa-position
                codonRecordList.append(wtCodon) # reference codon
                codonRecordList.append(wtAA) # reference aa
                codonRecordList.append(mtCodon) # mutant codon
                codonRecordList.append(mtAA) # mutant aa
                if(LogFrequencies):
                    codonRecordList.append(math.log10(item[freqPosition] * 100))
                else:
                    codonRecordList.append(item[freqPosition]) # frequency of occurrence
                megaRecordList.append(codonRecordList)
    
    return megaRecordList

#========================"frameShift=========================#
#This subfunction computes nucleotide positions into respective amino acid position.
#Must provide the shiftValue, for example, if the user wants to read nucleotides with none (0), one or two
#shifted. 
#

def frameShift(position,  shiftValue, startCodon):
        
        codingPos = position - (startCodon - 1)
        shift1 = shiftValue % 3
        shift2 = shiftValue / 3
        q = codingPos / 3
        r = codingPos % 3
        
        if(shift1 == 0):
            if(r == 0):
                return (q - shift2)
            else:   
                return (q - shift2 + 1)
        elif(shift1 == 1):
            if(r == 2):
                return (q - shift2 + 1)
            else:   
                return q - shift2
        else:
                return q - shift2

#==========================="reparse"===========================#
#Function creates a temporary file if -rp is triggered to remove ^M characters coming from
#DOS files.

def reparse(file):
    
    import string
    tempFileName = "awefbiawwoynaw2q3095092w3.txt"
    f = open(file)
    temp = open(tempFileName,  'w')
    
    for line in f:        
        temp.write(line.translate(string.maketrans("\r",  "\n")))
    f.close()    
    temp.close()
    
    return tempFileName

#==========================="main function"===========================#

if __name__ == '__main__':  
    
    program_function = """
    *****   
    
    This script takes snps recorded in a file and sub the snp into reference sequence 
    for translation and records amino acid mutations.
    
    Ver. 0.01
           
    Extra Help:
        - For argument -c, the column count starts from 0. E.g. First column is 0, second column is 1 etc.
        
        - For argument -p, -m and -fq, the positions are in accordance to the columns selected. For example:
          if the columns in option -c are 1,5,6 and the 5th column is where nt positions are placed, then -p should
          be 1. If the 6th column is where mutations are placed, then -m should be 2 etc.
          
          *Note: This tool only takes in 3 columns of data. Any extra information will be ignored. 
        
        - If unsure if option -rp is needed, use:
            cat -A FILE
                or
            cat -v FILE
          to check whether ^M character exists, if seen, use -rp.
        
        - For argument -fup, the tool will ignore frequencies greater than inputed frequency. This is to ignore 
          outputs that fall under the case such as a wrong consensus reference is used to align and results in
          high instances of mutation recorded.

    *****
    """
    parser = argparse.ArgumentParser()
    #Required Arguments
    parser.add_argument('-f', '--file', help='File to extract snp information from.',  required = True)
    parser.add_argument('-r', '--referencefilename', help='Filename of reference file.',  required = True) 
    parser.add_argument('-d',  '--delimiter', help='Character(s) used to clearly separate columns.', required = True)
    parser.add_argument('-c',  '--columns',  help='Data desired from columns separated by user-defined delimiter.',  required = True)
    parser.add_argument('-p',  '--pos_column',  help='The column position storing the positions of mutations',  required = True,  type = int)
    parser.add_argument('-m',  '--pos_mutation',  help='The column position storing the mutations',  required = True,  type = int)
    parser.add_argument('-fq',  '--pos_freq',  help='The column position storing the frequencies',  required = True,  type = int)
    parser.add_argument('-sc', '--starting_codon',  help='starting codon position of the coding region',  required = True,  type = int)
    
    #Optional Arguments
    parser.add_argument('-e',  '--easyOutput',  help = 'Generate easyoutput for easy parsing.',  action = 'store_true')
    parser.add_argument('-o',  '--output_file',  help='Filename for the tool to write into. Otherwise default is stdout.')
    parser.add_argument('-i',  '--ignoreFirst',  help='Ignores the very first line of the input file.',  action ='store_true')
    parser.add_argument('-rp',  '--reparse',  help='Reparses input files to get rid of ^M characters from DOS files.',  action = 'store_true')    
    parser.add_argument('-lt',  '--logTen',  help = "Puts input occurrence frequencies in log base 10.",  action = 'store_true')
    parser.add_argument('-fup',  '--filter',  help="Filter frequencies above a certain value.",  type = float)
    parser.add_argument('-sf',  '--sortby_freq',  help="Sort by frequency, default sorts by nt position",  action = 'store_true')
    parser.add_argument('-is',  '--include_syn',  help='Include synonymous mutations',  action = 'store_true')
    
    if len(sys.argv) < 3:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()
    
    if(not okFile(args.file)):
        raise InputError(args.file,  'Cannot find file: ')    
    if(not okFile(args.referencefilename)):
        raise InputError(args.referencefilename,  "Invalid reference file: ")
    if(args.filter):
        if(args.filter <= 0):        
                raise InputError(args.filter,  "Filter value must be greater than 0.0. Given: ")
    elif(args.filter == None):        
        pass    
        #raise InputError(args.filter,  "Filter value must be greater than 0.0. Given: ")
    elif(args.filter <= 0):        
                raise InputError(args.filter,  "Filter value must be greater than 0.0. Given: ")    
    if(args.logTen):
        import math
        LogFrequencies = True
    if(args.columns):    
        if(re.search('-', args.columns)):            
            if(not re.search(r'^\d+-\d+$', args.columns)):        
                raise InputError(args.columns, "Invalid Columns: ")
            scope = re.search(r'^(\d+)-(\d+)$', args.columns) 
            if(int(scope.group(1)) > int(scope.group(2))):
                raise InputError(args.columns,  "Invalid Columns: ")
            columnList = args.columns.split('-')
            columnList = range(int(columnList[0]),  int(columnList[1]) + 1)
        elif(re.search(',',  args.columns)):            
            if(re.search('[^0-9,]', args.columns)):
                raise InputError(args.columns,  "Invalid Columns: ")
            columnList = args.columns.split(',')
            while('' in columnList):
                columnList.remove('')
            if(len(columnList) < 1):
                raise InputError(args.columns,  "Invalid Columns: ")            
        else:
            if(re.search('[^\d+]',  args.columns)):
                raise InputError(args.columns,  "Invalid Column: ")
            columnList = [args.columns]
    if(not args.output_file):
        WantOutput = False        
    else:
        WantOutput = True
        OutputFileName = args.output_file 
    if(args.ignoreFirst):
        IgnoreFirstLine = True
    if(args.include_syn):
        IncludeSyn = True
    if(args.reparse):
        Reparse = True
        reparsedFile = reparse(args.file)    
    
    if(Reparse):
        rawList = extractLine(reparsedFile, columnList,  args.delimiter)
    else:
        rawList = extractLine(args.file, columnList,  args.delimiter)        
   
    dataList = changeTypeThenSort(rawList, args.pos_column, args.pos_freq,  args.pos_mutation,  args.sortby_freq)    
    refFile = getFile(args.referencefilename)
    allDataList = translation(dataList,  refFile, args.starting_codon, args.pos_column, args.pos_freq,  args.pos_mutation,  args.filter)    
    printData(allDataList)    
    if(args.easyOutput):
        printEasy(allDataList)
    if(Reparse):
        os.remove(reparsedFile)
    
