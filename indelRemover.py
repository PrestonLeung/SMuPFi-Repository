#!/usr/local/bin/python
#Ver. 0.02
#Pre.L


import argparse,  re,  os 
import sys ,  itertools,  math,  gc
from Bio import SeqIO, AlignIO, Seq
from Bio.SeqRecord import SeqRecord


#========================"Global Variables"=======================#

FrequencyPosition = 0
OutputFileName = ''
WantOutput = None
IgnoreAllGaps = False

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

#==========================="getSeqLength"===========================#
#Simple function that opens a file and checks the sequence length of the fasta sequence. Assumes all
#sequence in the file are the same length (also assumes its an aligned sequence).

def getSeqLength(filename):
    try:
        OpenFastaFile = open(filename)    
    except IOError:
            print "Error: cannot find {}".format(filename)            
    else:
        record = SeqIO.parse(OpenFastaFile,  "fasta")
        length = 0
        for item in record:
            if(length != 0) :
                break            
            length = len(str(item.seq).lstrip())    
        
        OpenFastaFile.close        
        return length

#==========================="smartIgnoreCheck"===========================#
#A checker function to look +/- 2 spots from a gap character '-'. Returns a True if it's ok to ignore the gap and
#False if the it is not ok to ignore gap.

def smartIgnoreCheck(sequenceList,  position):        
    #assert (sequenceList[position] == '-')    
    seqLength = len(sequenceList)
    #print position,  seqLength

    if(position == seqLength-1):
        #print "Chk pt 1"
        #case 3 => e.g GG-
        if(sequenceList[position - 1] == sequenceList[position - 2]):
            return True
        else:
            return False
    elif(position == seqLength-2):
        #print "Chk pt 2"
        #case 1 => e.g. G-G
        if(sequenceList[position - 1] == sequenceList[position + 1]):
            return True
        #case 3 => e.g GG-
        elif(sequenceList[position - 1] == sequenceList[position - 2]):
            return True
        else:
            return False
    elif(position == 1):
        #print "Chk pt 3"
        #case 1 => e.g. G-G
        if(sequenceList[position - 1] == sequenceList[position + 1]):
            return True
        #case 2 => e.g -GG
        elif(sequenceList[position + 1] == sequenceList[position + 2]):
            return True
        else:
            return False
    elif(position == 0):
        #print "Chk pt 4"
        #case 2 => e.g -GG
        if(sequenceList[position + 1] == sequenceList[position + 2]):
            return True    
        else:
            return False
    else:
        #print 'position = ', position ,"\n"
        #case 1 => e.g. G-G
        if(sequenceList[position - 1] == sequenceList[position + 1]):
            #print "Chk pt 5.1"
            return True
        #case 2 => e.g -GG
        elif(sequenceList[position + 1] == sequenceList[position + 2]):
            #print "Chk pt 5.2"
            return True
            #case 3 => e.g GG-
        elif(sequenceList[position - 1] == sequenceList[position - 2]):
            #print "Chk pt 5.3"
            return True
        else:
            #print "Chk pt 5.4"
            return False
        

#==========================="indelRemover"===========================#
#Removes '-' characters and replaces with the reference character.

def indelRemove(ref,  fas, start,  end):
        
    refSeq = ''
    allSeq = {}
    
    for item in ref:        
       refSeq = item.seq     
    for item in fas:        
        fasSeq = list(item.seq)     
        
        gapPosList = [x for x, y in enumerate(fasSeq) if y == '-']        
        #print item.id, "\t",  gapPosList,  "\n\n"
        
        if(not IgnoreAllGaps):
            (single, double, triple) = distributeGaps(sorted(gapPosList))
            
            fasSeq = manageSingle(single,  fasSeq,  refSeq,  start-1)
            fasSeq = manageDouble(double, fasSeq, refSeq, start-1)
            fasSeq = manageTriple(triple, fasSeq,  refSeq, start-1)
            
        else:        
           for position in gapPosList:
                fasSeq[position] = refSeq[position + (start-1)]
                
        if(FrequencyPosition):            
            gotFreq = re.search('(\d+\.\d+e?-?\d+)', item.description.split('_')[FrequencyPosition-1])
            if(not gotFreq):
                raise InputError(FrequencyPosition,  "Incorrect frequency position: ")                
            freq = float(gotFreq.group(1))
        else:    
            freq = 1.0
        #print ''.join(fasSeq)
        keySeq = re.sub("[\[\],' ]", '',str(fasSeq))
        
        if (keySeq in allSeq):
            allSeq[keySeq] += freq
        else:
            allSeq[keySeq] = freq
       
    printAllTheStuff(allSeq)      
 

#==========================="manageSingle"===========================#
 #Manages single deletions and substitutes them when it's adjacent to homopolymer region
 #Note variable 'start' is the array starting position of reference sequence. 
def manageSingle(singleList, seq, refseq,  start):
    seqLength = len(seq)
    okToReplace = False
    
    for singlet in singleList:        
        pos = singlet[0]        
        
        if(smartIgnoreCheck(seq, pos) or smartIgnoreCheck(refseq, (pos + start))):
            seq[pos ] = refseq[pos + start].upper()
            
    #print ''.join(seq)
    return seq

#==========================="manageDouble"===========================#
 #Manages double deletions and substitutes them when it's adjacent to homopolymer region
 #Note variable 'start' is the array starting position of reference sequence. 
def manageDouble(doubleList, seq, refseq,  start):
    seqLength = len(seq)
    
    
    for duplet in doubleList:
        okToReplace = False
        firstPos = duplet[0]
        lastPos = duplet[1]
        #check the left hand side first
        if(not firstPos < 2):
            #E.g GG--
            if(seq[firstPos - 1] == seq[firstPos - 2]):
                okToReplace = True        
        #then the right hand side of deletion
        if(not lastPos > (len(seq) - 3)):
            #E.g. --GG
            if(seq[lastPos + 1] == seq[lastPos + 2]):
                okToReplace = True            
        #check if deletion lies in a homopolymer region on the reference as well
        if(smartIgnoreCheck(refseq, (firstPos + start)) or smartIgnoreCheck(refseq, (lastPos + start))):
                okToReplace = True            
        if(okToReplace):            
            for i in range(firstPos - 1,  lastPos + 2):
                if(i == lastPos + 1 and seq[i] == '-'):
                    continue
                    #or break?
                seq[i] = refseq[i + start].upper()
    return seq
    


#==========================="manageTriple"===========================#
#Manages triple deletions and substitutes them when it's adjacent to homopolymer region
#Note variable 'start' is the array starting position of reference sequence. 
def manageTriple(tripleList, seq, refseq,  start):
    seqLength = len(seq)    
    
    for triplet in tripleList:
        okToReplace = False
        firstPos = triplet[0]
        lastPos = triplet[2]
        #check the left hand side first
        if(not firstPos < 2):
            #E.g GG---
            if(seq[firstPos - 1] == seq[firstPos - 2]):
                okToReplace = True
        
        #then the right hand side of deletion
        if(not lastPos > (len(seq) - 3)):
            #E.g. ---GG
            if(seq[lastPos + 1] == seq[lastPos + 2]):
                okToReplace = True
    
        #check if deletion lies in a homopolymer region on the reference as well
        if(smartIgnoreCheck(refseq, (firstPos + start)) or smartIgnoreCheck(refseq, (lastPos + start))):
                okToReplace = True  
        
        
        if(okToReplace):            
            for i in range(firstPos - 1,  lastPos + 2):
                if(i == lastPos + 1 and seq[i] == '-'):
                    continue
                    #or break?
                seq[i] = refseq[i + start].upper()
    return seq
 
 
 
#==========================="distributeGaps"===========================#
#Distributes gaps into '-' gap, '--' gaps and '---' gaps.
#Note to self, if gaps occur in '----', then it will split to '---' and '-'. If it's 5, then it would split 3:2 and so on.

def distributeGaps(gapList):
    comboCounter = 1
    currentGapPos = None
    comboList = []
    singleList, doubleList, tripleList = ([], [], [])
    
    for item in gapList:        
        if(currentGapPos == None):
            currentGapPos = item
            comboList.append(item)
            continue
        if(item == currentGapPos + 1 and comboCounter < 3):
            comboCounter += 1            
        else:
            if(comboCounter == 1):                
                singleList.append(comboList)
            elif(comboCounter == 2):               
                doubleList.append(comboList)
            elif(comboCounter == 3):                
                tripleList.append(comboList)
            comboList = []
            comboCounter = 1            
        comboList.append(item)
        currentGapPos = item
    
    if(len(comboList) == 1):
        singleList.append(comboList)
    elif(len(comboList) == 2):        
        doubleList.append(comboList)
    elif(len(comboList) == 3):        
        tripleList.append(comboList)
    
    return singleList, doubleList, tripleList
    

#==========================="printAllTheStuff"===========================#
#Prints all the stuff inside allSeq (name says it all >:]).

def printAllTheStuff(allSeq):

    #Reminder, sorts allSeq into tuples
    allSeq = sorted(allSeq.items(), key=lambda x: x[1], reverse=True)
    hapNum = 0
    
    if(WantOutput):
        openFile = open(OutputFileName, 'w')
    else:
        openFile = sys.stdout
    
    for seq in allSeq:        
        seq_record = SeqRecord(Seq.Seq(seq[0]))        
        seq_record.id = 'HAP' + str(hapNum) + '_' + str(seq[1])
        seq_record.description = ''
        SeqIO.write(seq_record, openFile, 'fasta')
        hapNum += 1
    
    if(WantOutput):
        openFile.close()
    
#==========================="processInput"===========================#
#Handles the regions where the fasta file lie within the reference file. Function takes the smaller value
#between the length of the user-specified region and the length of the fasta sequence length as the 
#ending position.
#Reference sequence is assumed to be always longer or equal the length of fasta sequence, if it is not case,
#it will throw an error.

def processInput(refFile,  fasFile,  start,  end):
    
    refLength = getSeqLength(refFile)
    fasLength = getSeqLength(fasFile)
    regionLength = (end - start) + 1
    ref = getFile(refFile)
    fas = getFile(fasFile)
    minLength = min(regionLength,  fasLength)
    
    if(not minLength <= refLength):
        raise InputError(minLength, "Reference sequence length is smaller than query length: ")
    
    indelRemove(ref,  fas,  start,  (start + minLength - 1))
    
#==========================="main function"===========================#

if __name__ == '__main__':   
    
    program_function = """
    *****
    IndelRemover:
        
        A tool to remove indels from fasta sequences and replace it with reference sequence
        characters and merges fasta sequences that end up being identical.
        
        Ver: 0.02
            - Refined homopolymer region indel removal algorithm.
                -> Double and Triple indels are now considered.
                -> Anything above Triple indel are split into sets 
                    of 1,2 or 3. E.g. a '----' indel will be split into 
                    '---' and '-'.
                -> Uses reference file to check if an indel lies
                    on homopolymer region in addition to the
                    old algorithm that only looks at fasta sequence file.
        
        Extra Help:
        
            The frequency value in the Fasta name (e.g. >FastaName_0.123)
            is split via underscores '_'. In the example given in the brackets, 
            the frequency index will be in the 2nd position. If your Fasta name 
            looks like ">FastaName0.123", then the frequency index will be in 
            1st position.
            
            If there is no frequency in the Fasta name to use 
            (leave -fr option blank), all sequences will be counted as 1.        

        Requires:
            - Python2.7
            - Biopython 

    *****
    """
    parser = argparse.ArgumentParser()
    #Required Arguments
    parser.add_argument('-f', '--fastafilename', help='Filename and aligned FASTA file.',  required = True)
    parser.add_argument('-r', '--referencefilename', help='Filename of reference file.',  required = True) 
    parser.add_argument('-s',  '--scope',  
                                      help ='Scope area of the reference genome. E.g. Nucleotide region 650-900 of reference sequence.', 
                                      required = True)
    #Optional Arguments
    parser.add_argument('-fr',  '--freq_index',  help='position of frequency in the tag, seperated by "_".',  type = int)    
    parser.add_argument('-o', '--output-file', help='output filename of result.')
    parser.add_argument('-a', '--replace_all', help='any indel "-" characters are replaced with reference when found', action='store_true')

    if len(sys.argv) < 2:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()
    
    cleanScope = re.sub('\s+',  '',  args.scope)
    refScopeTest = re.search(r'^\d+-\d+$', cleanScope)
    if(not refScopeTest):        
        raise InputError(args.scope, "Invalid Region: ")
    refScope = re.search(r'^(\d+)-(\d+)$', cleanScope) 
    if(int(refScope.group(1)) > int(refScope.group(2))):
        raise InputError(args.scope,  "Invalid Region: ")
    if(int(refScope.group(1)) == 0):
        raise InputError(args.scope,  "Invalid start: ")
    if(not okFile(args.fastafilename)):
        raise InputError(args.fastafilename,  "Invalid fasta file: ")
    if(not okFile(args.referencefilename)):
        raise InputError(args.referencefilename,  "Invalid reference file: ")
    if(args.freq_index):
        FrequencyPosition = int(args.freq_index)
    if(not args.output_file):
        WantOutput = False        
    else:
        WantOutput = True
        OutputFileName = args.output_file        
    if(args.replace_all):
        IgnoreAllGaps = True
    processInput(args.referencefilename,  args.fastafilename,  int(refScope.group(1)),  int(refScope.group(2)))
    
