#!/usr/local/bin/python
#Version 0.869AA -

# This script version deals with ShoRAH6. 

import argparse,  re,  os 
import sys ,  itertools,  math,  gc
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalwCommandline
from sets import Set

#========================"Global DataStructures and Variables"=======================#
#Default posterior threshold. 
Threshold = 0.85
#Default global frequency cutoff.
GThreshold = 0.005
#Default mutation combination range.
MCombinationMax = 0
MCombinationMin = 1
#Save directory path. All outputs will be saved to this path.
SaveDirectory = ''
#Easy ouput. Yes or No. Defaults at N. Ony Y if user specifies in argument.
EasyOutput = 'N'
#Global flag. Defaults at 'N' unless user specifies input files to be global sequences
GlobalRun = 'N'
#Variable to notify SMuPFi to ignore gaps if it's 'Y'. Otherwise defaults at 'N'
IgnoreGap = 'N'
#Sum of all Avg Reads
TotalAvgRead= 0.0
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


#========================"initialSearch"========================#
#This loop iterates through a bioseq record list. This loop will output the name of the query 
#sequence and the reference sequence with markers to identify where and what the difference is. 
#Fills up global structure diffHash, re-defines global MCombinationMax if needed and defines
#global frequency sum (TotalAvgRead).
#
#Input:
#   - Reference sequence list (seqRef) containing the record of the reference sequence, with
#     relevant boundaries (refLowBound, refHighBound) linking the reference to the query sequences.
#   - Query sequence list (seqRec) containing multiple records of sequences the user wants to compare.
#     The boundaries (lowerBound, upperBound) link the query sequences to the reference sequence.
#   - The largest length possible (maxScope) indicated by the user to search through.
#Output:
#   - Stdout will print the difference between reference and query.
#   - Returns a dictionary filled with found differences and the filename used to store results.

def initialSearch(seqRef, seqRec, lowerBound, upperBound, refLowBound,  refHighBound, maxScope,  startingCodon):
    global TotalAvgRead
    global MCombinationMax
    diffHash = {}
    region = [str(refLowBound),  str(refHighBound)]
    fileName = '-'.join(region)
    startingAA = frameShift(refLowBound,  0,  startingCodon)    
    unalignedAA = {}
    AAFreqTable = {}
    AAHapTable = {}
    AASimplified = {}
    beforeMut = {}
    dynamicSpace = 0
    
    try:        
        if(SaveDirectory):
            output = open(SaveDirectory + '/'+ fileName + '.txt',  'w')
        else:    
            output = open(fileName+'.txt',  'w')
        try:
            output.write("Looking at a maximum of {} nucleotides. ".format(maxScope))
            if(GlobalRun == "Y"):
                output.write("Comparing position(s) {} to {} with threshold {}.\n".format(lowerBound,  upperBound,  GThreshold))   
            else:
                output.write("Comparing position(s) {} to {} with threshold {}.\n".format(lowerBound,  upperBound,  Threshold))   
            output.write("at nucleotide reference regions {} to {}.\n".format(refLowBound,  refHighBound))            
            output.write("Amino acid region {} to {}.\n\n".format(startingAA,  frameShift(refHighBound,  0, startingCodon)))
            output.write("Reference:\n")
            
            for seqReference in seqRef:
                pre_AARef = seqReference.seq[startingCodon-1 : len(seqReference)].translate()                
                my_AARef = pre_AARef[(startingAA - 1) : frameShift(refHighBound,  0, startingCodon)]

                output.write("%s" %my_AARef + '\n')                
                
                for singleRecord in seqRec:                    
                    my_posString = ''
                    aaSeqHash = {}
                    refMutHash = {}
                    aaMutPos = set()                    
                    
                    my_AATranslation = codonProcess(str(singleRecord.seq),  lowerBound, startingCodon)                    
                    
                    if(GlobalRun == 'Y'):
                        stringMatch = re.search(r'(.*)_(.*)',  singleRecord.description,  re.M)                        
                    else:
                        stringMatch = re.search(r'.*posterior=(.*) ave_reads=(.*)',  singleRecord.description,  re.M)
                    
                    if (GlobalRun == 'Y'):                        
                        if(float(stringMatch.group(2)) <= GThreshold):                            
                            continue
                    else:
                        if(float(stringMatch.group(1)) <= Threshold):
                            continue                             
                    
                    #This is conservative method
                    #assert (len(my_AATranslation) == len(my_AARef))
                    if(len(my_AATranslation) > dynamicSpace):                        
                        dynamicSpace = len(my_AATranslation)
                    
                    
                    for i in range (0,  len(my_AARef)):                        
                        if(my_AATranslation[i] == '?'):
                            my_AATranslation = my_AATranslation.replace('?', my_AARef[i],  1)
                        if(my_AATranslation[i] != my_AARef[i]):
                            my_posString += my_AATranslation[i]
                            aaSeqHash[i + startingAA] =  my_AATranslation[i]
                            refMutHash[i + startingAA] = my_AARef[i]
                        else:                            
                            my_posString += '.'
                    
                    if(AAFreqTable.has_key(my_AATranslation)):                            
                        AAFreqTable[my_AATranslation] += float(stringMatch.group(2))
                        if(GlobalRun == 'Y'):
                            AAHapTable[my_AATranslation] += (', '+ re.sub('_\d+\.?\d*.*$', '', singleRecord.id))                        
                        else:
                            AAHapTable[my_AATranslation] += (', '+ re.sub('\|posterior.*', '', singleRecord.id))                        
                    else:    
                        AAFreqTable[my_AATranslation] = float(stringMatch.group(2))                        
                        
                        if(GlobalRun == 'Y'):
                            AAHapTable[my_AATranslation] = re.sub('_\d+\.?\d*.*$', '', singleRecord.id)
                        else:
                            AAHapTable[my_AATranslation] = re.sub('\|posterior.*', '', singleRecord.id)
                        AASimplified[my_AATranslation] = my_posString

                    diffHash[singleRecord.description] = aaSeqHash
                    beforeMut[singleRecord.description] = refMutHash
                    
                    if(not args.mutationcombination):                        
                        if(len(diffHash[singleRecord.description]) > Max):
                            MCombinationMax = len(diffHash[singleRecord.description])
                    TotalAvgRead= TotalAvgRead + float(stringMatch.group(2))
                    
                
                for keys in sorted(AAFreqTable,  key = len,  reverse = True):                         
                        output.write(str(AASimplified[keys]).ljust(dynamicSpace + 10))
                        output.write(str(math.ceil(AAFreqTable[keys] / TotalAvgRead* 1000.0) / 1000.0).ljust(10))                        
                        output.write(str(AAHapTable[keys]) + '\n')


        finally:
            output.write('\n\n')            
            output.close()
            return fileName,  diffHash,  beforeMut
    except IOError:
        print "Cannot write file {}.txt".format(fileName)


#========================"secondSearch"========================#
#This loop iterates through diffHash and searches for common mutations found through initialSearch.
#
#Input:
#   - a start and end point for file naming
#Output:
#   - Stdout a summary table of the shared patterns found and patterns observed in the pool of reads. .

def secondSearch(diffHash,  fileName, refHash):         
    allNames,  allAvgSum,  allRef  = ({}  ,  {},  {})    
    uniqueFrequency, uniquePattern= ({},  {})    
    uniqueName,  uniqueRefPat = ({},  {})     
    relationCheckSet = set()
    global TotalAvgRead
    totalHapCount = len(diffHash)
   
    for haploKey1 in sorted(diffHash.iterkeys()):
        
        #a hash table storing the positions (key) of the combination where mutations can occur (value, n times).
        #cMutation and commonNT uses the same set of keys           
        cMutationSet = {}
        avgSum = {}
        #a hashTable linking the pattern of mutation(value) to their position(key). 
        #E.g. Pos 678,900,1000 -> C,A,T
        commonNT = {}
        #a hash table linking names of haplotypes (the value is a set of haplo ids) 
        #carrying the same mutation pattern (key)
        patternAndHap = {}
        shRefPat = {}
        
        if(GlobalRun == 'Y'):            
            hapName = re.search('(.*)_.*', haploKey1)
            hapAvg1 = re.search('.*_(.*)', haploKey1)            
        else:        
            hapName = re.search(r'^(hap_\d+)', haploKey1)
            hapAvg1 = re.search(r'(\d+\.?\d*)$', haploKey1)
        
        uniquePosKey = tuple(sorted(diffHash[haploKey1].iterkeys()))
        #print type(uniqueKey)
        testAmino = getLetters(haploKey1,  uniquePosKey,  diffHash)
        uniqueKey = (tuple(sorted(diffHash[haploKey1].iterkeys())), tuple(testAmino))        
        if(len(uniqueKey) > 0):            
            if(not uniqueFrequency.has_key(uniqueKey)):
                uniqueFrequency[uniqueKey] = float(hapAvg1.group(1))
                uniquePattern[uniqueKey] = getLetters(haploKey1, uniquePosKey, diffHash)
                uniqueName[uniqueKey] = [hapName.group(1)]
                uniqueRefPat[uniqueKey] = getLetters(haploKey1,  uniquePosKey,  refHash)
            else:    
                uniqueFrequency[uniqueKey] += float(hapAvg1.group(1))
                uniqueName[uniqueKey].append(hapName.group(1))
        
        for haploKey2 in sorted(diffHash.iterkeys()):            
            #making sure not to compare values with the same key since they will be identical
            if (cmp(haploKey1,  haploKey2) == 0):
                continue
            sortedRelation = tuple(sorted((haploKey1, haploKey2)))
            if(sortedRelation in relationCheckSet):                
                continue
            relationCheckSet.add(sortedRelation)
            
            #the common mutation set between two haplotypes using haploKey1 and haploKey2
            #this returns the common keys (which contain the position of mutations) of the hash stored
            #inside diffHash. 
            commonSet = set(diffHash[haploKey1]) & set(diffHash[haploKey2]) 
            for cLength in xrange (MCombinationMax, (MCombinationMin - 1),  -1):
                #checking the possible mutations from the commonSet
                for subset in itertools.combinations(sorted(commonSet),  cLength):
                    if(len(subset) > 0):                        
                        okCombination = 1
                        for position in subset:
                            #checking if the mutations in the same positions are common mutations
                            #if not then it is not an okCombination
                            if(not diffHash[haploKey1][position] == diffHash[haploKey2][position]):
                                okCombination = 0
                                break
                        
                        if(okCombination == 1): 
                            if(GlobalRun == 'Y'):
                                hapAvg2 = re.search('.*_(.*)', haploKey2)
                                hapName2 = re.search('(.*)_.*', haploKey2)
                            else:    
                                hapAvg2 = re.search(r'(\d+\.?\d*)$', haploKey2)
                                hapName2 = re.search(r'^(hap_\d+)', haploKey2)
                            
                            if(cMutationSet.has_key(subset)):
                                cMutationSet[subset] += 1
                                avgSum[subset] += float(hapAvg2.group(1))       
                                patternAndHap[subset].append(hapName2.group(1))
                            else:
                                cMutationSet[subset] = 2                                
                                avgSum[subset] = float(hapAvg1.group(1)) + float(hapAvg2.group(1))
                                patternAndHap[subset] = []
                                patternAndHap[subset].append(hapName.group(1))
                                patternAndHap[subset].append(hapName2.group(1))
                            #retrieving the mutation nucleotides using the positions that have been filtered
                            #print type(subset)
                            shRefPat[subset] = getLetters(haploKey2,  subset,  refHash)                            
                            commonNT[subset] = getLetters(haploKey1,  subset,  diffHash)
        keyList = cMutationSet.iterkeys()
        sortedKeyList = sorted(keyList,  key = len,  reverse = True)        
        #Saving the common mutations we found
        #Only saves the common mutations that cannot be found in the hashtables
        #E.g. if it is there, then it will not write anything
        for key in sortedKeyList:            
            if(cMutationSet[key] > 1): 
                bigKey = (key, tuple(commonNT[key]))
                if(bigKey not in allNames):    
                    allNames[bigKey] = patternAndHap[key] 
                if(bigKey not in allAvgSum):
                    allAvgSum[bigKey] = avgSum[key]
                if(bigKey not in allRef):
                    allRef[bigKey] = shRefPat[key]
        cMutationSet.clear()
        avgSum.clear()
        patternAndHap.clear()
        commonNT.clear()
        shRefPat.clear()
    relationCheckSet.clear()
    #Printing the common mutations we've just found
    sortedKeyList2 = sorted(allNames.iterkeys())    
    sortedKeyList3 = sorted(uniqueFrequency.iterkeys())

    
    if(len(sortedKeyList2) > 0 or len(sortedKeyList3) > 0):
        try:
            if(SaveDirectory):
                output = open(SaveDirectory + '/'+ fileName + '_Summary.txt',  'w')
            else:    
                output = open(fileName+'_Summary.txt',  'w')
            try:                
                output.write("Summary of common mutations with length from "+str(MCombinationMin)+ '-'+ str(MCombinationMax)+ ":\n\n")
                output.write("Position(s)\tPattern\tShared Occurrence\n\n")

            
                
                if(len(sortedKeyList2) > 0):
                    for key in sortedKeyList2:                    
                                            
                        if(len(key[0]) < 2):
                            string = re.sub('[,()]',  '', str(key[0]))                        
                            output.write(string + "\t")
                        else:
                            string = re.sub('[,()]',  '', str(key[0]))                        
                            output.write(string + "\t")
                        
                                          
                        if(len(key[1]) < 2):
                            string = re.sub('[(,)\']',  '', str(key[1]))                        
                            output.write(string + "\t")
                        else:
                            string = re.sub('[(,)\' ]',  '', str(key[1]))                        
                            output.write(string + "\t")
                        
                        calcAvgSum = float(allAvgSum[key]) / TotalAvgRead
                        calcAvgSum = math.ceil(calcAvgSum * 1000.0) / 1000.0                    
                        output.write(str(calcAvgSum) + "\t")
                        
                        string = re.sub('[\],()\[\']',  '', str(allNames[key]))  
                        output.write(string + "\t\n")
                                     
            finally:
                output.close()                
                print "File {}.txt done.".format(fileName)
                print "File {}_Summary.txt done.".format(fileName)
        except IOError:
                print "Cannot write file {}.txt".format(fileName)
    else:
        try:
            if(SaveDirectory):
                output = open(SaveDirectory + '/'+ fileName + '.txt',  'a')
            else:
                output = open(fileName + '.txt',  'a') 
            try:
                output.write("Common differences not found 378.")
            finally:
                output.close()
                print "File {}.txt done".format(fileName)
        except IOError:
            print "Cannot write file {}.txt".format(fileName)
    
    if(EasyOutput == 'Y'):
        printEasyOutputShared(allNames, allAvgSum,  sortedKeyList2, fileName,  totalHapCount, allRef)
        printEasyOutputUnique(sortedKeyList3,  uniqueName,  uniqueFrequency,  uniquePattern,  fileName,  totalHapCount,  uniqueRefPat)   
    allNames.clear()  
    allAvgSum.clear()
    uniqueFrequency.clear()
    uniquePattern.clear()
    uniqueName.clear()
    del sortedKeyList2[:]
    del sortedKeyList3[:]
    allNames = allAvgSum = uniqueName = uniqueFrequency = uniquePattern = None
    sortedKeyList2 = sortedKeyList3 = None


#========================"constructGap"========================#
#Helper function to create the right spacing for printing using sortedKeyList2.

def constructGap(aList):
    
    currentGap = 0
    for item in aList:
        if(len(str(item[0][-1])) > currentGap):
            currentGap = len(str(item[0][0]))
    return currentGap        


#========================"getLetters"========================#
#Helper functions to turn retrieve patterns from given positions.
#Input:
#   - Key (dictKey) to unlock the haplotype
#   -  a list/set/tuple of positions (positionSet)
#   - a dictionary (refHash) for the function to retrieve pattern 
#     using the key and set of positions.
#Output:
#   - returns a list containing the pattern (stringList)

def getLetters(dictKey,  positionSet, refHash):
    positionList = list(positionSet)
    stringList = []
    
    for item in positionList:
        stringList.append(refHash[dictKey][item])
    
    return stringList

#========================"codonProcess"========================#
#Includes deletions and translates its as a '?' if it see's any codons with '-'

def codonProcess(seq,  position,  startCodon):    
    aminoAcid = '' 
    codePos = position - (startCodon - 1)
    r = codePos % 3
    
    if(r == 0):
        seq = '--' + seq        
    if(r == 2):
        seq = '-'+ seq 
    
    codonList = re.findall('.?.?.?',  seq)    
    for triplet in codonList:        
        if(len(triplet) < 3 and len(triplet) > 0):            
            aminoAcid+= '?'
        else:            
            if('-' in triplet):
                aminoAcid+= '?'
            else:
                aminoAcid += str(Seq(triplet).translate()) 
    return aminoAcid        

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

#========================"getFileLength" ==========================#
#This function opens up a file and retrieves the amount of lines in that file

def getFileLength(path):
    with open(path) as file:
        for i,  l in enumerate(file):
            pass
    return i + 1

#========================"printEasyOutputUnique" ==========================#
#This function prints out a simplified summary table ONLY if the user specifies it.

def printEasyOutputUnique(keyList,  uniqueName,  uniqueFrequency,  uniquePattern,  fileName,  totalHapCount,  uniqueRef):    
    try:
        if(SaveDirectory):                            
            eOutput = open(SaveDirectory + '/'+ fileName + '_EasyOutputUnique.txt',  'w')
        else:                   
            eOutput = open(fileName + '_EasyOutputUnique.txt',  'w')
        try:
            for key in keyList:
                keyString0 = re.sub('[(),\']',  '', str(key[0]))
                keyString1 = re.sub('[\[(),\'\]]',  '', str(uniquePattern[key]))
                keyString2 = re.sub('[\[(),\'\]]',  '', str(uniqueRef[key]))
                eOutput.write(keyString0 + "|")
                eOutput.write(keyString1 + "|")                    
                calcAvgSum = float(uniqueFrequency[key]) / TotalAvgRead
                calcAvgSum = math.ceil(calcAvgSum * 1000.0) / 1000.0
                eOutput.write(str(calcAvgSum) + "|")                        
                eOutput.write(re.sub('[\[\',\]]', '', str(uniqueName[key])) + "|")
                eOutput.write(str(totalHapCount) + "|") 
                eOutput.write(keyString2 + "\n")
        finally:
            eOutput.close()
            print "File "+ fileName + "_EasyOutputUnique.txt done."
    except IOError:
        print "Cannot write easy output {}_EasyOutputUnique.txt".format(fileName)

    except IOError:
        print "Cannot write easy output {}_EasyOutputUnique.txt".format(fileName)     

#========================"printEasyOutputShared" ==========================#
#This function prints out a simplified summary table ONLY if the user specifies it.

def printEasyOutputShared(nameshash,  counthash,  keyList,  fileName,  totalHapCount,  refpat):
    if(len(keyList) > 0): 
        try:
            if(SaveDirectory):                            
                eOutput = open(SaveDirectory + '/'+ fileName + '_EasyOutputShared.txt',  'w')
            else:                   
                eOutput = open(fileName + '_EasyOutputShared.txt',  'w')
            try:
                for key in keyList:    
                    keyString0 = re.sub('[(),\']',  '', str(key[0]))
                    keyString1 = re.sub('[(),\']',  '', str(key[1]))
                    keyString2 = re.sub('[\[(),\'\]]',  '', str(refpat[key]))
                    eOutput.write(keyString0 + "|")
                    eOutput.write(keyString1 + "|")                    
                    calcAvgSum = float(counthash[key]) / TotalAvgRead
                    calcAvgSum = math.ceil(calcAvgSum * 1000.0) / 1000.0
                    eOutput.write(str(calcAvgSum) + "|")                        
                    eOutput.write(re.sub('[\[\',\]]', '', str(nameshash[key])) + "|")
                    eOutput.write(str(totalHapCount) + "|")
                    eOutput.write(keyString2 + '\n')
                    
            finally:
                eOutput.close()
                print "File "+ fileName + "_EasyOutputShared.txt done."
        except IOError:
            print "Cannot write easy output {}_EasyOutputShared.txt".format(fileName)
    else:
        try:
            if(SaveDirectory):
                eOutput = open(SaveDirectory + '/'+ fileName + '_EasyOutputShared.txt',  'w')
            else:
                eOutput = open(fileName + '_EasyOutputShared.txt',  'w')    
            try:                    
                eOutput.write("Common differences not found 646.")
            finally:
                eOutput.close()
                print "File "+ fileName + "_EasyOutputShared.txt done."
        except IOError:
            print "Cannot write easy output {}_EasyOutputShared.txt".format(fileName) 


#========================"checkRegion"========================#
#Checks if the focus region entered by user is legit. If it is NOT
#will edit the range and does the minimum of either the fasta sequence length
#or user specified length. E.g. Focus region is 1-400 but the fasta length is only 300. It will change to 1-300.
#Conversely, if the user specifies a shorter one E.g. 1-100 while fasta length is 300. It will change to 1-100.

def checkRegion(refScopeLow, refScopeHigh, recLow, recHigh):        
    refScopeLength = (refScopeHigh - refScopeLow) + 1
    recScopeLength = (recHigh - recLow) + 1
    
    if (refScopeHigh < refScopeLow):
        raise InputError(args.scope, "Invalid reference scope: ")
    #if the given reference scope region is wider than the query scope region
    #then the limits will be automatically switched to maximum of query scope region
    if (not refScopeLength <= recScopeLength):
       refScopeLength = recScopeLength    
   
    return (recLow,  (recLow + refScopeLength - 1), refScopeLow, (refScopeLow + refScopeLength - 1),  refScopeLength)

#=========================="tryRecFile"================================#
#A subfunction to try and open the given file. If the file cannot be opened, throws an error.

def tryRecFile(filename):
    try:
        OpenFastaFile = open(filename)    
    except IOError:
            print "Error: cannot find {}".format(filename)            
    else:
        record = SeqIO.parse(OpenFastaFile,  "fasta")
        OpenFastaFile.close        
        return record
        
#=========================="getGlobalLength"================================#
#A subfunction to try and open the given file. If the file cannot be opened, throws an error. 
#Function differs to the above function because it returns the length of the record file.
#This because there is no other way to retrieve the length besides opening it up and manually
#check the length. Unlike local windowed sequence files of ShoRAH, the file names do not contain the regions.

def getGlobalLength(filename):
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


#==========================="processDir"========================#
#processDir does the distribution of work after the main function has taken in all the user inputs.
#This sub function will process all the data given the conditions provided  by the user.
#It calls initialSearch() and secondSearch() functions after it has checked that user inputs are OK.

def processDir (pathToFasta, pathToReference, scope,  startingCodon):      
    #Check the scope input that the form is correct and throws error if user tries to put 89d0 - )*^8 as a region    
    cleanScope = re.sub('\s+',  '',  scope)    
    refScopeTest = re.search(r'^\d+-\d+$', cleanScope)
    global TotalAvgRead
    fileName = ''
    startingAA = 0
    cHash = {}
    ncHash = {}
    ref = {}
    
    if(not refScopeTest):        
        raise InputError(scope, "Invalid Region: ")
    refScope = re.search(r'^(\d+)-(\d+)$', cleanScope) 
    if(int(refScope.group(1)) > int(refScope.group(2))):
        raise InputError(scope,  "Invalid Region: ") 
    if(not os.path.isfile(pathToReference)):
        raise InputError(pathToReference, "Invalid reference file: ")    
    if(os.path.isdir(pathToFasta)):
        fileList = os.listdir(pathToFasta)        
        for fastaFile in sorted(fileList):            
            TotalAvgRead= 0.0
            if(not re.search('reads-support.fas$',  fastaFile)):
                continue                                        
            
            searchLimit = re.search('(\d+)-(\d+)\.reads-support\.fas$',  fastaFile) 
            
            if((int(searchLimit.group(1)) -1) > int(refScope.group(2))):    
                raise InputError(searchLimit.group(1),  "Out of Bound at starting position: ")            
            
            recordLength = int(searchLimit.group(2)) - int(searchLimit.group(1)) + 1
            start = int(searchLimit.group(1))
            end = min(int(refScope.group(2)), int(searchLimit.group(2)))
            
            (recLow,  recHigh, 
             refLow,  refHigh,  maxLength) = checkRegion(start, end, int(searchLimit.group(1)),  int(searchLimit.group(2)))
            
            seqRec = tryRecFile(os.path.join(pathToFasta, fastaFile))
            seqRef = tryRecFile(pathToReference)
           
            fileName, cHash,  ref = initialSearch(seqRef,  seqRec,  int(recLow),  int(recHigh),  int(refLow),  int(refHigh),  int(maxLength),  startingCodon)
            
            secondSearch(cHash,  fileName,  ref)
                
    elif(os.path.isfile(pathToFasta)):
        seqRec = tryRecFile(pathToFasta)

        seqRef = tryRecFile(pathToReference)
        
        pathFasta,  fastaFile =  os.path.split(args.fastafilename)
        
        if(GlobalRun == 'Y'):
#                if(not re.search('.popl$',  fastaFile)):
#                    raise InputError(args.fastafilename,  "File is not in correct .popl format. Given argument: ")    
#                    sys.exit()                
                recLength = getGlobalLength(pathToFasta)                    
                refLength = int(refScope.group(2)) - int(refScope.group(1)) + 1
                recLow, refLow = int(refScope.group(1)),  int(refScope.group(1))
                refHigh = min(recLength,  refLength) + refLow - 1
                recHigh = min(recLength,  refLength) + recLow - 1
                maxLength = min(recLength,  refLength)                

                fileName, cHash,  ref = initialSearch(seqRef,  seqRec,  recLow, recHigh,  refLow,  refHigh,  maxLength,  startingCodon)
                
                secondSearch(cHash,  fileName,  ref)                
        else:    
            searchLimit = re.search('(\d+)-(\d+)\.reads-support\.fas$',  fastaFile)               
        
            (recLow,  recHigh, 
             refLow,  refHigh,  maxLength) = checkRegion(int(refScope.group(1)), int(refScope.group(2)),  
                                                                       int(searchLimit.group(1)),  int(searchLimit.group(2))) 
        
            fileName,  cHash,  ref = initialSearch(seqRef,  seqRec,  int(recLow),  int(recHigh),  int(refLow),  int(refHigh),  int(maxLength),  startingCodon)
            
            secondSearch(cHash,  fileName, ref) 
                        
    else:
        raise InputError(args.fastafilename,  "Is not a directory or file. Given argument: ")
        sys.exit()
    

#==========================="main function"===========================#
#Main function does the parsing of the arguments given by the user
#It will then open the files and store it as a BioSeq record using biopython and then call the 
#sub functions to process.

if __name__ == '__main__':   
    
    program_function = """
    *****
    ~a-SMuPFi (amino acid - Shared Mutation Pattern Finder)~

 	A tool that takes in sequences and searches for shared mutations!	 

	Version 0.869AA

    *****
    """
    parser = argparse.ArgumentParser()
    #required arguments
    parser.add_argument('-f', '--fastafilename', help='Filename or directory of sequence(s) in FASTA format.',  required = True)
    parser.add_argument('-r', '--referencefilename', help='Filename of reference sequence.',  required = True)    
    parser.add_argument('-s',  '--scope',  
                                      help ='Scope area of the reference genome. E.g. Nucleotide region 650-900 of reference sequence.', 
                                      required = True)
    parser.add_argument('-sc',  '--starting_codon',  
                                        help='''The first position of the nucleotide of the start codon in the coding region in query. E.g. Need the
                                                    position of A in ATG (starting codon).''',  required = True,  type = int)
    #optional arguments
    parser.add_argument('-m', '--mutationcombination', 
                                       help = '''Number of mutation combination desired. If not specified, will take the largest 
                                                    mutation length found as the value, starting from 1. Takes in a range as input. E.g. 1-2 means looking for
                                                    min of 1 up to max of 2 mutation combination. If only looking for one specific combination, enter,  
                                                    for example, 5-5 to get min of 5 and max of 5.''')
    parser.add_argument('-t', '--threshold', help='Minimum posterior requirement. Default is 0.85.',  type=float)
    parser.add_argument('-gf', '--globalfrequency',  
                                        help="Threshold value for frequency of occurrence of haplotypes in global files.",  type = float)
    parser.add_argument('-d', '--directory',  help='Directory path to save the files. If none specified, saves in local directory.')
    parser.add_argument('-e','--easyoutput',  help='Produces a simple output file. Made for easy read-ins for other programs.' ,  action = 'store_true')
    parser.add_argument('-g','--globalseq',  help='Parsing in global files (.popl files from ShoRAH output) or generic fasta sequences.',  action = 'store_true')

    if len(sys.argv) <= 3:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()       
 
    #Check the scope input that the form is correct and throws error if user tries to put 89d0 - )*^8 as a region    
    cleanScope = re.sub('\s+',  '',  args.scope)
    refScopeTest = re.search(r'^\d+-\d+$', cleanScope)
    if(not refScopeTest):        
        raise InputError(args.scope, "Invalid Region: ")
    refScope = re.search(r'^(\d+)-(\d+)$', cleanScope) 
    if(int(refScope.group(1)) > int(refScope.group(2))):
        raise InputError(args.scope,  "Invalid Region: ")
    if(args.starting_codon < 1):
        raise InputError(args.starting_codon,  "Start codon position less than 1: ")    
    if(args.starting_codon > int(refScope.group(1))):
        raise InputError(args.starting_codon,  "Invalid Coding region and start codon: ")
    #Assigning values to variables. If the user specifies, then it will change
    #If it is not specified, it is left on it's default value
    if(args.threshold):
        Threshold = args.threshold
    if(args.mutationcombination):
        if(not re.search(r'^\d+-\d+$', args.mutationcombination)):
            raise InputError(args.mutationcombination,  "Invalid Range: ")
        span = re.search(r'^(\d+)-(\d+)$', args.mutationcombination)
        if(int(span.group(1)) > int(span.group(2))):
            raise InputError(args.mutationcombination,  "Reversed Range: ")                
        MCombinationMin = int(span.group(1))
        MCombinationMax = int(span.group(2))        
    if(args.directory):
        if(os.path.isdir(args.directory)):            
            SaveDirectory = re.sub('/+$', '',  args.directory)
        else:
            raise InputError(args.directory,  "Invalid directory: ") 
    if(args.easyoutput):            
        EasyOutput = 'Y'
        #print EasyOutput
    if(args.globalseq):        
        GlobalRun = 'Y'
        #print GlobalRun
    if(args.globalfrequency):
        GThreshold = float(args.globalfrequency)

    processDir(args.fastafilename,  args.referencefilename,  args.scope,  args.starting_codon)    

