#!/usr/local/bin/python
#
#JaVSSimDupplet:
#   Version 0.06

import random
import math,  sys, gc
import re,  os,  argparse
import numpy.random as nprd
import numpy as np
import scipy.stats as scistat

#========================"Global Variables"=======================#

#File name for the jaccard file.
FileName = ''
#Directory path for saving jaccard file.
SavePath = ''
#Easyoutput flag
MakeEasy = False

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

#========================"everyDayImShuffling"========================#
#Input:
#   - A two dimensional array.
#Proccess:
#   - Shuffles the elements inside the vectors of mDList and returns mDList.

def everyDayImShuffling(mDList):
    
    for i in range(0,  len(mDList)):
        mDList[i] = random.sample(mDList[i],  len(mDList[i]))
    return mDList    

#========================"DoJaccard"========================#
#Calculates Jaccard value.
#
#Input:
#   - a two dimensional list (mDlist). A list containing lists of vectors (also lists).
#   - the length of one vector (vectorLen).
#Process:
#   - Iterates through mDlist and calculates the Jaccard value of how many times
#     it observes bothXY happening in ratio to all combinations of onlyX, onlyY and bothXY.
#   - Returns the jaccard value.

def doJaccard(mDList,  vectorLen):    
    
    jValue = 0.0
    onlyX,  onlyY, bothXY= (0.0,  0.0,  0.0)
    
    for column in range(0,  vectorLen):        
        if(mDList[0][column] == 1 and mDList[1][column] == 1):
           bothXY += 1.0 / float(vectorLen)
        elif(mDList[0][column] == 1 and not mDList[1][column] == 1):
            onlyY += 1.0 / float(vectorLen)
        elif(not mDList[0][column] == 1 and mDList[1][column] == 1):
            onlyX += 1.0 / float(vectorLen)
    
    
    assert(bothXY >= 0.0)
    assert(onlyY >= 0.0)
    assert(onlyX >= 0.0)
    if(not(bothXY + onlyY + onlyX) == 0.0):
        jValue = bothXY / (bothXY + onlyY + onlyX)
    else:
        jValue = 0.0

    return jValue

#========================"jackKnife"========================#
#This function performs the jackKnife procedure and produces a standard error of 
#the mean (SEM).
#
#Input: 
#   - a two dimensional list (mDlist). A list containing lists of vectors (also lists).
#   - the length of one vector (vectorLen).
#   - the original observed Jaccard value before shuffling (obsJ).
#Process:
#   - Iterates through mDlist and sequentially omit one pair of values in mDlist to 
#     calculate new psuedo Jaccard values. 
#   - These new Jaccard values are then used to calculate standard deviation and 
#     returns the SEM.

def jackKnife(mDList,  vectorLen,  obsJ):
    
    pList = []
    vectorJ = []
    p_mean,  pSqSum, pStd_Deviation,  pStd_Err = (0.0,  0.0,  0.0,  0.0)
    jValue = 0.0
    onlyX,  onlyY, bothXY= (0.0,  0.0,  0.0)
        
    gc.disable()
    for i in range(0,  vectorLen):
        for column in range(0,  vectorLen):
            if (column == i):
                continue
            else:
                if(mDList[0][column] == 1 and mDList[1][column] == 1):
                    bothXY += 1.0 / float(vectorLen)
                elif(mDList[0][column] == 1 and not mDList[1][column] == 1):
                    onlyY += 1.0 / float(vectorLen)
                elif(not mDList[0][column] == 1 and mDList[1][column] == 1):
                    onlyX += 1.0 / float(vectorLen)
        
        #print "bothXY:{}, onlyY:{}, onlyX:{}".format(bothXY, onlyY, onlyX)
        
        assert(bothXY >= 0.0)
        assert(onlyY >= 0.0)
        assert(onlyX >= 0.0)
        if(not (bothXY + onlyY + onlyX) == 0.0):
            jValue = bothXY / (bothXY + onlyY + onlyX)
        else:
            jValue = 0.0
        p = (vectorLen * obsJ) - (vectorLen - 1) * jValue        
        pList.append(p)
    gc.enable()
    assert(len(pList) == vectorLen)    
    pStd_Deviation = np.std(pList,  ddof=1)    
    pStd_Err = pStd_Deviation / math.sqrt(vectorLen)     
    return pStd_Err

#========================"fillVectors"========================#
#This function fills the two vectors with 1 (yes mutation) and 0 (no mutation)
#to the size of vectorLen.
#
#Input:
#   - frequencies x, y and xy. e.g. 0.07.
#   - vectorLen is the size the vector used for sampling.
#Process:
#   - if xy is 0.07, in a vectorLen of 1000, both vectorX and vectorY
#     will contain 1s for 7 slots in the array then returns the array.  

#   vectorInfo = [onlyXFreq,  onlyYFreq, mutXYZFreq, vectorSize]

def fillVectors(vList):
    vectorX,  vectorY = ([], [])    
    vectorLen = int(vList[3])
    modiX = int(vList[0] * vectorLen)
    modiY = int(vList[1] * vectorLen)    
    modiXY = int(vList[2] * vectorLen)
    
    start = 0
    
    gc.disable()
    
        
    for i in range(start,  start+modiXY):
        vectorX.append(1)
        vectorY.append(1)        
    start += modiXY
    
    for i in range(start, start+modiX):
        vectorX.append(1)
        vectorY.append(0)        
    start += modiX    
    
    for i in range(start, start+modiY):
        vectorX.append(0)
        vectorY.append(1)       
    start += modiY   

    for i in range (start,  vectorLen):
        vectorX.append(0)
        vectorY.append(0)
    
    gc.enable()
    

    return [vectorX,  vectorY]    

#========================"openFile"========================#
#This function opens the file and selects which lines to proccess
#using given arguments.
#
#Input:
#   - a list containing file names (fileList).
#   - a directory path of where to find the folder containing
#     the file (folderPath).
#   - the length of pattern required for doing jaccard (mCombination)
#Proccess:
#   - opens file, extracting required information and returns the 
#     extracted information (lineList).

def openFile(fileName, mCombination):       
    
    lineList = []    
        
    try:
        openFile = open(fileName)    
    except IOError:
        print "Error: cannot find {}".format(fileName)
    else:
        for line in openFile:                    
            if(re.search('^Common differences not found.$',  line)):                            
                break
            split = re.split('\|',  line)
            position = re.split(' ',  split[0])
            if(not (len(position) <= mCombination)):
                continue
            lineList.append(line)
    openFile.close()
    return lineList

#========================"dataExtraction"========================#
#This function picks out individual information and stores them for easy accessing.
#
#Input:
#   - A

def dataExtraction(lineList, mCombination):
    
    single, double = ({}, {})
    
    for line in lineList:
        split = re.split('\|',  line)
        position = re.split(' ',  split[0])
        pattern = re.split(' ',  split[1])
        refPat = re.split(' ', split[5])
        freq = "%0.3f" % float(split[2])      
        
        if(len(position) == mCombination):
            key = refPat[0].rstrip('\n') + str(position[0]) + pattern[0] + "|" + refPat[1].rstrip('\n') + str(position[1]) + pattern[1]
            double[key] = freq
        elif(len(position) == 1):
            key = refPat[0].rstrip('\n') + str(position[0]) + pattern[0]
            single[key] = freq
    
    return single, double

#========================"printStats"========================#
def printStats(indicies, mutXList, onlyXList, mutYList, onlyYList,
                         xyList, jList, jrList, jsList, pvList, zsList): 
    
    if(SavePath):
        filepath = os.path.join(SavePath, FileName)
    else:
        filepath = FileName 
    
    try:
        output = open(filepath + '.txt', 'w')
        output.write("MutX\tX-OnlyFreq\tMutY\tY-OnlyFreq\tXYFreq\tJ\tJ_rand\tJ_se\tP\tZ\n")
        
        for i in range(len(indicies)-1, -1, -1):            
        
            output.write(str(mutXList[indicies[i]]) + "\t" + str(onlyXList[indicies[i]]) + "\t")
            output.write(str(mutYList[indicies[i]]) + "\t" + str(onlyYList[indicies[i]]) + "\t")            
            output.write(str(xyList[indicies[i]]) + "\t" + str(jList[indicies[i]]) + "\t")
            output.write(str(jrList[indicies[i]]) + "\t" + str(jsList[indicies[i]]) + "\t")
            output.write(str(pvList[indicies[i]]) + "\t" + str(zsList[indicies[i]]) + "\n")

    
        output.close()
        print "{} done.".format(FileName)
    except IOError:
        print "Cannot write file {}.txt".format(FileName)

#========================"printEasyOutput"========================#

def printEasyOutput(indicies, mutXList, onlyXList, mutYList, onlyYList,
           xyList, jList, jrList, jsList, pvList, zsList):

    if(SavePath):
        filepath = os.path.join(SavePath, FileName)
    else:
        filepath = FileName
    try:
        output = open(filepath + '_EasyOutput.txt', 'w')
        for i in range(len(indicies)-1, -1, -1):
            
    
            output.write(str(mutXList[indicies[i]]) + '|')
            output.write(str(onlyXList[indicies[i]]) + '|')
        
            output.write(str(mutYList[indicies[i]]) + '|')
            output.write(str(onlyYList[indicies[i]]) + '|')
        
            output.write(str(xyList[indicies[i]]) + '|')
            output.write(str(jList[indicies[i]]) + '|')
        
            output.write(str(jrList[indicies[i]]) + '|')
            output.write(str(jsList[indicies[i]]) + '|')
        
            output.write(str(pvList[indicies[i]]) + '|')
            output.write(str(zsList[indicies[i]]) + '\n')
    
        output.close()
        print "{} done.".format(FileName + '_EasyOutput.txt')
    except IOError:
        print "Cannot write file {}".format(FileName + '_EasyOutput.txt')



#========================"Main"========================#

if __name__ == '__main__':
    
    vectorSize = 1000
    mCombination = 2
    allVector, copyVector, vectorInfo = ([],[],[])
    mutXList, mutYList, jList = ([],[],[])
    jrList,jsList, zsList, pvList = ([],[],[],[])
    onlyXList, onlyYList, XYList = ([],[],[])
    rearrangement_num = 2000
    j_rand ,  j_se,  z_score = (0.0,  0.0,  0.0)
    single, double, triple = ({},{},{})
    relationCheck = set()
    program_function = """
    
            Jaccard Value Statistic Simulator (JaVSSim)            
            
            - Duplet Version
            
            - A simple simulator utilising easyOutputShared files to 
              conduct statistical test on paired mutation.            
            
            - Ver. 0.06
                -> Requires Numpy and Scipy and Python2.7.
                -> 
    """
      
    parser = argparse.ArgumentParser()
    #required argument
    parser.add_argument('-f', '--file', help = 'EasyOutputShared file to be used.',  required = True)    
    parser.add_argument('-n', '--name', help = 'Name of the output file.', required = True)
    
    #optional arguments
    parser.add_argument('-d', '--destination', help = 'Saving location.')   
    parser.add_argument('-e', '--easyoutput', help = 'Output an easy parsing file delimited by "|" and header removed.', action ='store_true')
    
        
    if len(sys.argv) <= 1:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()
    
    if(os.path.isfile(args.file)):
        if(not re.search('EasyOutputShared.txt$',args.file)):
            raise InputError(args.file, "File is not an EasyOutputShared file: ")        
        tempList = openFile(args.file, mCombination)
        single, double = dataExtraction(tempList, mCombination)
        del tempList
    else:
        raise InputError(args.file, "Invalid file: ")   
    if(args.destination):
        if(not os.path.isdir(args.destination)):
            raise InputError(args.destination, "Invalid saving destination: ")
        else:
            SavePath = args.destination     
    FileName = re.sub('\..*$', '', args.name)
    if(args.easyoutput):        
        MakeEasy = True
    
    
    for doubleMutXY in double.iterkeys():
        singleMutX = doubleMutXY.split('|')[0]
        singleMutY = doubleMutXY.split('|')[1]        

        mutXYFreq = float(double[doubleMutXY])
        onlyXFreq = float(single[singleMutX]) - mutXYFreq
        onlyYFreq = float(single[singleMutY]) - mutXYFreq        

        vectorInfo = [onlyXFreq,  onlyYFreq, mutXYFreq, vectorSize]        

        allVector = fillVectors(vectorInfo)
        #print len(allVector[0])
        copyVector = allVector[:]
        obsJ = doJaccard(allVector,  vectorSize)
        
        jVectorList = []
        for i in range (0,  rearrangement_num):
            allVector = everyDayImShuffling(allVector)
            jVectorList.append(doJaccard(allVector,  vectorSize))
            
        j_rand = np.mean(jVectorList)
        j_se = jackKnife(copyVector,  vectorSize,  obsJ)
       
       #Making j_se not 0, but a tiny number
        if(j_se == 0):                    
            z_score = (obsJ - j_rand) / 0.000000000000001
        else:
            z_score = (obsJ - j_rand) / j_se       
        
        p_value = scistat.norm.sf(z_score) #one-sided
        
        mutXList.append(singleMutX)
        onlyXList.append(onlyXFreq)
        mutYList.append(singleMutY)
        onlyYList.append(onlyYFreq)        
        XYList.append(mutXYFreq)
        jList.append(math.ceil(obsJ * 1000.0)/1000.0)
        jrList.append(math.ceil(j_rand * 1000.0)/1000.0)
        jsList.append(math.ceil(j_se * 1000.0)/1000.0)
        zsList.append(math.ceil(z_score * 1000.0)/1000.0)                
        pvList.append(math.ceil(p_value * 1000.0)/1000.0)
        
    indicies = np.lexsort(keys = (pvList, zsList))

    printStats(indicies, mutXList, onlyXList, mutYList, onlyYList,
                    XYList, jList, jrList, jsList, pvList, zsList)
    if(MakeEasy):
        printEasyOutput(indicies, mutXList, onlyXList, mutYList, onlyYList,
                    XYList, jList, jrList, jsList, pvList, zsList)
        
