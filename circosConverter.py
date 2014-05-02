#!/usr/local/bin/python

#Version 0.001-



import re,  os,  math,  argparse
import sys


#========================"Global Variables and Data Structures"=======================#

AvgReadHash = {}
SaveDirectory = '.'
Shift = 0

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

#================================"printLinks"================================#
#Prints in Link Circos format.

def printLinks(identity,  list, start,  end,  fileName):    

    spaceSwitch = 1
    writeFile = open(SaveDirectory + '/' + fileName+".links",  'w')
    for item in list:
        for x in range (0, len(item)):
            y = x + 1
            if(y >= len(item)):
                if(spaceSwitch == 0):
                    writeFile.write( '\n')
                    spaceSwitch = 1
                break
            if(item[x] >= start and item[y] <= end):    
                #print identity,  item[x],  item[x],  identity,  item[y],  item[y]
                #print identity, int(item[x])-Shift, int(item[x])-Shift, identity, int(item[y])-Shift, int(item[y])-Shift
                writeFile.write(identity + " ")
                writeFile.write(str(int(item[x])-Shift) + " ")
                writeFile.write(str(int(item[x])-Shift) + " ")
                writeFile.write(identity + " ")
                writeFile.write(str(int(item[y])-Shift) + " ")
                writeFile.write(str(int(item[y])-Shift) + "\n")
                spaceSwitch = 0
    writeFile.close()
    
#================================"printHistogram"================================#
#Prints in histogram Circos format.

def printHistogram(identity,  hashTable, start,  end,  fileName):    
    writeFile = open(SaveDirectory + '/' + fileName+".histogram",  "w")
    for keys in hashTable.iterkeys():
        if(int(keys) >= start and int(keys) <= end):
            #print identity, int(keys)-Shift,  keys,  hashTable[keys]
            writeFile.write(identity + " ")
            writeFile.write(str(int(keys)-Shift) + " ")
            writeFile.write(str(int(keys)-Shift) + " ")
            writeFile.write(str(hashTable[keys]) + "\n")
    writeFile.close()
    
#================================"Main Function"================================#
#Parses EasyOutput files and splits them to extract data.
#

def main():    
    global SaveDirectory
    global Shift
    megaList = []
    histogramHash = {}
    regionList = []
    
    program_function = """
    
    *****
    
        CircosConverter
            
            Ver. 0.01
            
            
            This tool takes in EasyOutput files and puts them into Circos format 
            for creating linkages within a single "Chromosome".
            
            
            Handles only EasyOutput files from SMuPFi.
            
            
    *****
            
    """
    parser = argparse.ArgumentParser()
    #required arguments
    parser.add_argument('-f', '--folder', help = 'Folder of EasyOutput files.',  required = True)
    parser.add_argument('-r',  '--region', help = 'The start and stop position of region of interest. E.g. 500-970',  required = True)
    parser.add_argument('-i', '--identity',  help = 'Identity relevent to circos identiy for genome.',  required = True)
    parser.add_argument('-n', '--output_name',  help = 'Name for output file.',  required = True)
    
    #optional
    parser.add_argument('-d',  '--directory',  help = 'Directory path for saving location.')
    parser.add_argument('-sh', '--shift', help='Deduct the position by number.', type = int)
    parser.add_argument('-shared',  '--easyShare', help = 'Switch to look for EasyOutputShared files. [Default looks for EasyOutputUnique]' , 
                                        action = 'store_true')
    parser.add_argument('-li', '--links',  help = 'Produce circos link format file.',  action = 'store_true')
    parser.add_argument('-hi',  '--histogram',  help = 'Produce circos histogram format file.',  action = 'store_true')
    parser.add_argument('-lt',  '--logten',  help = 'Produces log ten of frequency value',  action = 'store_true')
    if len(sys.argv) <= 1:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()
    
    if(args.directory):        
        if(os.path.isdir(args.directory)):            
            SaveDirectory = re.sub('/+$', '',  args.directory)            
        else:
            raise InputError(args.directory,  "Invalid directory: ") 
    if(args.region):
        regionList = re.split('-',  args.region)
    if(args.shift):
        Shift = args.shift - 1
    if(args.easyShare):
        wantFile = 'EasyOutputShared.txt'
    else:
        wantFile = 'EasyOutputUnique.txt'

    if(os.path.isdir(args.folder)):
        fileList = os.listdir(args.folder)        
        for item in sorted(fileList):
            
            reg = wantFile +'$'
            
            if(not re.search(reg,  item)):                
                continue

            openFile = open(os.path.join(args.folder,  item))   
            for line in openFile:                    
                if(re.search('^Common differences not found.$',  line)):
                    continue
                intSplit2 = []                
                split = re.split('\|',  line)
                positions = re.split(' ',  split[0])                    
                if(len(positions) == 1 and positions[0] == ''):
                    continue
                if(args.links):
                    for item in positions:
                        intSplit2.append(int(item))                
                    megaList.append(intSplit2)
                if(args.histogram):
                    for item in positions:
                        if(item not in histogramHash):
                            if(args.logten):
                                histogramHash[item] = math.log10(float(split[2]))
                            else:
                                histogramHash[item] = float(split[2])
                        else:
                            if(args.logten):
                                if(histogramHash[item] + math.log10(float(split[2])) > 2.0):
                                    histogramHash[item] = 2
                                else:
                                    histogramHash[item] += math.log10(float(split[2]))
                            else:
                                if(histogramHash[item] + float(split[2]) > 1.0):
                                    histogramHash[item] = 1.0
                                else:
                                    histogramHash[item] += float(split[2])
            openFile.close
            
            
        if(args.links):        
            printLinks(args.identity,  megaList,  int(regionList[0]),  int(regionList[1]), args.output_name)
            print "Links done."
        if(args.histogram):
            printHistogram(args.identity, histogramHash, int(regionList[0]), int(regionList[1]),  args.output_name)
            print "Histograms done."
    else:
        raise InputError(args.folder,  "Is not a valid directory or file. Given argument: ")
        
    
#==========================="Start of Script"=================================#

main()




