#!/usr/bin/python

#Version 0.001-



import re,  os,  math,  argparse
import sys


#========================"Global Variables and Data Structures"=======================#

AvgReadHash = {}
#OccurrenceHash = {}
SaveDirectory = ''
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

#================================"printSummary"================================#
#Prints a summary of the parsed files.

def printSummary(identity,  list, start,  end):    

    spaceSwitch = 1
    for item in list:
        for x in range (0, len(item)):
            y = x + 1
            if(y >= len(item)):
                if(spaceSwitch == 0):
                    print '\n'
                    spaceSwitch = 1
                break
            if(item[x] >= start and item[y] <= end):    
                #print identity,  item[x],  item[x],  identity,  item[y],  item[y]
                print identity, int(item[x])-Shift, int(item[x])-Shift, identity, int(item[y])-Shift, int(item[y])-Shift
                spaceSwitch = 0


#================================"Main Function"================================#
#Parses EasyOutput files
#The split list are always 4 elements long. ALWAYS.

def main():    
    global SaveDirectory
    global Shift
    megaList = []
    regionList = []
    
    program_function = """
    *****This tool takes in EasyOutput files and summarises. Handles only EasyOutput files from SMuPFi.*****
    """
    parser = argparse.ArgumentParser()
    #required arguments
    parser.add_argument('-f', '--folder', help = 'Folder of EasyOutput files.',  required = True)
    parser.add_argument('-r',  '--region', help = 'The genomic region of interest. E.g. 500-970',  required = True)
    parser.add_argument('-i', '--identity',  help = 'Identity relevent to circos identiy for genome',  required = True)
    
    
    #optional
    parser.add_argument('-d',  '--directory',  help = 'Directory path for saving location')
    parser.add_argument('-sh', '--shift', help='Deduct the position by number', type = int)
    if len(sys.argv) <= 0:
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
    
    if(os.path.isdir(args.folder)):
        fileList = os.listdir(args.folder)
        for item in sorted(fileList):
            ###                                                                              ###
            #   Switch this manually for now depending if you're using shared or unique files  #
            ###                                                                              ###
            if(not re.search('EasyOutputUnique.txt$',  item)):
                continue
            try:
                openFile = open(os.path.join(args.folder,  item))    
            except IOError:
                print "Error: cannot find {}".format(item)
            else:     
                for line in openFile:                    
                    if(re.search('^Common differences not found.$',  line)):
                        continue
                    intSplit2 = []
                    split = re.split('\|',  line)
                    if(len(split[0]) > 0):
                        positions = re.split(' ',  split[0])                    
                        for item in positions:
                            intSplit2.append(int(item))                    
                    
                    megaList.append(intSplit2)                    
                openFile.close        
        printSummary(args.identity,  megaList,  int(regionList[0]),  int(regionList[1]))         
    else:
        raise InputError(args.folder,  "Is not a valid directory or file. Given argument: ")
        
    
#==========================="Start of Script"=================================#

main()




