#!/usr/local/bin/python

import re, os, argparse
import sys, math

#======================"Main"======================#

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #required argument
    parser.add_argument('-f', '--file', help = 'EasyOutputShared file to be used.',  required = True)
    parser.add_argument('-o',  '--output_name',  help = 'Output file name.',  required = True)
    parser.add_argument('-t', '--threshold',  help ='Cut-off threshold',  required = True,  type = float)
    #optional
    parser.add_argument('-d',  '--directory',  help = 'Save directory path.')
    
    
    if len(sys.argv) <= 1:
            #print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()
    
    if(os.path.isfile(args.file)):
        fileHandle = open(args.file)
        fileOutput = args.output_name + ".sif"
        fileOutputHandle = open(fileOutput,  'w')
        for line in fileHandle:
            if(float(line.split()[2]) < args.threshold):
                continue
            node1 = line.split()[0]
            node2 = line.split()[1]
            fileOutputHandle.write(node1 + "\tlinks\t" + node2 + "\n")
#            if(float(line.split()[2]) < args.threshold):             
#                fileOutputHandle.write(node1 + "\n")
#                fileOutputHandle.write(node2 + "\n")
#            else:
#                fileOutputHandle.write(node1 + "\tlinks\t" + node2 + "\n")
            
        fileHandle.close()
        fileOutputHandle.close()
    else:
        print "Not valid file."
        sys.exit(1)
        
