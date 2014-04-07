#!/usr/bin/python
# Ver. 0.002
#renamed from jaccardToCircos2 -> jaccardToCircos

import re,  os,  math,  argparse
import sys

SaveDir = ''
#Z_Threshold = 2.325
P_Threshold = 0.05
HasRegion = 0
DoLogTen = False
#================================"filterAll"================================#

def filterAll(mutX,  mutXFreq,  xy_only_freq):
    
    mutX = re.sub('\s',  '',  mutX)
    
    #print mutX
    if (DoLogTen):
        freq = math.log10((xy_only_freq + mutXFreq)*100)
    else:
        freq = xy_only_freq + mutXFreq
    
    return (re.search('\s*[a-z*-]\d+([a-z*-])',  mutX,  re.I).group(1), #Mutant
                 freq,  
                 re.search('\s*[a-z*-](\d+)[a-z*-]',  mutX,  re.I).group(1)) #Position
    
    

#================================"makeLinks"================================#

def makeLinks(openFile, identity,  outputFileName, region = 'N'):
    
    outHandle = open(SaveDir + '/' + outputFileName+".links",  'w')
    for line in openFile:
        
        linesplit = re.split('\|',  line)        
        
        if(float(linesplit[8]) > P_Threshold):
            continue
                
        posX = re.search('\s*[a-z*-](\d+)[a-z*-]',  linesplit[0].strip(),  re.I)
        posY = re.search('\s*[a-z*-](\d+)[a-z*-]',  linesplit[2].strip(),  re.I)
        if(HasRegion):
            if(int(posX.group(1))>= int(region.group(1)) and int(posY.group(1)) <=int(region.group(2))):
                outHandle.write(identity + " ")
                outHandle.write(posX.group(1) + " " + posX.group(1) + " ")
                outHandle.write(identity + " ")
                outHandle.write( posY.group(1)+" "+ posY.group(1)+"\n")
        else:
            outHandle.write(identity + " ")
            outHandle.write(posX.group(1) + " " + posX.group(1) + " ")
            outHandle.write(identity + " ")
            outHandle.write( posY.group(1)+" "+ posY.group(1)+"\n")
    
    outHandle.close()        
    print "Links Done."    
    
#================================"makeLinks_T"================================#

def makeLinks_T(openFile, identity,  outputFileName, region = 'N'):
    
    outHandle = open(outputFileName+".links",  'w')
    for line in openFile:
        
        linesplit = re.split('\|',  line)        
        
        if(float(linesplit[10]) > P_Threshold):
            continue
                
        posX = re.search('\s*[a-z*-](\d+)[a-z*-]',  linesplit[0].strip(),  re.I)
        posY = re.search('\s*[a-z*-](\d+)[a-z*-]',  linesplit[2].strip(),  re.I)
        posZ = re.search('\s*[a-z*-](\d+)[a-z*-]',  linesplit[4].strip(),  re.I)
        if(HasRegion):
            if(int(posX.group(1))>= int(region.group(1)) and int(posZ.group(1)) <=int(region.group(2))):                
                outHandle.write(identity + " ")
                outHandle.write(posX.group(1) + " " + posX.group(1) + " ")
                outHandle.write(identity + " ")
                outHandle.write(posY.group(1)+" "+ posY.group(1)+"\n")
                outHandle.write(identity + " ")
                outHandle.write(posY.group(1) + " " + posY.group(1) + " ")
                outHandle.write(identity + " ")
                outHandle.write(posZ.group(1)+" "+ posZ.group(1)+"\n\n")
        else:            
                outHandle.write(identity + " ")
                outHandle.write(posX.group(1) + " " + posX.group(1) + " ")
                outHandle.write(identity + " ")
                outHandle.write(posY.group(1)+" "+ posY.group(1)+"\n")
                outHandle.write(identity + " ")
                outHandle.write(posY.group(1) + " " + posY.group(1) + " ")
                outHandle.write(identity + " ")
                outHandle.write(posZ.group(1)+" "+ posZ.group(1)+"\n\n")
    
    outHandle.close()        
    print "Links Done."    

#================================"Main Function"================================#

def main():    
    global SaveDir    
    global HasRegion
    global P_Threshold
    global DoLogTen
        
    program_function = """
    *****
        JaccardToCircos
            - Ver. 2
            - Deals with jaccard easyoutput files and converts it to a format 
              that circos can read and make arcs.
            - 
    
    *****
    """
    parser = argparse.ArgumentParser()
    #required arguments
    parser.add_argument('-f', '--file', help = 'jaccard easy output files.',  required = True)    
    parser.add_argument('-i', '--identity',  help = 'Identity relevent to circos identiy for genome',  required = True)
    parser.add_argument('-o',  '--output_file',  help = "Name for the output files.",  required = True)
    
    #optional
    parser.add_argument('-r',  '--region', help = 'The genomic region of interest. E.g. 500-970')
    parser.add_argument('-d',  '--directory',  help = 'Directory path for saving location')    
    parser.add_argument('-lt',  '--logten',  help='Histogram values are log ten',  action = 'store_true')
    parser.add_argument('-li', '--links',  help='Directs tool to produce links format.',  action = 'store_true')
    parser.add_argument('-p',  '--pValue_threshold',  help = 'P-value cut off. Default 0.05.',  type = float)
    parser.add_argument('-T',  '--triplet_file',  help = 'Input file is a jaccard triplet output',  action = 'store_true')
    
    if len(sys.argv) <= 1:
            print(program_function)
            parser.print_help()
            sys.exit()
    args=parser.parse_args()
    
    if(args.directory):        
        if(os.path.isdir(args.directory)):            
            SaveDir = re.sub('/+$', '',  args.directory)            
        else:
            raise InputError(args.directory,  "Invalid directory: ") 
    if(args.pValue_threshold):
        P_Threshold = args.pValue_threshold
    if(args.region):
        regionTest = re.search(r'^\d+-\d+$', args.region)
        if(not regionTest):        
            raise InputError(args.region, "Invalid Region: ")
        region = re.search(r'^(\d+)-(\d+)$', args.region) 
        if(int(region.group(1)) > int(region.group(2))):
            raise InputError(args.region,  "Invalid Region: ")
        HasRegion = 1    
    if(args.logten):
        DoLogTen = True
    
    if(os.path.isfile(args.file)):
        
        if(args.triplet_file):
            
            if(args.links):                
                print "Processing link format for Triplet."
                openFile = open(args.file)        
                if(HasRegion):
                    makeLinks_T(openFile,  args.identity, args.output_file,  region)
                else:
                    makeLinks_T(openFile,  args.identity, args.output_file)
                openFile.close()
        else:   
            
            
            if(args.links):
                openFile = open(args.file)        
                if(HasRegion):
                    makeLinks(openFile, args.identity, args.output_file, region)
                else:
                    makeLinks(openFile, args.identity,  args.output_file)      
                openFile.close()
        
       
        print "Done."
        
    else:
        raise InputError(args.file,  "Is not a valid file. Given argument: ")

    
#==========================="Start of Script"=================================#

main()
