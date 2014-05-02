#README - Version 0.02

This is a readme file for the SMuPFi Package that contains guidelines and important notes about
all the tools in this repository.

For the pipeline to work, these dependencies are required beforehand:

1. Python2.7

2. Biopython

3. Numpy

4. Scipy

5. Perl 5

**All these tools are freely available and users can modify and redistribute tools for their own purposes. The source code of the tools are provided within the repository**

##Tools

######==a_smupfi.py==

a-smupfi.py is the Shared Mutation Pattern Finder tool that takes in viral sequences in fasta format. It iterates through sequences and identifies all changes against the reference sequence and searches for mutations that are shared among two or more unique
sequences.

**Usage and Options:**

There are many options implemented in the tool, and these can be viewed by running:

`a_smupfi.py -h`

The tool has **-f FASTAFILENAME**, **-r REFERENCEFILENAME** **-s SCOPE** and **-sc START_CODON** as
required inputs where these are the basic details for searching shared mutations.

General usage of a-smupfi will require the options **-g -gf FREQUENCY -m RANGE**, where -g signifies it is either global ShoRAH outputs which have .popl endings or generic fasta sequences. -gf is a decimal number tells a-smupfi the threshold at which sequences below the threshold are not included in the shared mutation search and -m is a mutation combination range (e.g 1-2 or 4-4) meaning it will search for combination of patterns from length 1 to 2 or length 4 to 4 respectively. 

Usage:

	a_smupfi.py [-h] -f FASTAFILENAME -r REFERENCEFILENAME -s SCOPE -sc STARTING_CODON
	            [-m MUTATIONCOMBINATION] [-t THRESHOLD] [-gf GLOBALFREQUENCY]
	            [-d DIRECTORY] [-e] [-g]
**Important:**

-m option is crucial if the data set has a very diverse sequence and very long in length. It will limit the amount of processing power and memory required. If a very long sequence (>1500 amino acids or 4500 nucleotides), it is recommended to either decrease -m or use a-smupfi with only sections of the data set.

If the fasta sequence files that going to be used as input data does not resemble:

	>sequenceName_0.02
	ATGGTCTCATCTAGGCCAATG

where 0.02 signifies the frequency, use uniqvariant.py (further explained below in this README file) prior to a-smupfi.py.

***Outputs***

a-smupfi provides four files for output (XXXX-YYYY is the region specified by -s).

1. XXXX-YYYY.txt => shows a graphical representation of where the mutation patterns occur

2. XXXX-YYYY_Summary.txt => shows only shared mutations between unique sequences

(if -e is given):

3. XXXX-YYYY_EasyOutputShared.txt => a file made for easy parsing, same as Summary.txt

4. XXXX-YYYY_EasyOutputUnique.txt => a file made for easy parsing, but selects all mutations on any sequence

###### ==choplqb.py==

choplqb.py is a simple cleaning tool that uses average phred-score over a user given window length to proccess raw fasta sequences. The tool trims the reads when the average phred-score falls below the user given threshold It requires two files, a .fasta or .fna file and a .qual file. This is usually used with Next Generation Sequencing (NGS) data where the accuracy and trusthworthiness of individual reads may decrease at the near the end of the read.

Usage:

	choplqb.py [-h] -f FASTAFILENAME -q QUALITYFILENAME -w WINDOW_SIZE 
	            -t QUALITY_THRESHOLD [-d] [-D OUTPUT_DIR]

######==circosConverter.py==

circosConverter.py is a format converter tool that simply changes easyoutput files from a-smupfy.py to a format that Circos can understand. This tool is used for creating the histogram `-hi`as well as links `li` around the circos plots.

Usage:
	
	circosConverter.py [-h] -f FOLDER -r REGION -i IDENTITY -n OUTPUT_NAME
                       [-d DIRECTORY] [-sh SHIFT] [-shared] [-li] [-hi]
                       [-lt]
	


######==indelRemover.py==

A tool to remove indels from fasta sequences and replace it with reference sequence characters and merges fasta sequences that end up being identical.

In detail, due to 454 technology having a tendency to have high errors in the form of '-' characters, also known as 'indels', around homopolymer regions, this tool specifically targets indels around those regions and substitutes the indels with the character at the same position from the reference sequence.

In this tool a homopolymer region is any residue appearing 3 or more times consecutively. For example, 'GGG'. Hence if an indel is detected to be around a homopolymer region like 'G-G', '-GG' or 'GG-', the tool will proceed to clean that indel.

Usage:

	indelremover.py [-h] -f FASTAFILENAME -r REFERENCEFILENAME -s SCOPE 
			        [-fr FREQ_INDEX] [-o OUTPUT_FILE] [-a]

######==jaccardToCircos.py==

jaccardToCircos.py is a simple format converter tool that deals with jaccard easyoutput files and converts it to a format, which Circos can read and make arcs.

Usage:

	jaccardToCircos.py [-h] -f FILE -i IDENTITY -o OUTPUT_FILE [-r REGION]
                       [-d DIRECTORY] [-lt] [-li] [-p PVALUE_THRESHOLD]
                       [-T]

######==javssim.py==

javssim.py (Jaccard Value Statistical Simulator) is a simulator utilising easyOutputShared files from a-smupfi.py to conduct statistical test on paired mutation and returns a ranked list of mutation that identifies which pairs are occurring more than random and which pairs are occurring less than random.

Usage:
	
	javssim.py [-h] -f FILE -n NAME [-d DESTINATION] [-e]

where providing the option `-e` will tell javssim.py to output an easy to parse file similar to the easyoutput files created by a_smupfy.py. All the data are delimited by `|` character.

######==omes.py==

omes.py is an algorithm that uses the observed minus expected squared formula to calculate which pair of of positions are co-varying. It outputs a list of pairs which have been calculated to hold a score greater than 0.

Usage:

	omes.py [-h] -f FILE [-st START_POS]

######==qualfa2fq.pl/qualfa2fq_bt.pl==

qualfa2fq.pl is a merger tool that intergrates a .fasta / .fna file and a .qual file into a .fastq/.fq file.

qualfa2fq_bt.pl is a version doing the same process, but made for aligners (e.g. bowtie2) that dislike the 60 character limit in the .fastq /.fq files.

Script taken from BWA

Citation:
 
* Li H. and Durbin R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler Transform. Bioinformatics, Epub. [PMID: 20080505] 

Download link: [BWA Tool](http://sourceforge.net/projects/bio-bwa/files/)

Home Page: [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/)

**Note:**  Please download qualfa2fq.pl from above web link. You can find qualfa2fq.pl inside the compressed file once extraction is done. qualfa2fq_bt.pl has been very slightly modified by me (commented out some lines) and made available here.

######==sequenceExpand.py==

Uses fasta sequences with frequency information in the fasta tag as input and expands according to a base number. Using ShoRAH output as an example:

	>sequence1_0.3542
	AGTCGTATTTAC

Usage:
	
	sequenceExpand.py [-h] -f FILE -b BASE

If base is given as -b 100, then sequence 1 will be expanded to 35 reads. If -b 1000, then sequence 1 will be expanded as 354.

######==sif_converter.py==

Takes in outputs from omes.py and converts the user defined pairs into .sif format which can be taken in by Cytoscape to draw networks. The tool has the option -t which allows users to define their own threshold for selecting data. 

Usage:

	sif_converter.py [-h] -f FILE -o OUTPUT_NAME -t THRESHOLD
                     [-d DIRECTORY]

######==snpExtractor.py==
This script takes SNPS recorded in an output file from the software LoFreq. It substitute the SNP into a given reference sequence and performs translations to display amino acid mutations.

Download: [LoFreq Tool](http://sourceforge.net/projects/lofreq/)

Home Page: [Lofreq Wiki](http://sourceforge.net/p/lofreq/wiki/Home/)    

**Note:** This tool only takes in 3 columns of data. Any extra information will be ignored. 

If unsure if option -rp is needed, use:

	cat -A FILE 

or 
	
	cat -v FILE

to check whether ^M character exists, if seen, use -rp.

For argument -fup, the tool will ignore frequencies greater than inputed frequency. This is to ignore outputs that fall under the case such as a wrong consensus reference is used to align and results in high instances of mutation recorded.

Usage:

	snpExtractor.py [-h] -f FILE -r REFERENCEFILENAME -d DELIMITER -p POS_COLUMN
	                 -m MUT_COLUMN -fq FREQ_COLUMN -sc STARTING_CODON [-e]
                    [-o OUTPUT_FILE] [-i] [-rp] [-lt] [-fup FILTER] [-sf] [-is]


######==snpExtractConverter.py==

snpExtractConverter.py changes the format of the output of the snpExtractor.py into histogram format that can be understood by Circos.

Usage:

	snpExtractConverter.py [-h] -f FILE -r REGION -i IDENTITY [-d DIRECTORY] 
	                       [-sh SHIFT]

######==uniqvariant.py==

uniqvariant.py takes in fasta files that have (potentially) repeated sequences and merges them into a unique sequence and provides that sequence with a frequency of occurrence number. For example with 3 sequences:

	>Sequence1
	AGTCGTG
	>Sequence2
	AGTCGTG
	>Sequence3
	AGTCGTC

will become

	>Sequence1_0.66
	AGTCGTG
	>Sequence2_0.33
	AGTCGTC

If your sequences have frequency information, such as:

	>Sequence1_0.50
	AGTCGTG
	>Sequence2_0.25
	AGTCGTG
	>Sequence3_0.25
	AGTCGTC

then it can add it to become:

	>Sequence1_0.5
	AGTCGTG
	>Sequence2_0.5
	AGTCGTC

Usage:

	uniqvariant.py [-h] -f FASTAFILENAME [-i FREQ_INDEX] [-o OUTPUT_FILE]
                   [-p PREFIX]

##Data Files

######==Circos TEST DATA== 

To construct Circos plots, you need Circos installed.

Download: [Circos Tool](http://circos.ca/software/download/circos/)

Home Page: [Circos Website](http://circos.ca/)


There is a Circos Test data folder in the repository name `Ch240_Time4-CircosTESTDATA` which holds the necessary files to construct a sample Circos plot. The files that requires editing are:

	config_240_t4.conf
		
In `config_240_t4.conf`, you will need to provide the correct path to where it can find the data files. These files are located in the `Data` folder, where it holds all files that are generated using the tools in the pipeline with some manual refinement to select specific desired data.	
	
Other files within the `Ch240_Time4-CircosTESTDATA` are cosmetic parameters that can be altered to personalise how the Circos plot is drawn. See [Circos Tutorials](http://circos.ca/tutorials/lessons/) to see how to edit these parameters.


######==TEST DATA==

Sample viral sequences reconstructed through ShoRAH and SNP data generated from Lofreq has been provided to act as a sample run for the pipeline. They are located in the folder `TESTDATA`, with subfolders `LoFreq_Output` and `ShoRAH_Viral_Sequences`. 


###Authors
:man_with_gua_pi_mao: Preston Leung (preston@student.unsw.edu.au)

:man: Fabio Luciani (luciani@unsw.edu.au) 
