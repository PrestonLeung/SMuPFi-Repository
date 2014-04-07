README - Version 0.01

This is a readme file for the SMuPFi Package that contains guidelines and important notes about
all the tools in this repository.

For the pipeline to work, these dependencies are required beforehand:

1) Python2.7
2) Biopython
3) Numpy
4) Scipy
5) Perl 5

******All these tools are freely available and users can modify tools. ******

================================a_smupfi.py================================

a-smupfi.py is the Shared Mutation Pattern Finder tool that takes in viral sequences in
fasta format. It iterates through sequences and identifies all changes against the 
reference sequence and searches for mutations that are shared among two or more unique
sequences.

***Usage and Options:***

There are many options implemented in the tool, and these can be viewed by running 
'python2.7 a_smupfi.py -h'.

The tool has -f FASTAFILENAME, -r REFERENCEFILENAME -s SCOPE and -sc START_CODON as
required inputs where these are the basic details for searching shared mutations.

General usage of a-smupfi will require the options -g -gf FREQUENCY -m RANGE, where
-g signifies it is either global ShoRAH outputs which have .popl endings or generic fasta 
sequences. -gf is a decimal number tells a-smupfi the threshold at which sequences 
below the threshold are not included in the shared mutation search and -m is a 
mutation combination range (e.g 1-2 or 4-4) meaning it will search for combination 
of patterns from length 1 to 2 or length 4 to 4 respectively. 

---Important---
-m option is crucial if the data set has a very diverse sequence and very long in 
length. It will limit the amount of processing power and memory required. If a very long
sequence (>1500 amino acids or 4500 nucleotides), it is recommended to either decrease -m
or use a-smupfi with only sections of the data set.

If the fasta sequence files that going to be used as input data does not resemble:

>sequenceName_0.02
ATGGTCTCATCTAGGCCAATG

where 0.02 signifies the frequency, use uniqvariant (further explained below in this
README file) prior to a-smupfi.py.

***Outputs***

a-smupfi provides four files for output.

1) XXXX-YYYY.txt => shows a graphical representation of where the mutation patterns occur
2) XXXX-YYYY_Summary.txt => shows only shared mutations between unique sequences

(if -e is given):
3) XXXX-YYYY_EasyOutputShared.txt => a file made for easy parsing, same as Summary.txt
4) XXXX-YYYY_EasyOutputUnique.txt => a file made for easy parsing, but selects all mutations
                                     on any sequence

================================choplqb.py================================

choplqb.py is a simple cleaning tool that uses average phred-score to proccess raw fasta
sequences. It requires two files, a .fasta or .fna file and a .qual file. This is usually
used with Next Generation Sequencing (NGS) data where the accuracy and trusthworthiness
of individual reads may decrease at the near the end of the read.

================================circosConverter.py================================

circosConverter.py is a format converter tool that simply changes easyoutput files from
a-smupfy.py to a format that Circos can understand. This tool is used for creating the
histogram around the circos plots.

================================indelRemover.py================================

A tool to remove indels from fasta sequences and replace it with reference sequence
characters and merges fasta sequences that end up being identical.

In detail, due to 454 technology having a tendency to have high errors in the form of 
'-' characters, also known as 'indels', around homopolymer regions, this tool specifically
targets indels around those regions and substitutes the indels with the character at 
the same position from the reference sequence.

In this tool a homopolymer region is any residue appearing 3 or more times consecutively.
For example, 'GGG'. Hence if an indel is detected to be around a homopolymer region like
'G-G', '-GG' or 'GG-', the tool will proceed to clean that indel.

================================jaccardToCircos.py================================

jaccardToCircos.py is a simple format converter tool that deals with jaccard easyoutput 
files and converts it to a format, which Circos can read and make arcs

================================javssim.py================================

javssim.py (Jaccard Value Statistical Simulator) is a simulator utilising easyOutputShared 
files from a-smupfi.py to conduct statistical test on paired mutation and returns a ranked list
of mutation that identifies which pairs are occurring more than random and which pairs are
occurring less than random.

================================omes.py================================

omes.py is an algorithm that uses the observed minus expected squared formula to calculate which
pair of of positions are co-varying.

================================qualfa2fq.pl/qualfa2fq_bt.pl================================

qualfa2fq.pl is a merger tool that intergrates a .fasta / .fna file and a .qual file into a
.fastq/.fq file.

qualfa2fq_bt.pl is a version doing the same process, but made for aligners (e.g. bowtie2) 
that dislike the 60 character limit in the .fastq /.fq files.

================================sequenceExpand.py================================

Uses fasta sequences with frequency information in the fasta tag as input and expands
according to a base number. Using ShoRAH output as an example:

>sequence1_0.3542
AGTCGTATTTAC

If base is given as -b 100, then sequence 1 will be expanded to 35 reads. If -b 1000,
then sequence 1 will be expanded as 354.

================================sif_converter.py================================

Takes in outputs from omes.py and converts the user defined pairs into .sif format which can be
taken in by Cytoscape to draw networks. The tool has the option -t which allows users to define
their own threshold for selecting data. 

================================snpExtractor.py================================
This script takes snps recorded in a LoFreq file and substitute the snp into reference sequence 
for translation and records amino acid mutations.    

*Note: This tool only takes in 3 columns of data. Any extra information will be ignored. 

If unsure if option -rp is needed, use:
cat -A FILE or cat -v FILE to check whether ^M character exists, if seen, use -rp.

For argument -fup, the tool will ignore frequencies greater than inputed frequency. This is to
ignore outputs that fall under the case such as a wrong consensus reference is used to align 
and results in high instances of mutation recorded.


================================snpExtractConverter.py================================

snpExtractConverter.py changes the format of the output of the snpExtractor.py into histogram
format that can be understood by Circos.

================================uniqvariant.py================================

uniqvariant.py takes in fasta files that have (potentially) repeated sequences and merges them into
a unique sequence and provides that sequence with a frequency of occurrence number. For example with
3 sequences:

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
 