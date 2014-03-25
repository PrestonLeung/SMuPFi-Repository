from __future__ import print_function

if __name__ == '__main__':
    prog_function = '''

Merge the same sequences and sum up the frequencies.

Note:
    - The frequence in the FASTA name should be seperated from others by underscore (_). 
    - If you provide the PREFIX, an underscore (_) will be added automatically between 
      PREFIX and the original FASTA name.

    For example if FASTA tag is:
    >seq1_0.5
    AGT
    
    the FREQ_INDEX is 2. If the PREFIX is 'run1', then in the result, the FASTA name becomes
    
    >run1_seq1_0.5
    AGT

    - If freq-index (option -i) is not provided, every sequence will be counted as 1 and the 
      percentage of occurrence of one particular sequenced will be the sum of identical 
      sequence divided by total number of sequences in the file.

'''
    import sys
    import argparse
    totalSeq = 0
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastafilename', help='filename of sequence in FASTA format', required=True)
    parser.add_argument('-i', '--freq-index', type=int, help='position of frequency in the FASTA name')
    parser.add_argument('-o', '--output-file', help='output filename of result')
    parser.add_argument('-p', '--prefix', help='prefix added to the FASTA name')
    if len(sys.argv) <= 1:
        print(prog_function)
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    output_handler = None
    if not args.output_file:
        output_handler = sys.stdout
    else:
        output_handler = open(args.output_file, 'w')

    from Bio import SeqIO
    all_seq = {}
    for seq_record in SeqIO.parse(args.fastafilename, 'fasta'):
        
        if(args.freq_index): #for sequences that look like shorah
            desc = seq_record.id.split('_')
            freq = float(desc[args.freq_index-1])
        else: # for sequences that look like sga
            freq = 1.0
            totalSeq += 1.0
        seq = str(seq_record.seq)
        if seq in all_seq:
            all_seq[seq] += freq
        else:
            all_seq[seq] = freq
            

    
    all_seq = sorted(all_seq.items(), key=lambda x: x[1], reverse=True) #changes a dictionary to a list of sorted tuples (sort by key)    
    n = 0    
    
    from Bio.SeqRecord import SeqRecord
    from Bio import Seq
    
    for seq in all_seq:
        n += 1        
        
        seq_record = SeqRecord(Seq.Seq(seq[0]))
        if(args.prefix == None):
            args.prefix = ''
        if(args.freq_index):
            seq_record.id = args.prefix + 'HAP' + str(n) + '_' + str(seq[1])
            seq_record.description = ''
        else:
            seq_record.id = args.prefix + 'HAP' + str(n) + '_' + str(round(seq[1] / totalSeq, 3))
            seq_record.description = ''
        SeqIO.write(seq_record, output_handler, 'fasta')
