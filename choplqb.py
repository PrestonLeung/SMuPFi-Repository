#!/usr/bin/python
from __future__ import print_function

def window_qual(qual_list, window_size):
    '''Return a new list of quality by sliding a window on original quality.

    If the window_size parameter is equal to zero or greater than the length (L)
    of original quality, the window_size is set to L.

    When the window slides along the list of original quality, the average value
    is assinged to the element of new list of quality in the same position to the
    first element of the window. The elements in the last window are assigned to
    the same value.
    '''
    if window_size == 1:
        return qual_list

    if window_size == 0 or len(qual_list) < window_size:
        window_size = len(qual_list)

    new_qual = qual_list[:]
    size = float(window_size)
    window_sum = sum(qual_list[0:window_size])
    new_qual[0] = window_sum/size
    n = 1
    last_element_index = n-1
    new_element_index = n+window_size-1
    while n <= len(qual_list) - window_size:
        window_sum += qual_list[new_element_index]-qual_list[last_element_index]
        new_qual[n] = window_sum/size
        last_element_index += 1
        new_element_index += 1
        n += 1
    while n < len(qual_list):
        new_qual[n] = new_qual[n-1]
        n += 1
    return new_qual

def find_chop_position(qual_list, qual_threshold):
    '''Return a tuple of chop position in the left end and right end.

    The elements which is less than qual_threshold in the left and right end of quality
    will be chopped.

    The returned position tuple (left, right) is zero-based. So use qual[left:right] can
    get a copy of remained quality.
    '''
    # looking for left side
    left = 0
    while left < len(qual_list):
        if qual_list[left] >= qual_threshold:
            break
        left += 1
    if left == len(qual_list):
        return left, left
    # looking for right side
    right = len(qual_list)-1
    while right >= 0:
        if qual_list[right] >= qual_threshold:
            break
        right -= 1
    return left, right+1

if __name__ == '__main__':
    prog_function = """
The average quality is calculated when a window with size WINDOW_SIZE slides
along each read sequence. Then chop low quality (lower than QUALITY_THRESHOLD)
bases at the left and right ends of 454 reads.

Special case: if the WINDOW_SIZE is set to zero, it removes the read whose
average quality is lower than the threshold.

The result will be output in new files whose names include WINDOW_SIZE and QUALITY_THRESHOLD
information. For example, if the sequence file and quality file are named filename.fna and
filename.qual, respectively, and suppose WINDOW_SIZE is 8 and QUALITY_THRESHOLD is 20, then
the result will be saved to filenameT20W8.fna and filenameT20W8.qual, respectively.
"""
    # Parse the command line arguments
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastafilename', help='filename of sequence in FASTA format', required=True)
    parser.add_argument('-q', '--qualityfilename', help='filename of quality corresponding to the sequence', required=True)
    parser.add_argument('-w', '--window-size', type=int, help='size of sliding window', required=True)
    parser.add_argument('-t', '--quality-threshold', type=int, help='threshold for removing the low quality bases', required=True)
    #parser.add_argument('-n', '--number-N', type=int, help='number of base N', required=True)
    parser.add_argument('-d', '--debug', action='store_true', help='for each read, output if it is chopped or not')
    parser.add_argument('-D', '--output-dir', default='.', help='the directory within which the result will keep, default is current working directory')
    if len(sys.argv) <= 1:
        print(prog_function)
        parser.print_help()
        sys.exit()
    args=parser.parse_args()

    # Set the output filename for the chopped sequence
    import os.path
    path, filename = os.path.split(args.fastafilename)
    filename_base, filename_ext = os.path.splitext(filename)
    output_readsfile = os.path.normpath(args.output_dir) + '/' + filename_base + 'T' + str(args.quality_threshold) + 'W' + str(args.window_size)
    #if args.number_N != sys.maxint:
    #    output_readsfile += str(args.number_N)
    output_readsfile += filename_ext
    output_reads_handler = open(output_readsfile, 'w')

    # Set the output filename for the chopped quality
    path, filename = os.path.split(args.qualityfilename)
    filename_base, filename_ext = os.path.splitext(filename)
    output_qualityfile = os.path.normpath(args.output_dir) + '/' + filename_base + 'T' + str(args.quality_threshold) + 'W' + str(args.window_size)
    #if args.number_N != sys.maxint:
    #    output_qualityfile += str(args.number_N)
    output_qualityfile += filename_ext
    output_quality_handler = open(output_qualityfile, 'w')

    threshold = float(args.quality_threshold)
    #number_N = args.number_N

    # Chop the sequence and quality
    from Bio.SeqIO.QualityIO import PairedFastaQualIterator
    from Bio import SeqIO
    import re
    re_pattern = re.compile(r'(length=)\d+(.*)')
    for seq_qual_record in PairedFastaQualIterator(open(args.fastafilename), open(args.qualityfilename)):
        qual_list = seq_qual_record.letter_annotations['phred_quality']
        lhs, rhs = find_chop_position(window_qual(qual_list, args.window_size), args.quality_threshold)
        if args.debug:
            print('{0},{1},{2}'.format(seq_qual_record.id, lhs, rhs), file=sys.stderr)

        if lhs == len(qual_list):
            if args.debug:
                print('Sequence ' + seq_qual_record.id + ' is abandoned', file=sys.stderr)
            continue
        elif lhs != 0 or rhs != len(qual_list):
            new_qual = qual_list[lhs:rhs]
            description = seq_qual_record.description
            new_desc = re_pattern.split(description)
            seq_qual_record.description = new_desc[0] + new_desc[1] + str(len(new_qual)) + new_desc[2]
            seq_qual_record.letter_annotations = {}
            seq_qual_record.seq = seq_qual_record.seq[lhs:rhs]
            seq_qual_record.letter_annotations['phred_quality'] = new_qual
            if args.debug:
                print('Sequence ' + seq_qual_record.id + ' is chopped', file=sys.stderr)
        else:
            if args.debug:
                print('Sequence ' + seq_qual_record.id + ' is kept', file=sys.stderr)
        SeqIO.write(seq_qual_record, output_reads_handler, 'fasta')
        SeqIO.write(seq_qual_record, output_quality_handler, 'qual')

