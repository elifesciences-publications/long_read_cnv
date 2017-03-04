# Copyright (c) 2017 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import pysam
import argparse


CIGAR_TUPLES_DICT ={ 0: "M",
                     1: "I",
                     2: "D",
        }

def get_cigar_start(cigar_tuple):

    if cigar_tuple[0][0] == 4:
        return cigar_tuple[0][1]
    return 0 


def get_cigar_end(cigar_tuple):

    if cigar_tuple[len(cigar_tuple)-1][0] == 4:
        return -cigar_tuple[len(cigar_tuple) - 1][1]
    return 0

def cig_start_clip(cigar_tuple):
    """     

    """
    if cigar_tuple[0][0] == 4 or cigar_tuple[0][1] == 5:
        return True
    return False


class Variant(object):

    def __init__(self, chrom, pos, ref,alt, read_id):
        return None

    def add_supporting(self, chrom, pos, ref, alt, read_id):
        return None


def _scan_read_for_variants(fastq_subset, read_subset, cigartuples, start_reference):
    i = 0
    read_index = 0
    cigar_tuple_index = 1
    next_index = 0
    stop = 0
    matches = False
    insertion = False
    deletion = False 
    while i < len(fastq_subset):
        if next_index == (i) or insertion: 
            matches = False
            insertion = False
            deletion = False 
            if cigartuples[cigar_tuple_index][0] == 0: 
                next_index += cigartuples[cigar_tuple_index][1] 
                matches = True
            elif cigartuples[cigar_tuple_index][0] == 1:
                # Insertion
                ins_length = cigartuples[cigar_tuple_index][1] 
                insertion = True
            elif cigartuples[cigar_tuple_index][0] == 2:
                # Deletion
                del_length = cigartuples[cigar_tuple_index][1] 
                deletion = True  
                next_index += del_length
            cigar_tuple_index += 1
        if matches:
            if fastq_subset[i] != read_subset[read_index]:
                print(start_reference + i +1)
                print(fastq_subset[i] + "->" + read_subset[read_index])
            read_index += 1
            i += 1
        elif insertion:
            print("INS")
            print(start_reference + i + 1)
            print(fastq_subset[i-1] + "->" + read_subset[(read_index-1):(read_index+ ins_length )])
            read_index += (ins_length)
            print(ins_length)
            print(read_subset[(read_index):(read_index+4)])
        elif deletion:
            print("DEL")
            print(start_reference +i + 1)
            print(fastq_subset[(i-1):i+del_length] + "->" + read_subset[read_index+1])
            i += del_length 
            print(read_subset[(read_index):(read_index+10)])
            print(fastq_subset[(i):(i+10)])


def extract_sequence_from_alignments(bam_file, reference_genome):

    ref_genome = pysam.FastaFile(reference_genome)
    bam_in = pysam.Samfile(bam_file)
    ref_sequences  = (ref_genome.references) 
    for ref_seq in ref_sequences:
        reads_in_ref = bam_in.fetch(reference=ref_seq)
        fastq_tmp = ref_genome.fetch(reference=ref_seq)
        for read in reads_in_ref:
            if read.mapping_quality > 30:
                start_read = get_cigar_start(read.cigartuples)
                end_read = read.query_length + get_cigar_end(read.cigartuples)
                start_reference = read.reference_start 
                end_reference = start_reference + read.reference_length
                fastq_subset = (fastq_tmp[start_reference:end_reference])
                read_subset = read.query_sequence[start_read:end_read]
                _scan_read_for_variants(fastq_subset, read_subset, read.cigartuples, start_reference) 
            # variant scan, scan the long reads for variants


            
            



def main():
    parser = argparse.ArgumentParser(description="Extract indels and SNPs")
    parser.add_argument("-r", "--reference", dest="reference_genome")
    parser.add_argument("-b", "--bam-file", dest="bam_file")
    args = parser.parse_args()
    bam_file = args.bam_file
    reference = args.reference_genome
    extract_sequence_from_alignments(bam_file , reference)

if __name__ == "__main__":
    main()
