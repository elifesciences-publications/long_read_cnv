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


def get_cigar_start(cigar_tuple):

    if cigar_tuple[0][0] == 4:
        return cigar_tuple[0][1]
    return 0 


def get_cigar_end(cigar_tuple):

    if cigar_tuple[len(cigar_tuple)-1][0] == 4:
        return -cigar_tuple[0][1]
    return 0

def extract_sequence_from_alignments(bam_file, reference_genome):

    ref_genome = pysam.FastaFile(reference_genome)
    bam_in = pysam.Samfile(bam_file)
    ref_sequences  = (ref_genome.references) 
    for ref_seq in ref_sequences:
        reads_in_ref = bam_in.fetch(reference=ref_seq)
        fastq_tmp = ref_genome.fetch(reference=ref_seq)
        for read in reads_in_ref:
            start = 0
            start += (get_cigar_start(read.cigartuples))
            end = start + read.query_length + get_cigar_end(read.cigartuples)
            
            # variant scan


            
            



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
