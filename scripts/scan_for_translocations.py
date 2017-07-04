#!/usr/bin/env python
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
import tempfile
import subprocess

def get_clipped_seq_length(bam_row):
    """
        Clipped sequence length
    """
    if bam_row.cigartuples[0] == 4 or bam_row.cigartuples[0] == 5:
        return bam_row.cigartuples[0][1]
    return 0

def get_inter_chromosomal_arrangements(reference,bam_file):

    # Sort bam by name
    read_dict = {}
    bam_py = pysam.AlignmentFile(bam_file)
    for read in bam_py:
        if read.mapping_quality > 30 and read.reference_length > 10000:
            try:
                read_dict[read.query_name].append(read)
            except KeyError:
                read_dict[read.query_name] = [read]
    # Now process each query
    # Chromosomes = hello.
    for key, values in read_dict.items():
        values_sorted = (sorted(values, key=lambda x: get_clipped_seq_length(x)))
        chromosomes = set()
        total_reference_length = 0 
        for value in values_sorted:
            chromosomes.add(value.reference_name)
    #        if len(chromosomes) > 1:
    #            print(value.query_name)
            total_reference_length += value.reference_length
        if len(chromosomes) > 1:
            chr_id = ("_".join(chromosomes))
            print(key)
            print("\t".join([chr_id, str(total_reference_length)]))
def main():
    parser = argparse.ArgumentParser(description="Interchromosomal events")
    parser.add_argument("-r","--reference", dest="reference")
    parser.add_argument("bam")
    args = parser.parse_args() 
    get_inter_chromosomal_arrangements(args.reference, args.bam)

if __name__ == "__main__":
    main()
