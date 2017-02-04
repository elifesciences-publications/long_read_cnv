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

import os
import pysam
import argparse

def split_fasta(fasta_input, split_size):
    """
    Split fasta files

    """
    fasta_input = pysam.FastaFile(fasta_input)
    references = fasta_input.references 
    for reference in references:
        record = fasta_input.fetch(reference=reference)
        length = len(record)
        n_splits = length//split_size
        start = 0
        end = split_size 
        for i in range(n_splits):
            temp_seq = record[start:end]
            print(">" + reference + "_" + str(i))
            print(temp_seq)
            start = end
            end = start + split_size
        print(">" + reference + "_" + str(i+1))
        temp_seq = record[start:end]
        print(temp_seq)
def main():
    parser = argparse.ArgumentParser(description="Split fasta input into smaller file")
    parser.add_argument("fasta_file")
    parser.add_argument("-s", "-split_length", dest="split_length", default=2000)
    args = parser.parse_args()
    args.split_length = int(args.split_length)
    split_fasta(args.fasta_file, args.split_length)

if __name__ =="__main__":
    main()

