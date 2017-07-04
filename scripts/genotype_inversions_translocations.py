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


import argparse
from long_read_pipeline import input_file


import sys

def genotype_inversions(input_file, input_directory, output_directory, inversions):
    """
        Genotyping inversions using short-read data.. 
    """
    in_file = input_file.InputFile(args.input_file, args.input_directory)
    for sample in in_file:
        inversion_file = os.path.join(output_directory, sample.samples_name + ".inv")
        print(inversion_file)
        align.align_reads(sample,args.temp_dir, reference_file, align_to_ref=True,skip=True)
        print(sample.bam_file)
        # Read in file. 

        # Call each CNV using short-read data. 

    sys.exit(1)

def identify_inversions(args.input_file, args.input_directory, args.output_directory):
    """
        
    """

def genotype_translocations(transinput):
    """
        Genotyping translocations, and identify reciprocal ones.
    """ 

def main():
    parser = argparse.ArgumentParser(description="Genotypes inversions and translocations using short-read datasets") 
    parser.add_argument("-l","--translocations", dest="translocation_input", help="Translocation input")
    parser.add_argument("-i","--input-file", dest="input_file", help="Input file")
    parser.add_argument("-d", "--directory", dest="input_directory", help="Input directory", required=True)
    parser.add_argument("-o", "--output-directory", dest="output_directory", help="Output directory", default="out_dir")
    args = parser.parse_args()
    genotype_inversions(args.input_file, args.input_directory, args.output_directory)
    # Take each sample name . 
    genotype_translocations(args.input_file, args.input_directory, args.output_directory)
