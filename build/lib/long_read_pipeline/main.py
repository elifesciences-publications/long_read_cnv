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
#

import os
import argparse
import sys
import logging
import subprocess

from long_read_pipeline import input_file 
from long_read_pipeline.speedseq import align
from long_read_pipeline.utils import scaffold
from long_read_pipeline.pilon import pilon 
from long_read_pipeline.cnv import cnv

def setup_log(log_file):
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")

    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
    
    fileHandler = logging.FileHandler(log_file,mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)


def parse_args(): 
    parser=argparse.ArgumentParser(description="Causal estimation from summary statistics")

    subparsers = parser.add_subparsers(help="sub-command help")
    parser.add_argument("-l", "--log-file", dest="log_file", help="Log file",default="lr_cnv.log")
    pilon_align_parser= subparsers.add_parser("pilonalign", help="assembly_to_var")
    pilon_align_parser.add_argument("-i","--input-file", dest="input_file", help="Input file", required=True)
    pilon_align_parser.add_argument("-r","--reference-file", dest="reference_file", help="Reference file", required=True)
    pilon_align_parser.add_argument("-t","--temp-working-dir", dest="temp_dir", help="Reference file", default="tmp_dir") 
    pilon_align_parser.add_argument("-d", "--directory", dest="input_directory", help="Input directory", required=True)
    pilon_align_parser.add_argument("-o", "--output-directory", dest="output_directory", help="Output directory", default="out_dir") 
    pilon_align_parser.set_defaults(func=pilon_align_wrap) 
    # Causal betas
    cnv_parser= subparsers.add_parser("callcnvs", help="assembly_to_var")
    cnv_parser.add_argument("-i", "--input-file", dest="input_file", help="Input file", required=True)
    cnv_parser.add_argument("-t", "--temp-working-dir", dest="temp_dir", help="Reference file", default="tmp_dir") 
    cnv_parser.add_argument("-d", "--directory", dest="input_directory", help="Input directory", required=True)
    cnv_parser.add_argument("-o", "--output-directory", dest="output_directory", help="Output directory", default="out_dir") 
    cnv_parser.add_argument("-q", "--map-quality", dest="mapping_quality", help="Mapping Quality", default=50)
    cnv_parser.set_defaults(func=cnv_call_wrap)

    args = parser.parse_args()
    return args

def cnv_call_wrap(args):
    """
        CNV call
    """
    in_file = input_file.InputFile(args.input_file, args.input_directory)
    try:
        os.mkdir(args.temp_dir)
    except OSError:
        pass
    try:
        os.mkdir(args.output_directory)
    except OSError:
        pass
    for sample in in_file:
        bwa_input = os.path.join(args.input_directory, sample.samples_name +".align_to_ref.sorted.bam")
        fasta_reference = os.path.join(args.input_directory, sample.samples_name +".fasta")
        cnv.call_cnvs(sample, args.mapping_quality, bwa_input, args.output_directory, fasta_reference)

def pilon_align_wrap(args):

    in_file = input_file.InputFile(args.input_file, args.input_directory)
    try:
        os.mkdir(args.temp_dir)
    except OSError:
        pass
    try:
        os.mkdir(args.output_directory)
    except OSError:
        pass
    # Align reads only if they don't already exist.
    # Then,run pilon to fix the reads removing contigs that are bad.
    # If PE and SE files are present run pilon to fix the contigs otherwise just filter.
    for sample in in_file:
        temp_fasta = scaffold.scaffold_filter(sample.scaffolds, args.temp_dir, sample.min_contig_length)
        align.align_reads(sample,args.temp_dir, temp_fasta) 
        pilon.run_pilon_contigs(sample, args.temp_dir, temp_fasta)
        bwa_output = align.simple_bwa(sample, args.temp_dir, args.reference_file)
        os.rename(bwa_output, os.path.join(args.output_directory, os.path.basename(bwa_output)))
        os.rename(bwa_output +".bai", os.path.join(args.output_directory, os.path.basename(bwa_output) +".bai"))
        os.rename(sample.pilon_fasta,os.path.join(args.output_directory, os.path.basename(sample.pilon_fasta) ))
        samtools_faidx = " samtools faidx {0}".format(os.path.join(args.output_directory, os.path.basename(sample.pilon_fasta) ))
        subprocess.check_call(samtools_faidx, shell=True)
    # Then align the reads to reference genome.  
    #speedseq_align = speeqseq.speedseq_align(in_file)
    # 

def main():
    """
        Main 
    """
    logging.info("Welcome to long-read CNV calling (a tool for genotyping long reads)")
    args = parse_args()
    setup_log(args.log_file)
    args.func(args)

