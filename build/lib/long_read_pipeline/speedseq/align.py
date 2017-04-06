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

import logging
import subprocess
import os
from long_read_pipeline.config import *

def simple_bwa(input_file, temp_dir, reference_file):
    """
        Simple BWA
    """
    bwa_output = os.path.join(temp_dir, input_file.samples_name + ".align_to_ref.bam")
    bwa_command = " bwa mem {0} {1} | samtools view -bS - > {2}".format(reference_file, input_file.pilon_fasta, bwa_output)
    subprocess.check_call(bwa_command, shell=True)
    bwa_sorted_output= os.path.join(temp_dir, input_file.samples_name + ".align_to_ref.sorted.bam")
    bwa_command = " samtools sort -o {0} {1}".format(bwa_sorted_output, bwa_output)
    #bwa_sorted_output += ".bam"
    subprocess.check_call(bwa_command, shell=True)
    bwa_index  = "samtools index {0}".format(bwa_sorted_output)
    subprocess.check_call(bwa_index, shell=True)
    return bwa_sorted_output

def align_reads(input_file, temp_dir, reference_genome, skip=False):
    """
        Align reads using speedseq
    """
    logging.info("Aligning reads using speedseq")
    # Index reference

    if not os.path.isfile(reference_genome + ".pac"):
        bwa_command = BWA + " index " + reference_genome 
        logging.info("Indexing reference genome {0}".format(reference_genome))
        if not skip:
            subprocess.check_call(bwa_command,shell=True)
    # Run speedseq align
    for pairs in input_file.paired_end_list:
        p1 = pairs[0]
        p2 = pairs[1]
        read_group = "@RG\\tID:{0}\\tSM:{0}\\tLB:lib".format(input_file.samples_name)
        logging.info("Running speedseq align on {0}".format(input_file.samples_name))
        speed_seq_align = SPEEDSEQ_EXEC + " align " + " -o " + os.path.join(temp_dir, os.path.basename(p1))  + " -R \"{0}\"".format(read_group)+ " " +reference_genome + " " + p1 + " " + p2 
        if not skip:
            subprocess.check_call(speed_seq_align,shell=True)
        bam_file = os.path.join(temp_dir, os.path.basename(p1)) +".bam"
        input_file.set_bam_file(bam_file)
