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
import argparse
import sys
import subprocess
from long_read_pipeline.utils import sequence_manipulation

class FastaInput(object):
    """
        FastaInput

        simple class representing FASTA format

        Only does two line fasta to begin with
    """
    def __init__(self,r1_input_string, r2_input_string, complement_r1, complement_r2, pos_r1, pos_r2, read1_pos_mapq, read2_pos_mapq):
        self.sequences = {}
        j = 0 
        for i, line in enumerate(r1_input_string.split("\n")):
            if i % 2 == 0:
                if line.strip() == "":
                    break
                get_id = line.strip().split(">")[1]
            else:
                seq = line.strip()
                if complement_r1[j]:
                    self.sequences[get_id] = (sequence_manipulation.rev_complement(seq), pos_r1[j], read1_pos_mapq[j])
                else:
                    self.sequences[get_id] = (seq, pos_r1[j], read1_pos_mapq[j]) 
                j += 1
        j = 0 
        for i, line in enumerate(r2_input_string.split("\n")):
            if i % 2 == 0:
                if line.strip() == "":
                    break
                get_id = line.strip().split(">")[1]
            else:
                seq = line.strip()
                if complement_r2[j]:
                    self.sequences[get_id] = (sequence_manipulation.rev_complement(seq), pos_r2[j],read2_pos_mapq[j]) 
                else:
                    self.sequences[get_id] = (seq, pos_r2[j],read2_pos_mapq[j])
                j+=1 


def index_fasta(input_file, temp_dir, skip=False):
    """
        Fastaq to fasta.

        Uses sed to create fastas from fastq (ignoring quality)

        Then index these files with cdyank 
    """
    zsed_conversion = "zcat {0} | sed -n '1~4s/^@/>/p;2~4p' > {1} "
    sed_conversion = "cat {0} | sed -n '1~4s/^@/>/p;2~4p' > {1} "
    zsed_conversion_a= "zcat {0} | sed -n '1~4s/^@/>/p;2~4p' >> {1} "
    sed_conversion_a = "cat {0} | sed -n '1~4s/^@/>/p;2~4p' >> {1} "
    cdbfasta = "cdbfasta {0}" 
    sed_conversion = "cat {0} > {1}"
    for i, pairs in enumerate(input_file.paired_end_list):
        p1 = pairs[0]
        p2 = pairs[1]
        if i > 0:
            if ".gz" in p1:
                #p1_out =  os.path.join(temp_dir, os.path.basename(p1.split('.gz')[0])) + ".fasta"
                subprocess.check_call(zsed_conversion_a.format(p1, p1_out), shell=True)
            else:
                #p1_out =  os.path.join(temp_dir, os.path.splitext(p1)[0]) + ".fasta"
                subprocess.check_call(sed_conversion_a.format(p1, p1_out),shell=True)
            if ".gz" in p2:
                #p2_out =  os.path.join(temp_dir, os.path.basename(p2.split('.gz')[0])) + ".fasta"
                subprocess.check_call(zsed_conversion_a.format(p2, p2_out),shell=True)
            else:
                #p2_out =  os.path.join(temp_dir, os.path.splitext(p2)[0]) + ".fasta"
                subprocess.check_call(sed_conversion_a.format(p2, p2_out),shell=True)
        else:
            if ".gz" in p1:
                p1_out =  os.path.join(temp_dir, os.path.basename(p1.split('.gz')[0])) + ".fasta"
                #print(p1_out)
                #subprocess.check_call(zsed_conversion.format(p1, p1_out), shell=True)
            else:
                p1_out =  os.path.join(temp_dir, os.path.splitext(p1)[0]) + ".fasta"
                #subprocess.check_call(sed_conversion.format(p1, p1_out),shell=True)
            if ".gz" in p2:
                p2_out =  os.path.join(temp_dir, os.path.basename(p2.split('.gz')[0])) + ".fasta"
                #subprocess.check_call(zsed_conversion.format(p2, p2_out),shell=True)
            else:
                p2_out =  os.path.join(temp_dir, os.path.splitext(p2)[0]) + ".fasta"
                #subprocess.check_call(sed_conversion.format(p2, p2_out),shell=True)
    #subprocess.check_call(cdbfasta.format(p1_out),shell=True)
    #subprocess.check_call(cdbfasta.format(p2_out),shell=True)
    input_file.set_paired_fastas(p1_out, p2_out)

def extract_reads(fasta_input, temp_read_list_file):
    """
        Extract reads. 
    """

    cdx_input = fasta_input+ ".cidx" 
    cdbyank = "cdbyank {0} < {1}"
    yank_command = cdbyank.format(cdx_input, temp_read_list_file)
    fasta = subprocess.check_output(yank_command, shell=True)
    return fasta 
