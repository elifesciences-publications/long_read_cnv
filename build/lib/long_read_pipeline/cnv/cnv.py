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

from long_read_pipeline.config import *
import os
import subprocess

def call_cnvs(input_file, mapping_quality, bam_file, output_directory, fasta_reference):
    """
       Call CNVS questions 
    """
    cnv_output_script = os.path.join(output_directory, input_file.samples_name) 
    cnv_analysis_script = "samtools view -q {0} -S {1} | " + SPLITREADBEDTOPE + " -i stdin | grep -v telo | " +SPLITTERTOBREAKPOINT + " -i stdin -s 0 -r {2} -o {3}"
    cnv_analysis_script =  cnv_analysis_script.format(mapping_quality, bam_file, fasta_reference, cnv_output_script) 
    subprocess.check_call(cnv_analysis_script,shell=True) 
