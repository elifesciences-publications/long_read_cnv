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
import logging
import subprocess
import os

def run_pilon_contigs(input_file, temp_dir, temp_fasta=None): 
    """
        java -jar -Xmx 32G pilon.jar  
    """
    logging.info("Running PILON to fix contigs")
    if temp_fasta is not None:
        pilon_jar = " java -Xmx4G -jar {0} --genome {1} --frags {2} --outdir {3} --output {4}".format(PILON_PATH, temp_fasta, input_file.bam_file, temp_dir, input_file.samples_name) 
    else:
        pilon_jar = " java -Xmx4G -jar {0} --genome {1} --frags {2} --outdir {3} --output {4}".format(PILON_PATH, temp_fasta, input_file.bam_file, temp_dir, input_file.samples_name) 
    output_file = os.path.join(temp_dir, input_file.samples_name + ".fasta")
    input_file.set_pilon_fasta(output_file)
    subprocess.check_call(pilon_jar, shell=True)

