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
import logging
import os

def scaffold_filter(scaffolds, temp_dir, min_contig_length):
    """
       Filter the scaffolds. 
    """
    
    assembly = pysam.Fastafile(scaffolds)
    idx_keep = [i for i, lengths in enumerate(assembly.lengths) if lengths > int(min_contig_length)] 
    references_to_keep = [reference  for i, reference in enumerate(assembly.references) if i in idx_keep]
    temp_fname = os.path.join(temp_dir,os.path.basename(scaffolds)) 
    with open(temp_fname,"w") as temp_scaffolds:
        for references in references_to_keep:
            seq = assembly.fetch(region=references)
            temp_scaffolds.write(">" +references + "\n" + seq + "\n")
    return temp_fname
    
