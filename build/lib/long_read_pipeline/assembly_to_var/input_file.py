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


HEADER = ['SAMPLES_NAME', 'SCAFFOLDS', 'PE', 'SE', 'MIN_CONTIG_LENGTH']


class InputRow(object):
    """
        Input row.
    """
    def __init__(self, samples_name, scaffolds, paired_end, single_end=None,min_contig_length=10000):
        self.samples_name = samples_name
        self.scaffolds = scaffolds 
        self.paired_end = paired_end
        self.single_end = single_end
        self.min_contig_length = min_contig_length

    def __str__(self):
        return(self.samples_name +"\t" + self.scaffolds +"\t" + ",".join(self.paired_end) + "\t" + ",".join(self.paired_end)  + "\t"+ str(self.min_contig_length))

class InputFile(object):
    """
        Sample input file.  

    """
    def __init__(self, input_file, reference):
        """
            Initialise input files. 
        """
        self.input_rows = [] 
        idxs_in_header = []
        with open(input_file) as in_file:
            for i, line in enumerate(in_file):
                if i > 0:
                    row = line.strip().split()
                    paired_end_reads =  row[2].split(",")
                    single_end_reads =  row[3].split(",")
                    min_contig_length = int(row[4])
                    in_row = InputRow(row[0], row[1], paired_end_reads, single_end_reads, min_contig_length)
                    self.input_rows.append(in_row)
    def __iter__(self):
        return iter(self.input_rows) 
