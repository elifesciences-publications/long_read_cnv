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
import os
import sys
import logging

import pysam


CIGAR_TUPLES={

        }


"


class CNVrow(object):

    def __init__(self, chrom1, breakStart1, breakEnd1, chrom2, breakStart2, breakEnd2, ID, score, strand1, strand2, queryStart1, queryEnd1, queryStart2, queryEnd2, minNonOverlap, queryLength, qualScores, variant_type, unaccounted_for_sequence, event_length, event_id, bam_file):
        self._chrom1 = chrom1
        self._breakStart1 = int(breakStart1)
        self._breakEnd1 = int(breakEnd1) 
        self._chrom2 = chrom2
        self._breakStart2 = int(breakStart2)
        self._breakEnd2 = int(breakEnd2)
        self._ID = ID
        self._score = score
        self._strand1 = strand1
        self._strand2 = strand2
        self._queryStart1 =queryStart1
        self._queryEnd1 = queryEnd1
        self._queryStart2 = queryStart2 
        self._queryEnd2 = queryEnd2 
        self._minNonOverlap = minNonOverlap 
        self._queryLength = queryLength 
        self._qualScores = qualScores 
        self._variant_type = variant_type 
        self._unaccounted_for_sequence = unaccounted_for_sequence 
        self._event_length = event_length
        self._event_id = event_id
        print(bam_file)
        self._bam_file = pysam.AlignmentFile(bam_file, "rb")

    @property
    def chrom(self):
        return self._chrom1

    @property
    def breakStart1(self):
        return self._breakStart1
    @property 
    def breakStart2(self):
        return self._breakStart2
    @property
    def bam_file(self):
        return self._bam_file

    def __str__(self):
        return '\t'.join(map(str, [self.chrom1, self.breakStart1, self.breakEnd1, self.chrom2, self.breakStart2, self.breakEnd2, self.ID, self.score, self.strand1, self.strand2, self.queryStart1, self.queryEnd1, self.queryStart2, self.queryEnd2  , self.minNonOverlap, self.queryLength, self.qualScores, self.variant_type, self.unaccounted_for_sequence, self.event_length,self.chrom1+":"+self.breakEnd1 +"-" +self.breakStart2]))

    def extract_windowed_reads(self, slop): 
        """
            Extract windowed reads for both R1 and R2.
        """
        reads_tmp = self.bam_file.fetch(self.chrom, self.breakStart1 - slop, self.breakStart1 + slop)
        for reads in reads_tmp:
            
            print(reads)
            print(reads.cigartuples)
        return reads_tmp

class CNVs(object):

    def __init__(self, cnv_input_file, bam_file):
        """
            Initialise event IDs.
        """
        self.input_rows = []
        with open(cnv_input_file) as cnv_in:
            for i, cnv_input_row in enumerate(cnv_in):
               if i == 0:
                   continue
               cnv_input_row_split = cnv_input_row.split()
               self.input_rows.append(CNVrow(cnv_input_row_split[0],cnv_input_row_split[1],cnv_input_row_split[2],cnv_input_row_split[3],
                      cnv_input_row_split[4],cnv_input_row_split[5],cnv_input_row_split[6],cnv_input_row_split[7],
                      cnv_input_row_split[8],cnv_input_row_split[9],cnv_input_row_split[10],cnv_input_row_split[11],
                      cnv_input_row_split[12],cnv_input_row_split[13],cnv_input_row_split[14],cnv_input_row_split[15],
                      cnv_input_row_split[16],cnv_input_row_split[17],cnv_input_row_split[18],cnv_input_row_split[19],
                      cnv_input_row_split[20],bam_file))

    def extract_windowed_bam_reads(self,i, slop=1000):
        temp_reads = self.input_rows[i].extract_windowed_reads(slop) 

    def __len__(self):
        """
            Number of CNV events
        """
        return len(self.input_rows)
