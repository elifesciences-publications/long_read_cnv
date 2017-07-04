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
import tempfile
import pysam

from long_read_pipeline.utils import fasta
from long_read_pipeline.utils import waterman
from long_read_pipeline.utils import gt_likelihoods



class CNVrow(object):

    def __init__(self, chrom1, breakStart1, breakEnd1, chrom2, breakStart2, breakEnd2, ID, score, strand1, strand2, queryStart1, queryEnd1, queryStart2, queryEnd2, minNonOverlap, queryLength, qualScores, variant_type, unaccounted_for_sequence, event_length, event_id, bam_file, r1, r2, break_point_folder):
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
        self._bam_file = pysam.AlignmentFile(bam_file, "rb")
        self._fasta_one = r1
        self._fasta_two = r2
        self._break_point_folder = break_point_folder 

    @property
    def strand1(self):
        return self._strand1
    @property
    def strand2(self):
        return self._strand2
    
    @property
    def svlen(self):
        return self._event_length
    @property 
    def unaccounted_for_sequence(self):
        return self._unaccounted_for_sequence
    @property
    def breakstart1(self):
        return self._breakStart1
    @property
    def breakstart2(self):
        return self._breakStart2
    @property 
    def breakend1(self):
        return self._breakEnd1
    @property 
    def breakend2(self):
        return self._breakEnd2

    @property
    def read_one(self):
        return self._r1
    @property
    def read_two(self):
        return self._r2

    @property
    def chrom(self):
        return self._chrom1

    @property 
    def variant_type(self):
        return self._variant_type

    @property
    def breakStart1(self):
        return self._breakStart1
    @property 
    def breakStart2(self):
        return self._breakStart2
    @property
    def breakEnd1(self):
        return self._breakEnd1
    @property 
    def breakEnd2(self):
        return self._breakEnd2
    @property
    def event_id(self):
        return self._event_id
    @property
    def bam_file(self):
        return self._bam_file
    

    def __str__(self):
        return '\t'.join(map(str, [self.chrom, self.breakStart1, self.breakEnd1, self.chrom, self.breakStart2, self.breakEnd2, self.ID, self.score, self.strand1, self.strand2, self.queryStart1, self.queryEnd1, self.queryStart2, self.queryEnd2  , self.minNonOverlap, self.queryLength, self.qualScores, self.variant_type, self.unaccounted_for_sequence, self.event_length,self.chrom1+":"+self.breakEnd1 +"-" +self.breakStart2]))

    def _prob_mapq(self, mapq):
        """
            Probability that a read maps
        """
        return 1 - 10**(-mapq/10)

    def _filt_reads(self, mapping_quality, mapped_bases, clip_size):
        """
            Filter reads.

            @author James Boocock

            @date 10 Feb 2017
        """
        
        return None

    def _get_left_bp(self, slop, mapping_quality, duplication=False, clip_size=5, add_slash=False, allow_slop=False, inversion=False): 
        """
            Get the left break point
        """
        if duplication:
            reads_tmp = self.bam_file.fetch(self.chrom, self.breakStart2 - slop, self.breakStart2 + slop)
        else:
            reads_tmp = self.bam_file.fetch(self.chrom, self.breakStart1 - slop, self.breakStart1 + slop)
        supporting = 0
        not_supporting = 0
        read_one_list = []
        read_two_list = []
        read_one_complement =[]
        read_two_complement =[]
        read_pos_r1 = []
        read_pos_r2 = []
        read1_pos_mapq=[]
        read2_pos_mapq=[]
        if duplication:
            breakStart1 = self._breakStart2
        else:
            breakStart1 = self._breakStart1
        total_bad = 0
        for reads in reads_tmp:
            if reads.mapping_quality <= mapping_quality or reads.is_duplicate:
                continue
            read_split = reads.cigartuples[len(reads.cigartuples)-1]
            if not (read_split[0] == 4 or read_split[0] == 5):
                if reads.reference_start > breakStart1:
                    continue
                if (reads.reference_start + reads.reference_length) > (breakStart1):
                    if duplication:
                        if reads.reference_start > (self._breakStart1 - 1):
                            continue
                    not_supporting += self._prob_mapq(reads.mapping_quality)
                continue
            type_clip = read_split[0]
            if reads.cigartuples[0][0] == 4 or reads.cigartuples[0][0] == 5:
                continue
            if read_split[1] < clip_size:
                continue
            if (reads.reference_start + reads.reference_length + read_split[1]) < breakStart1:
                continue
            #print(reads)
            #print(reads.reference_length)
            #if type_clip == 4:
            #    break_point = reads.reference_start  + (reads.query_length - read_split[1])
            #elif type_clip == 5:
            #    break_point = reads.reference_start  + (reads.query_length)
            # Exact breakpoint
            break_point = reads.reference_start + reads.reference_length
            max_slop_break = 10
            break_point_match = False
            if allow_slop:
                if abs(breakStart1 - break_point) < 10: 
                    break_point_match = True 
            else:
                if breakStart1 == break_point:
                    break_point_match = True
            if break_point_match: 
                if reads.is_read1:
                    if add_slash:
                        read_one_list.append(reads.query_name + "/1")
                    else:
                        read_one_list.append(reads.query_name)
                    read_one_complement.append(reads.is_reverse)
                    read_pos_r1.append(reads.reference_start) 
                    read1_pos_mapq.append(reads.mapping_quality)
                    #fasta.extract_reads(self._fasta_one, reads.query_name + "/1")
                else:
                    if add_slash:
                        read_two_list.append(reads.query_name + "/2")
                    else:
                        read_two_list.append(reads.query_name)
                    read_two_complement.append(reads.is_reverse)
                    read_pos_r2.append(reads.reference_start) 
                    read2_pos_mapq.append(reads.mapping_quality)
                    #fasta.extract_reads(self._fasta_two, reads.query_name +"/2")
                #  Check for matches
            elif break_point != breakStart1: 
                not_supporting += self._prob_mapq(reads.mapping_quality)
        temp_out = tempfile.NamedTemporaryFile(delete=False) 
        for read_one in read_one_list:
            temp_out.write(read_one +"\n")
        temp_out.close()
        tmp_fasta1 = fasta.extract_reads(self._fasta_one, temp_out.name) 
        os.remove(temp_out.name)
        temp_out = tempfile.NamedTemporaryFile(delete=False) 
        for read_two in read_two_list:
            temp_out.write(read_two + "\n")
        temp_out.close()
        tmp_fasta2 = fasta.extract_reads(self._fasta_two, temp_out.name) 
        os.remove(temp_out.name)
        fasta_check = fasta.FastaInput(tmp_fasta1, tmp_fasta2, read_one_complement, read_two_complement, read_pos_r1, read_pos_r2,read1_pos_mapq, read2_pos_mapq)
        with open(os.path.join(self._break_point_folder, self._event_id +  ".brk1.fasta")) as break_in:
            #with open(self._break_point_folder 
            break_in.readline()
            seq_break = break_in.readline().strip()
            for read_id, item in fasta_check.sequences.items():
                pos_in_brk1 = (item[1] - (breakStart1 -300))
                tmp_seq_brk = (seq_break[pos_in_brk1:])
                tmp_seq_read = item[0][:len(tmp_seq_brk)]
                tmp_seq_brk = tmp_seq_brk[:len(tmp_seq_read)]
                identity = waterman.water(tmp_seq_brk, tmp_seq_read)
                if identity > 0.9:
                    supporting += self._prob_mapq(item[2])
                else:
                    not_supporting += self._prob_mapq(item[2])
            #print(read_id)
        return supporting, not_supporting

    def _get_right_bp(self, slop, mapping_quality, duplication=False, clip_size = 5, add_slash=False, allow_slop=False):
        """
            Get the right breakpoint 
        """
        supporting = 0
        not_supporting = 0
        read_one_list = []
        read_two_list = []
        read_one_complement =[]
        read_two_complement =[]
        read_pos_r1 = []
        read_pos_r2 = []
        read1_pos_mapq=[]
        read2_pos_mapq=[]
        if duplication:
            breakStart2 = self._breakStart1
        else:
            breakStart2 = self._breakStart2
        if duplication:
            break_point_tmp = self._breakStart1
        else:
            break_point_tmp = self._breakEnd2
        reads_tmp = self.bam_file.fetch(self.chrom, breakStart2 - slop, breakStart2 + slop)
        for reads in reads_tmp:
            if reads.mapping_quality <= mapping_quality or reads.is_duplicate:
                continue
            read_split = reads.cigartuples[0]
            if not (read_split[0] == 4 or read_split[0] == 5):
                # Ensure end is not working
                if reads.reference_start > breakStart2:
                    continue
                if (reads.reference_start < breakStart2) and (reads.reference_start + reads.query_length) > (breakStart2):
                    if duplication:
                        continue
                    not_supporting +=1
                continue
            type_clip = read_split[0]
            if reads.cigartuples[len(reads.cigartuples) -1][0] == 4 or reads.cigartuples[len(reads.cigartuples)-1][0] == 5:
                continue
            if read_split[1] < clip_size:
                continue
            if type_clip == 4:
                #print(reads)
                #print(read_split)
                break_point = reads.reference_start 
            elif type_clip == 5:
                #print(reads)
                #print(read_split)
                break_point = reads.reference_start 
            break_point_match = False
            if allow_slop:
                if abs(break_point_tmp - break_point) < 10: 
                    break_point_match = True 
            else:
                if break_point_tmp == break_point:
                    break_point_match = True
            if break_point_match: 
                if reads.is_read1:
                    if add_slash:
                        read_one_list.append(reads.query_name + "/1")
                    else:
                        read_one_list.append(reads.query_name)
                    read_one_complement.append(reads.is_reverse)
                    read_pos_r1.append(reads.reference_start - read_split[1])
                    #print("HERE")
                    #print(reads.reference_start - read_split[1])
                    read1_pos_mapq.append(reads.mapping_quality)
                    #fasta.extract_reads(self._fasta_one, reads.query_name + "/1")
                else:
                    if add_slash:
                        read_two_list.append(reads.query_name + "/2")
                    else:
                        read_two_list.append(reads.query_name)
                    read_two_complement.append(reads.is_reverse)
                    read_pos_r2.append(reads.reference_start - read_split[1]) 
                    read2_pos_mapq.append(reads.mapping_quality)
                    #print(reads.reference_start - read_split[1])
                    #fasta.extract_reads(self._fasta_two, reads.query_name +"/2")
                #  Check for matches 
            elif break_point != self.breakStart1: 
                not_supporting +=1
        temp_out = tempfile.NamedTemporaryFile(delete=False) 
        for read_one in read_one_list:
            temp_out.write(read_one +"\n")
        temp_out.close()
        tmp_fasta1 = fasta.extract_reads(self._fasta_one, temp_out.name) 
        os.remove(temp_out.name)
        temp_out = tempfile.NamedTemporaryFile(delete=False) 
        for read_two in read_two_list:
            temp_out.write(read_two + "\n")
        temp_out.close()
        tmp_fasta2 = fasta.extract_reads(self._fasta_two, temp_out.name) 
        #os.remove(temp_out.name)
        fasta_check = fasta.FastaInput(tmp_fasta1, tmp_fasta2, read_one_complement, read_two_complement, read_pos_r1, read_pos_r2,read1_pos_mapq, read2_pos_mapq)
        with open(os.path.join(self._break_point_folder, self._event_id + ".brk2.fasta")) as break_in:
        #with open(self._break_point_folder 
            break_in.readline()
            seq_break = break_in.readline().strip()
            for read_id, item in fasta_check.sequences.items():
                pos_in_brk1 = (item[1] - (break_point_tmp-300))
                tmp_seq_brk = (seq_break[pos_in_brk1:])
                tmp_seq_read = item[0][:len(tmp_seq_brk)]
                tmp_seq_brk = tmp_seq_brk[:len(tmp_seq_read)]
                identity = waterman.water(tmp_seq_brk, tmp_seq_read)
                if identity > 0.90:
                    supporting +=1
                else:
                    not_supporting +=1
        return supporting, not_supporting

    def extract_windowed_reads_new_bam(self, bam_file, slop, mapping_quality=20): 

        tmp_bam_file = self._bam_file
        bam_fname = bam_file
        self._bam_file = pysam.AlignmentFile(bam_file, "rb")
        self.extract_windowed_reads(slop,tmp_bam_file=bam_fname, allow_slop=True)
        self._bam_file = tmp_bam_file 

    def extract_windowed_reads(self, slop, mapping_quality=20, add_slash=False, tmp_bam_file=None,allow_slop=False): 
        """
            Extract windowed reads for both R1 and R2.
        """
        if "duplication" in self._variant_type:
            duplication = True
        else:
            duplication = False
        # hack
        if "inversion" in self._variant_type:
            inversion=True
        else:
            inversion=False
        if tmp_bam_file is not None:
            sample = os.path.basename(tmp_bam_file).split(".")[0]
        else:
            sample = (os.path.dirname(self._break_point_folder).split("/")[len(os.path.dirname(self._break_point_folder).split("/"))-1])
        if sample == "RMx" or sample == "BYa" or sample == "RM" or sample == "BY":
            add_slash =True 
        left_support, left_no_support = self._get_left_bp(slop, mapping_quality, duplication=duplication, add_slash=add_slash,allow_slop=allow_slop) 
        right_support, right_no_support = self._get_right_bp(slop, mapping_quality, duplication=duplication, add_slash=add_slash,allow_slop=allow_slop)
        # Fisher's exact to ensure we are sampling both ends of the breakpoint
        (GQ, GT, SQ)= gt_likelihoods.genotype_likelihoods(left_support + right_support, left_no_support + right_no_support, duplication)
        self._GT = GT
        self._GQ = GQ
	self._SQ = SQ
        self._LS = left_support
        self._LNS = left_no_support
        self._RS = right_support 
        self._RNS = right_no_support 
        self._S = left_support + right_support 
        self._NS = left_no_support + right_no_support 
        if (self._S + self._NS) == 0:
            self._AB = 0
        else:
	    self._AB = float(left_support + right_support)/float(left_support + right_support + left_no_support + right_no_support)
        print("GQ = {0}, GT {1}, S {2}, NS {3}, EVENTID = {4}".format(GQ, GT, left_support + right_support, left_no_support + right_no_support, self._event_id))
	print("SQ = {0}, AB = {1}".format(SQ, self._AB))
        print("LS = {0}, LNS = {1}, RS = {2}, RNS = {3}".format(left_support, left_no_support, right_support, right_no_support))

    @property 
    def GT(self):
        return self._GT
    @property
    def GQ(self):
        return self._GQ
    @property
    def LS(self):
        return self._LS
    @property
    def RS(self):
        return self._RS
    @property
    def RNS(self):
        return self._RNS
    @property
    def LNS(self):
        return self._LNS
    @property
    def S(self):
        return self._S
    @property
    def NS(self):
        return self._NS

    @property 
    def AB(self):
        return self._AB
    @property
    def SQ(self):
        return self._SQ


class CNVFromVCF(object):

    def __init__(self, vcfs, bam_file, fasta_one, fasta_two, break_point_folder):
        """
            Create CNV object from CNV 
        """
        #def __init__(self, chrom1, breakStart1, breakEnd1, chrom2, breakStart2, breakEnd2, ID, score, strand1, strand2, queryStart1, queryEnd1, queryStart2, queryEnd2, minNonOverlap, queryLength, qualScores, variant_type, unaccounted_for_sequence, event_length, event_id, bam_file, r1, r2, break_point_folder):
        self.input_rows = []
        for ids, vcf_row in vcfs._vcf_dict.items():
            info_dict = (vcf_row.info_dict)
            breakstart1 = info_dict["BREAKSTART1"] 
            breakstart2 = info_dict["BREAKSTART2"] 
            breakend1 = info_dict["BREAKEND1"] 
            breakend2 = info_dict["BREAKEND2"] 
            strand1 = info_dict["STRANDS"][0]
            strand2 = info_dict["STRANDS"][1]
            ids = vcf_row.ids
            self.input_rows.append(CNVrow(vcf_row.chrom, breakstart1, breakend1,vcf_row.chrom,breakstart2,breakend2,".","",strand1,strand2,"","","","","","","",
                                   info_dict["SVTYPE"],50,0,ids,bam_file, fasta_one, fasta_two, break_point_folder))

    def extract_windowed_bam_reads(self,i, slop=200):
        self.input_rows[i].extract_windowed_reads(slop) 
    
    def extract_windowed_bam_reads_other_samples(self,i , bam_file, slop=200):
        """
            Extract other samples
        """
        self.input_rows[i].extract_windowed_reads_new_bam(bam_file,slop)

    def __len__(self):
        """
            Number of CNV events
        """
        return len(self.input_rows)
class CNVs(object):

    def __init__(self, cnv_input_file, bam_file,fasta_one,fasta_two, break_point_folder):
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
                      cnv_input_row_split[20],bam_file,fasta_one,fasta_two, break_point_folder))
#    def __init__(self, vcf):
#        """
#            Create CNV object from CNV 
#        """
#        #def __init__(self, chrom1, breakStart1, breakEnd1, chrom2, breakStart2, breakEnd2, ID, score, strand1, strand2, queryStart1, queryEnd1, queryStart2, queryEnd2, minNonOverlap, queryLength, qualScores, variant_type, unaccounted_for_sequence, event_length, event_id, bam_file, r1, r2, break_point_folder):
#        print(vcf) 
#


    def extract_windowed_bam_reads(self,i, slop=200):
        self.input_rows[i].extract_windowed_reads(slop) 

    def __len__(self):
        """
            Number of CNV events
        """
        return len(self.input_rows)
