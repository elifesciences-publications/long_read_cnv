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
import argparse
import operator
import os

CIGAR_TUPLES_DICT ={ 0: "M",
                     1: "I",
                     2: "D",
        }

VCF_HEADER="""##fileformat=VCFv4.1
##reference=file:///data/rr/Parents/Reference/sacCer3.fasta
##source=SelectVariants
##fileDate=20150512
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##contig=<ID=chrI,length=230218,assembly="sacCer3">
##contig=<ID=chrII,length=813184,assembly="sacCer3">
##contig=<ID=chrIII,length=316620,assembly="sacCer3">
##contig=<ID=chrIV,length=1531933,assembly="sacCer3">
##contig=<ID=chrV,length=576874,assembly="sacCer3">
##contig=<ID=chrVI,length=270161,assembly="sacCer3">
##contig=<ID=chrVII,length=1090940,assembly="sacCer3">
##contig=<ID=chrVIII,length=562643,assembly="sacCer3">
##contig=<ID=chrIX,length=439888,assembly="sacCer3">
##contig=<ID=chrX,length=745751,assembly="sacCer3">
##contig=<ID=chrXI,length=666816,assembly="sacCer3">
##contig=<ID=chrXII,length=1078177,assembly="sacCer3">
##contig=<ID=chrXIII,length=924431,assembly="sacCer3">
##contig=<ID=chrXIV,length=784333,assembly="sacCer3">
##contig=<ID=chrXV,length=1091291,assembly="sacCer3">
##contig=<ID=chrXVI,length=948066,assembly="sacCer3">
##contig=<ID=chrM,length=85779,assembly="sacCer3">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""


def get_cigar_start(cigar_tuple):

    if cigar_tuple[0][0] == 4:
        return cigar_tuple[0][1]
    return 0 


def get_cigar_end(cigar_tuple):

    if cigar_tuple[len(cigar_tuple)-1][0] == 4:
        return -cigar_tuple[len(cigar_tuple) - 1][1]
    return 0

def cig_start_clip(cigar_tuple):
    """     

    """
    if cigar_tuple[0][0] == 4 or cigar_tuple[0][1] == 5:
        return True
    return False


class Variants(object):

    def __init__(self, chrom, pos, ref,alt, variant_id, variant_type):
        self._chrom = chrom
        self._pos = pos
        self._ref = ref
        self._alt = alt
        self._read_ids = []
        self._read_count = 0
        self._not_supporting_count = 0
        self._variant_id = variant_id
        self._variant_type = variant_type
        if self.alt != self.ref:
            self._alts = [self.alt]
        else:
            self._alts = []

    @property
    def variant_type(self):
        return self._variant_type
    @property
    def variant_id(self):
        return self._variant_id
    @property
    def chrom(self):
        return self._chrom
    @property
    def pos(self):
        return self._pos
    @property
    def ref(self):
        return self._ref
    @property
    def alt(self):
        return self._alt
    @property
    def read_ids(self):
        return self._read_ids
    
    @property 
    def read_count(self):
        return self._read_count

    def set_alt(self, alt):
        if self.alt == self.ref:
            self._alt = alt
            self._alts.append(alt)
        elif self.alt != alt:
            self._alts.append(alt)
    
    def add_supporting(self, read_id):
        """
            Add supporting reads.
        """
        self._read_count +=1 
        self._read_ids.append(read_id)

    def add_not_supporting(self, read_id):
        """
            Add not supporting reads.
        """
        self._not_supporting_count +=1
        self._read_ids.append(read_id)

    def _get_gt_string(self):
        return self.alt

    def __str__(self):
        return self.chrom +"\t" + str(self.pos) + "\t" + self.variant_id + "\t" + self.ref + "\t" + ','.join(self._alts) + "\t.\t.\tVARIANT_TYPE="+self.variant_type +";READID=" + ",".join(self.read_ids) + "\tGT:DP:AD\t"+ "1/1:"+str(self.read_count +self._not_supporting_count) + ":" + str(self._not_supporting_count) + "," + str(self.read_count)

variants = {}

def _scan_read_for_variants(fastq_subset, read_subset, cigartuples, start_reference, chrom, read_id):
    i = 0
    read_index = 0
    if cig_start_clip(cigartuples): 
        cigar_tuple_index = 1
    else:
        cigar_tuple_index = 0
    next_index = 0
    stop = 0
    matches = False
    insertion = False
    deletion = False 
    global variants 
    while i < len(fastq_subset):
        if next_index == (i) or insertion: 
            matches = False
            insertion = False
            deletion = False
            #print(read_id)
            #print(cigartuples)
            #print(cigar_tuple_index)
            if cigartuples[cigar_tuple_index][0] == 0: 
                next_index += cigartuples[cigar_tuple_index][1] 
                matches = True
            elif cigartuples[cigar_tuple_index][0] == 1:
                # Insertion
                ins_length = cigartuples[cigar_tuple_index][1] 
                insertion = True
            elif cigartuples[cigar_tuple_index][0] == 2:
                # Deletion
                del_length = cigartuples[cigar_tuple_index][1] 
                deletion = True  
                next_index += del_length
            elif cigartuples[cigar_tuple_index][0] == 4 or cigartuples[cigar_tuple_index][0] == 5:
                break
            cigar_tuple_index += 1
            #print(len(cigartuples))
            if cigar_tuple_index >= (len(cigartuples)):
                break
        pos = start_reference + i + 1
        if matches:
            # Make sure that position doesn't already have a read!!!!! Cause it might have a non-supporting alignment.
            ref = fastq_subset[i]
            alt = read_subset[read_index]
            ids = chrom + ":" + str(pos) + "_" + ref 
            if fastq_subset[i] != read_subset[read_index]:
                try:
                    variants[ids].add_supporting(read_id)
                    variants[ids].set_alt(alt)
                except KeyError:
                    variants[ids] = Variants(chrom, pos, ref, alt, ids,"SNP") 
                    variants[ids].add_supporting(read_id)
            else:
                try:
                    variants[ids].add_not_supporting(read_id)
                except KeyError:
                    variants[ids] = Variants(chrom, pos, ref, alt, ids,"SNP") 
                    variants[ids].add_not_supporting(read_id)

            read_index += 1
            i += 1
        elif insertion:
            ref = (fastq_subset[i-1]) 
            alt = read_subset[(read_index-1):(read_index+ ins_length )]
            ids = chrom + ":" + str(pos) + "_" + ref + "/" + alt
            read_index += (ins_length)
            try:
                variants[ids].add_supporting(read_id) 
            except KeyError:
                variants[ids] = Variants(chrom, pos, ref, alt,  ids,"INDEL") 
                variants[ids].add_supporting(read_id)
        elif deletion:
            ref =  fastq_subset[(i-1):i+del_length]
            alt = read_subset[read_index+1]
            i += del_length 
            ids = chrom + ":" + str(pos) + "_" + ref + "/" + alt
            try:
                variants[ids].add_supporting(read_id) 
            except KeyError:
                variants[ids] = Variants(chrom, pos, ref, alt, ids,"INDEL") 
                variants[ids].add_supporting(read_id)

def extract_sequence_from_alignments(bam_file, reference_genome):

    ref_genome = pysam.FastaFile(reference_genome)
    bam_in = pysam.Samfile(bam_file)
    ref_sequences  = (ref_genome.references) 
    global variants
    sample_name = (os.path.basename(bam_file).split(".")[0])
    print(VCF_HEADER)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+ sample_name)
    for ref_seq in ref_sequences:
        reads_in_ref = bam_in.fetch(reference=ref_seq)
        fastq_tmp = ref_genome.fetch(reference=ref_seq)
        for read in reads_in_ref:
            # TODO change it so it isn't fixed at mapping quality 30 ....
            if read.mapping_quality > 30:
                start_read = get_cigar_start(read.cigartuples)
                end_read = read.query_length + get_cigar_end(read.cigartuples)
                start_reference = read.reference_start 
                end_reference = start_reference + read.reference_length
                fastq_subset = (fastq_tmp[start_reference:end_reference])
                read_subset = read.query_sequence[start_read:end_read]
                _scan_read_for_variants(fastq_subset, read_subset, read.cigartuples, start_reference, ref_seq, read.query_name) 
        tuple_list = (sorted(variants.items(), key=lambda x: x[1].pos))
        for ids, variant in tuple_list:
            if variant.read_count > 0:
                print(str(variant))
            #continue
        variants = {}
            # variant scan, scan the long reads for variants


            
            



def main():
    parser = argparse.ArgumentParser(description="Extract indels and SNPs")
    parser.add_argument("-r", "--reference", dest="reference_genome")
    parser.add_argument("-b", "--bam-file", dest="bam_file")
    args = parser.parse_args()
    bam_file = args.bam_file
    reference = args.reference_genome
    extract_sequence_from_alignments(bam_file , reference)

if __name__ == "__main__":
    main()
