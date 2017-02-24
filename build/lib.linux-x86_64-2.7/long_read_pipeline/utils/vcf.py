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

import collections

default_vcf_header="""##fileformat=VCFv4.2
##source=LRCNV
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=UNACCOUNTED_SEQUENCE,Number=.,Type=Integer,Description="Sequence that cannot be accounted for in the comparison of the alignments, +ve novel sequence, -ve deleted sequence overall">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=BREAKSTART1,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=BREAKEND1,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=BREAKSTART2,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=BREAKEND2,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=COMPLEXDEL,Description="Deletion">
##ALT=<ID=COMPLEXDUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=S,Number=1,Type=Integer,Description="Number of split reads supporting the variant">
##FORMAT=<ID=NS,Number=1,Type=Integer,Description="Number of split reads not supporting the variant">
##FORMAT=<ID=LS,Number=1,Type=Integer,Description="Number of split reads supporting the left breakpoint of this variant">
##FORMAT=<ID=LNS,Number=1,Type=Integer,Description="Number of split reads not supporting the left breakpoint of this variant">
##FORMAT=<ID=RS,Number=1,Type=Integer,Description="Number of split reads supporting the right breakpoint of this variant">
##FORMAT=<ID=RNS,Number=1,Type=Integer,Description="Number of split reads not supporting the right breakpoint of this variant">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Number of split reads not supporting the right breakpoint of this variant">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"""


def write_vcf_header(sample_name, file_name):
    """
        Write lrcnv vcf header to open file.
    """
    file_name.write(default_vcf_header + "\t" + sample_name + "\n")

def _convert_variant_type(variant_type):

    if variant_type == "simple_deletion":
        return "DEL"
    if variant_type == "simple_duplication":
        return "DUP"
    if variant_type == "complex_duplication":
        return "DUP"
    if variant_type == "complex_deletion":
        return "DEL"


def _create_info_string(cnv_row):
    """
    

    """
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
    ##INFO=<ID=UNACCOUNTED_SEQUENCE,Number=.,Type=Integer,Description="Sequence that cannot be accounted for in the comparison of the alignments, +ve novel sequence, -ve deleted sequence overall">
    ##INFO=<ID=BREAKSTART1,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##INFO=<ID=BREAKEND1,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##INFO=<ID=BREAKSTART2,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##INFO=<ID=BREAKEND2,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">

    info_dict = collections.OrderedDict() 
    info_dict["SVTYPE"] = cnv_row.variant_type
    info_dict["SVLEN"] = cnv_row.svlen
    info_dict["UNACCOUTED_SEQUENCE"] = cnv_row.unaccounted_for_sequence
    info_dict["BREAKSTART1"] = cnv_row.breakstart1
    info_dict["BREAKSTART2"] = cnv_row.breakstart2
    info_dict["BREAKEND1"] = cnv_row.breakend1
    info_dict["BREAKEND2"] = cnv_row.breakend2
    info_dict["STRANDS"] = cnv_row.strand1 + cnv_row.strand2
    out_string = "" 
    length_of_info = len(info_dict.values())
    for i, item in enumerate(info_dict.items()):
        key = item[0]
        value = str(item[1])
        if (i + 1) != length_of_info:
            out_string += key + "=" + value +";"
        else:
            out_string += key + "=" + value 
    return(out_string)

def _create_gt_string(cnv_row):
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=S,Number=1,Type=Integer,Description="Number of split reads supporting the variant">
    ##FORMAT=<ID=NS,Number=1,Type=Integer,Description="Number of split reads not supporting the variant">
    ##FORMAT=<ID=LS,Number=1,Type=Integer,Description="Number of split reads supporting the left breakpoint of this variant">
    ##FORMAT=<ID=LNS,Number=1,Type=Integer,Description="Number of split reads not supporting the left breakpoint of this variant">
    ##FORMAT=<ID=RS,Number=1,Type=Integer,Description="Number of split reads supporting the right breakpoint of this variant">
    ##FORMAT=<ID=RNS,Number=1,Type=Integer,Description="Number of split reads not supporting the right breakpoint of this variant">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Number of split reads not supporting the right breakpoint of this variant">
    gt_dict = collections.OrderedDict()
    gt_dict["GT"] = cnv_row.GT
    gt_dict["S"] = cnv_row.S
    gt_dict["NS"] = cnv_row.NS
    gt_dict["LS"] = cnv_row.LS
    gt_dict["LNS"] = cnv_row.LNS
    gt_dict["RS"] = cnv_row.RS
    gt_dict["RNS"] = cnv_row.RNS
    gt_dict["GQ"] = cnv_row.GQ
    out_string = "" 
    length_of_info = len(gt_dict.values())
    for i, item in enumerate(gt_dict.items()):
        value = str(item[1])
        if (i + 1) != length_of_info:
            out_string += value +":"
        else:
            out_string += value 
    print(out_string)
    return out_string

def write_vcf_row(cnv_row, file_name):
    """

        Write vcf row for each CNV
    """

    CHROM = cnv_row.chrom
    ID = cnv_row.event_id + "_" + cnv_row.variant_type
    # Position is breakStart1
    POS = str(cnv_row.breakStart1)
    REF = "."
    ALT = _convert_variant_type(cnv_row.variant_type)
    # FIX ME: uses sample GQ currently
    QUAL = float(cnv_row.GQ)
    if QUAL > 100:
        FILTER = "PASS"
    else:
        FILTER = "lowqual"
    QUAL = str(cnv_row.GQ)
    INFO_STRING = _create_info_string(cnv_row)
    FORMAT="GT:S:NS:LS:LNS:RS:RNS:GQ"
    GT_STRING = _create_gt_string(cnv_row)
    print("\t".join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO_STRING,FORMAT,GT_STRING]))
    file_name.write("\t".join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO_STRING,FORMAT,GT_STRING]))
    file_name.write("\n")