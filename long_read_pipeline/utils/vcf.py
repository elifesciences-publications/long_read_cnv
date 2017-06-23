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
import shutil
import gzip
import collections
import glob
import os

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
    """
        Converts the varinat type from the .cnv format to the format used in the VCFs.
    """
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
        Creates VCF info string from information saved in CNV row files.
    """
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
    """
        Creates VCF gt string for a single-sample VCF.
    """
    gt_dict = collections.OrderedDict()
    gt_dict["GT"] = cnv_row.GT
    gt_dict["S"] = cnv_row.S
    gt_dict["NS"] = cnv_row.NS
    gt_dict["LS"] = cnv_row.LS
    gt_dict["LNS"] = cnv_row.LNS
    gt_dict["RS"] = cnv_row.RS
    gt_dict["RNS"] = cnv_row.RNS
    gt_dict["GQ"] = cnv_row.GQ
    gt_dict["AB"] = cnv_row.AB
    gt_dict["SQ"] = cnv_row.SQ
    out_string = "" 
    length_of_info = len(gt_dict.values())
    for i, item in enumerate(gt_dict.items()):
        value = str(item[1])
        if (i + 1) != length_of_info:
            out_string += value +":"
        else:
            out_string += value 
    return out_string

def write_vcf_row(cnv_row, file_name,output_directory, sample_name):
    """

        Write vcf row for each CNV
    """
    novel_dir_filtered = (os.path.join(output_directory, sample_name, "event_fastas","filtered","novel"))
    reference_dir_filtered = (os.path.join(output_directory, sample_name, "event_fastas","filtered","reference"))
    novel_dir = (os.path.join(output_directory, sample_name, "event_fastas","novel"))
    reference_dir = (os.path.join(output_directory, sample_name, "event_fastas","reference"))
    try:
        os.makedirs(novel_dir_filtered)
    except OSError:
        pass
    try:
        os.makedirs(reference_dir_filtered)
    except OSError:
        pass
    CHROM = cnv_row.chrom
    # TODO: fix
    ID = cnv_row.event_id 
    # Position is breakStart1
    POS = str(cnv_row.breakStart1)
    REF = "."
    ALT = _convert_variant_type(cnv_row.variant_type)
    # FIX ME: uses sample GQ currently
    QUAL = float(cnv_row.GQ)
    if QUAL > 100:
        FILTER = "PASS"
        # Lets copy the FASTA files for both the breakpoints and the novel sequence
        # TODO: Fix hack to move folders
        novel_fasta = (os.path.join(novel_dir, cnv_row.event_id +'.fasta'))
        try:
            shutil.copy(novel_fasta, os.path.join(novel_dir_filtered, os.path.basename(novel_fasta))) 
        except IOError:
            pass
        reference_fasta = (os.path.join(reference_dir, cnv_row.event_id +'.fasta'))
        shutil.copy(reference_fasta, os.path.join(reference_dir_filtered, os.path.basename(reference_fasta))) 
    else:
        FILTER = "lowqual"
    QUAL = str(cnv_row.GQ)
    INFO_STRING = _create_info_string(cnv_row)
    FORMAT="GT:S:NS:LS:LNS:RS:RNS:GQ:AB:SQ"
    GT_STRING = _create_gt_string(cnv_row)
    print("\t".join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO_STRING,FORMAT,GT_STRING]))
    file_name.write("\t".join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO_STRING,FORMAT,GT_STRING]))
    file_name.write("\n")


class VCFSimpleRow(object):
    """
        Class representing the row of a VCF, used exclusively in the VCFSimple class.
    """

    def __init__(self, row, ids=None):
        row_s = row.strip().split("\t")
        self._CHROM = row_s[0]
        self._POS = row_s[1]
        if ids == None:
            self._ID = row_s[2]
        else:
            self._ID = ids 
        self._REF = row_s[3]
        self._ALT = row_s[4]
        self._QUAL = row_s[5]
        self._FILTER = row_s[6]
        self._INFO = row_s[7]
        self._FORMAT = row_s[8]
        self._GTSTRING = row_s[9]
        self._extra_rows = []
        self._row_count = 1
        self._info_to_dict()

    def _info_to_dict(self):
        """
            Info to dictionary
        """ 
        self._info_dict = {}
        tmp_split_info = self.info.split(";")
        for item in tmp_split_info:
            item_s = item.split("=")
            self._info_dict[item_s[0]] = item_s[1]

    @property
    def info_dict(self):
        return self._info_dict

    def add_row(self,line):
        self._extra_rows.append(VCFSimpleRow(line))
        self._row_count += 1

    @property
    def row_count(self):
        return self._row_count
    @property
    def gtstring(self):
        return self._GTSTRING
    @property
    def chrom(self):
        return self._CHROM
    @property
    def pos(self):
        return self._POS
    @property
    def ids(self):
        return self._ID
    @property
    def ref(self):
        return self._REF
    @property
    def alt(self):
        return self._ALT
    @property
    def qual(self):
        return self._QUAL
    @property 
    def filter(self):
        return self._FILTER
    @property
    def format(self):
        return self._FORMAT
    @property
    def gtstring(self):
        return self._GTSTRING
    @property
    def info(self):
        return self._INFO
    @property
    def flag(self):
        return self._FLAG
    def __str__(self):
        out_string = "\t".join(map(str,[self.chrom, self.pos,self.ids, self.ref, self.alt, self.qual, self.filter, self.info,self.format] + self._gt_string_list))
        return out_string

    def set_qual(self, qual):
        self._QUAL = int(qual)

    def set_filter(self, filter_str):
        self._FILTER = filter_str

    def set_format(self, FORMAT):
        self._FORMAT = FORMAT
    
    def set_flag(self, flag):
        self._FLAG = flag

    # Bad design this should be in a different class.
    def add_sample(self, sample, info, gt_string):
        self._INFO = info
        try:
            self._gt_string_list[0] 
        except AttributeError:
            self._gt_string_list = []
            self._sample_list = []
        self._gt_string_list.append(gt_string)
        self._sample_list.append(sample)


class VCFSimple(object):
    """
        Class representing an entire VCF.
    """
    def __init__(self,vcf):
        self._vcf_dict=collections.OrderedDict()
        self._vcf_header=""
        with gzip.open(vcf) as f:
            for line in f:
                if "#" not in line:
                    ids = line.split("\t")[2] 
                    strip_line = line.strip()
                    try:
                        ids  = self._vcf_dict[ids].ids
                        self._vcf_dict[ids].add_row(line)
                    except KeyError:
                        self._vcf_dict[ids] = VCFSimpleRow(line,ids)
                else:
                    self._vcf_header+=(line)

        
