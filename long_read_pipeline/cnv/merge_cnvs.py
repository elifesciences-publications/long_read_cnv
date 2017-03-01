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


import subprocess
import os
from long_read_pipeline.utils import vcf, fasta
from long_read_pipeline.cnv import cnv_calls 
from long_read_pipeline.speedseq import align 
import copy

def _remove_last_line_from_string(s):
        return s[:s.rfind('\n')]


# FLAG definition (used as some CNVs will be inaccurately reported.
# 0 no conflicts in sample = unique CNV
# 1 conflicts in sample at same position 
# 2 conflicts in sample at different positions
# 3 multiple conflicts at the same and different positions. 


def _process_clustering(cluster_file,all_vcfs):
    sample_list = (sorted(all_vcfs.keys()))
    header_tmp = all_vcfs[all_vcfs.keys()[0]]._vcf_header.split("\n")
    header_tmp = "\n".join(header_tmp[:(len(header_tmp)-2)])
    sample_list_simplified = [os.path.basename(x).split(".vcf")[0] for x in sample_list] 
    header_tmp += """\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT""" + "\t" + "\t".join(sample_list_simplified)
    FORMAT="GT:S:NS:LS:LNS:RS:RNS:GQ:POS:EVENTID:FLAG"
    with open(cluster_file) as f:
        prev_file = ""
        prev_id = ""
        create_row = True 
        vcf_dict_tmp = []
        flag_same_pos = False
        flag_different_pos = False
        first_sample = True
        for line in f:
            line_s = (line.strip().split("\t"))
            ids = line_s[3]
            vcf_f = line_s[4] 
            id_reads = int(line_s[5])
            if id_reads == prev_id:
                # We must be on the same cluster
                vcf_row = (all_vcfs[vcf_f]._vcf_dict[ids])
                if prev_file == vcf_f:
                    # We need to merge these INFO files, because shit is being weird.
                    if ids == ids_old:
                        tmp_vcf_info = (all_vcfs[vcf_f]._vcf_dict[ids])
                        flag_same_pos = True 
                    elif ids != ids_old:
                        tmp_vcf_info = (all_vcfs[vcf_f]._vcf_dict[ids])
                        flag_different_pos = True 
                else:
                    flag = 0
                    if flag_same_pos:
                        flag = 1 
                    if flag_different_pos:
                        flag = 2
                    if flag_same_pos and flag_different_pos:
                        flag = 3
                    all_vcfs[prev_file]._vcf_dict[ids_old].set_flag(flag) 
                    # Must be a new sample with some new info
                    # GT_STRING
                    tmp_vcf_row = all_vcfs[vcf_f]._vcf_dict[ids]
                    vcf_dict_tmp[vcf_f] = tmp_vcf_row
                    qual += float(tmp_vcf_row.qual)
                    flag_same_pos = False
                    flag_different_pos = False
            else:
                if prev_id != "":
                    flag = 0
                    if flag_same_pos:
                        flag = 1 
                    if flag_different_pos:
                        flag = 2
                    if flag_same_pos and flag_different_pos:
                        flag = 3
                    all_vcfs[prev_file]._vcf_dict[ids_old].set_flag(flag) 
                    # Finished with old VCF row
                    for sample in sample_list:
                        if sample in vcf_dict_tmp.keys():
                            vcf_tmp = vcf_dict_tmp[sample]
                            add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
                            vcf_row.add_sample(sample, vcf_dict_tmp[sample].info, vcf_dict_tmp[sample].gtstring + add_data)
                        else:
                            vcf_row.add_sample(sample, "", "0")
                    vcf_row.set_qual(qual)
                    vcf_row.set_format(FORMAT)
                    print(vcf_row)
                vcf_row = copy.deepcopy(all_vcfs[vcf_f]._vcf_dict[ids])
                vcf_dict_tmp = {} 
                vcf_dict_tmp[vcf_f] = (all_vcfs[vcf_f]._vcf_dict[ids])
                qual = float(vcf_row.qual)
                create_row = False
            prev_file  = vcf_f
            prev_id  = id_reads 
            ids_old = ids
        for sample in sample_list:
            if sample in vcf_dict_tmp.keys():
                vcf_tmp = vcf_dict_tmp[sample]
                add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids 
                vcf_row.add_sample(sample, vcf_dict_tmp[sample].info, vcf_dict_tmp[sample].gtstring + add_data)
            else:
                vcf_row.add_sample(sample, "", "0")
        vcf_row.set_qual(qual)
        vcf_row.set_format(FORMAT)
        print(vcf_row)

def merge_cnvs_indiv(vcfs, temp_dir, output_directory, in_file):
    """
        Merge CNV individual
    """
    all_vcfs = {}
    __vcf_sort__ ="vcf-sort {0} |  bgzip -c > {0}.gz && tabix -p vcf {0}.gz" 
    basenames = [os.path.basename(x).split(".vcf")[0] for x in vcfs]
    for j, vcf_f in enumerate(vcfs):
        vcf_sort_command = __vcf_sort__.format(vcf_f)
        subprocess.check_call(vcf_sort_command,shell=True)
        all_vcfs[vcf_f] = vcf.VCFSimple(vcf_f + ".gz")
            # This should give us an idea of what the variants are.
        for i in range(len(cnvs)):
            for sample in in_file:
                # Have to align
                break_point_folder = os.path.join(output_directory,basenames[j], "breakpoints")
                align.align_reads(sample,temp_dir, "test.fa", skip=True)
                #input_cnvs = (os.path.join(input_directory, sample.samples_name + ".cnv"))
                # Index and create fastas from fastq.
                fasta.index_fasta(sample, temp_dir, skip=True)
                cnv_calls.CNVFromVCF(all_vcfs[vcf_f],sample.bam_file, sample.fasta_one, sample.fasta_two, break_point_folder)
                cnvs = cnv_calls.CNVs(input_cnvs, sample.bam_file, sample.fasta_one, sample.fasta_two, break_point_folder)
                # For each locus in the sample perform lookup mapping in every other cross using this
                # samples breakpoint.


def merge_cnvs_clusters(vcfs,temp_dir, output_directory, cluster_merge_slop=0):
    """
        Merge CNVs, first step is to cluster the putatively identical deletions and duplications using bedtools. Then make a master VCF file.  
    """

    # Quality is lierally the sum of the previous VCF files.

    basenames = [os.path.basename(x) for x in vcfs]
    __bedtools_duplication_string__ = """awk '{{if ($0 !~ /^#/ && $6 > 10){{ split($8,a,";"); split(a[2],b,"=");print $1,$2,$2+b[2],$3,FILENAME}}}}' {0} | grep duplication | tr ' ' '\t' |  sort -k 1,1 -k 2,2g | bedtools cluster -i - -d {1} > tmp_clusters_duplication.txt""" 
    __bedtools_deletion_string__ = """awk '{{if ($0 !~ /^#/ && $6 > 10){{ split($8,a,";"); split(a[2],b,"=");print $1,$2,$2+b[2],$3,FILENAME}}}}' {0} | grep deletion | tr ' ' '\t' |  sort -k 1,1 -k 2,2g | bedtools cluster -i - -d {1} > tmp_clusters_deletion.txt""" 
    __bedtools_all__= """awk '{{if ($0 !~ /^#/ && $6 > 10 ){{ split($8,a,";"); split(a[2],b,"=");print $1,$2,$2+b[2],$3,FILENAME}}}}' {0} | awk '$6 > 10 {{print $0}} ' |tr ' ' '\t' |  sort -k 1,1 -k 2,2g | bedtools cluster -i - -d {1} > tmp_clusters.txt""" 
    bedtools_cmd = __bedtools_deletion_string__.format(" ".join(vcfs), cluster_merge_slop)
    subprocess.check_call(bedtools_cmd, shell=True)
    bedtools_cmd = __bedtools_duplication_string__.format(" ".join(vcfs), cluster_merge_slop)
    subprocess.check_call(bedtools_cmd, shell=True)
    bedtools_cmd = __bedtools_all__.format(" ".join(vcfs), cluster_merge_slop)
    subprocess.check_call(bedtools_cmd, shell=True)
    all_vcfs = {}
    __vcf_sort__ ="vcf-sort {0} |  bgzip -c > {0}.gz && tabix -p vcf {0}.gz" 
    for vcf_f in vcfs:
        vcf_sort_command = __vcf_sort__.format(vcf_f)
        subprocess.check_call(vcf_sort_command,shell=True)
        all_vcfs[vcf_f] = vcf.VCFSimple(vcf_f + ".gz")
    # Ok now we have all the VCFs in this format.
    # Fix final files
    #_process_clustering("tmp_clusters_duplication.txt",all_vcfs)
    _process_clustering("tmp_clusters_deletion.txt",all_vcfs)
    #_process_clustering("tmp_clusters.txt",all_vcfs)
