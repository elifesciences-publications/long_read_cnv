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

from collections import OrderedDict

import subprocess
import os
from long_read_pipeline.utils import vcf, fasta
from long_read_pipeline.cnv import cnv_calls 
from long_read_pipeline.speedseq import align 
import copy
import logging
def _remove_last_line_from_string(s):
        return s[:s.rfind('\n')]


# FLAG definition (used as some CNVs will be inaccurately reported.
# 0 no conflicts in sample = unique CNV
# 1 conflicts in sample at same position 
# 2 conflicts in sample at different positions
# 3 multiple conflicts at the same and different positions. 



def _process_clustering_pairs(cluster_file,all_vcfs, vcf_out, in_file, output_directory, temp_dir):
    # TODO: get this running
    vcf_output = open(vcf_out,"w")
    sample_list = (sorted(all_vcfs.keys()))
    header_tmp = all_vcfs[all_vcfs.keys()[0]]._vcf_header.split("\n")
    header_tmp = "\n".join(header_tmp[:(len(header_tmp)-2)])
    sample_list_simplified = [os.path.basename(x).split(".vcf")[0] for x in sample_list] 
    header_tmp += """\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT""" + "\t" + "\t".join(sample_list_simplified)
    basenames = [os.path.basename(x).split(".vcf")[0] for x in all_vcfs.keys()]
    # Need to add Sample Quality to our files
    # TODO fix format to update to using sample quality
    FORMAT="GT:S:NS:LS:LNS:RS:RNS:GQ:POS:EVENTID:FLAG:SQ"
    vcf_dict_tmp = OrderedDict() 
    for key_v in all_vcfs.keys(): 
        vcf_dict_tmp[key_v] =""
    vcf_output.write(header_tmp+"\n")
    for sample in in_file:
        fasta.index_fasta(sample, temp_dir, skip=True)
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
            vcf_f = line_s[4]  +".gz"
            id_reads = int(line_s[5])
            only_cnv = False 
            print(line_s)
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
                    only_cnv = True 
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
                    tmp_vcf_row = all_vcfs[prev_file]._vcf_dict[ids_old]
                    vcf_dict_tmp[prev_file] = tmp_vcf_row
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
                    vcf_dict_tmp[prev_file] = (all_vcfs[prev_file]._vcf_dict[ids_old])
                    all_vcfs[prev_file]._vcf_dict[ids_old].set_flag(flag) 
                    all_vcfs[prev_file]._vcf_dict[ids_old].flag
                    # Finished with old VCF row
                    qual = 0
                    tmp_vcf_row = copy.deepcopy(all_vcfs[prev_file]._vcf_dict[ids_old])
                    if vcf_dict_tmp[all_vcfs.keys()[0]] != ""  and vcf_dict_tmp[all_vcfs.keys()[1]] != "": 
                        sample = all_vcfs.keys()[0]
                        sample_other = all_vcfs.keys()[1]
                        sample_name = [basename for basename in basenames if (os.path.basename(sample).split('.')[0]) == basename][0] 
                        sample_this = [item for item in in_file if item.samples_name == sample_name][0] 
                        vcf_tmp = vcf_dict_tmp[sample]
                        add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
                        new_sample_name = [basename for basename in basenames if (os.path.basename(sample_other).split('.')[0]) == basename][0] 
                        break_point_folder_this = os.path.join(output_directory,new_sample_name,"breakpoints")
                        sample_this_new = [item for item in in_file if item.samples_name == new_sample_name][0] 
                        align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
                        cnvs = cnv_calls.CNVFromVCF(all_vcfs[sample_other],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
                        for i, cnv in enumerate(cnvs.input_rows):
                            if cnv.event_id == vcf_dict_tmp[sample_other].ids:
                                cnv._event_id = '_'.join((cnv.event_id.split("_")[:len(cnv.event_id.split("_"))-2]))
                                cnvs.extract_windowed_bam_reads_other_samples(i,sample_this.bam_file)
                                qual += float(cnvs.input_rows[i].GQ)
                                gt_string = vcf._create_gt_string(cnvs.input_rows[i])  
                                tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[0]].info, vcf_dict_tmp[sample].gtstring + add_data)
                        sample = all_vcfs.keys()[1]
                        sample_other = all_vcfs.keys()[0]
                        sample_name = [basename for basename in basenames if (os.path.basename(sample).split('.')[0]) == basename][0] 
                        sample_this = [item for item in in_file if item.samples_name == sample_name][0] 
                        vcf_tmp = vcf_dict_tmp[sample]
                        add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
                        new_sample_name = [basename for basename in basenames if (os.path.basename(sample_other).split('.')[0]) == basename][0] 
                        break_point_folder_this = os.path.join(output_directory,new_sample_name,"breakpoints")
                        sample_this_new = [item for item in in_file if item.samples_name == new_sample_name][0] 
                        align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
                        cnvs = cnv_calls.CNVFromVCF(all_vcfs[sample_other],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
                        for i, cnv in enumerate(cnvs.input_rows):
                            if cnv.event_id == vcf_dict_tmp[sample_other].ids:
                                cnv._event_id = '_'.join((cnv.event_id.split("_")[:len(cnv.event_id.split("_"))-2]))
                                cnvs.extract_windowed_bam_reads_other_samples(i,sample_this.bam_file)
                                qual += float(cnvs.input_rows[i].GQ)
                                gt_string = vcf._create_gt_string(cnvs.input_rows[i])  
                                tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[0]].info, vcf_dict_tmp[sample].gtstring + add_data)
                    elif vcf_dict_tmp[all_vcfs.keys()[1]] != "": 
                        sample = all_vcfs.keys()[1]
                        vcf_tmp = vcf_dict_tmp[sample]
                        sample_other = all_vcfs.keys()[0]
                        sample_name = [basename for basename in basenames if (os.path.basename(sample).split('.')[0]) == basename][0] 
                        sample_this = [item for item in in_file if item.samples_name == sample_name][0] 
                        align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
                        new_sample_name = [basename for basename in basenames if (os.path.basename(sample_other).split('.')[0]) == basename][0] 
                        break_point_folder_this = os.path.join(output_directory,sample_name,"breakpoints")
                        sample_this_new = [item for item in in_file if item.samples_name == new_sample_name][0] 
                        align.align_reads(sample_this_new,temp_dir, "test.fa", skip=True)
                        cnvs = cnv_calls.CNVFromVCF(all_vcfs[sample],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
                        gt_string_tmp = (vcf_dict_tmp[all_vcfs.keys()[1]].gtstring)
                        for i, cnv in enumerate(cnvs.input_rows):
                            if cnv.event_id == vcf_dict_tmp[sample].ids:
                                cnv._event_id = '_'.join((cnv.event_id.split("_")[:len(cnv.event_id.split("_"))-2]))
                                align.align_reads(sample_this_new,temp_dir, "test.fa", skip=True)
                                cnvs.extract_windowed_bam_reads_other_samples(i,sample_this_new.bam_file)
                                qual += float(cnvs.input_rows[i].GQ)
                                add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":0" 
                                gt_string = vcf._create_gt_string(cnvs.input_rows[i])
                                tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[1]].info, gt_string + add_data)
                        vcf_tmp = vcf_dict_tmp[sample]
                        qual += float(vcf_tmp.qual)
                        add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
                        tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[1]].info, gt_string_tmp + add_data)
                    elif vcf_dict_tmp[all_vcfs.keys()[0]] != "": 
                        sample = all_vcfs.keys()[0]
                        gt_string_tmp = (vcf_dict_tmp[all_vcfs.keys()[0]].gtstring)
                        sample_other = all_vcfs.keys()[1]
                        vcf_tmp = vcf_dict_tmp[sample]
                        qual += float(vcf_tmp.qual)
                        add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
                        sample_name = [basename for basename in basenames if (os.path.basename(sample).split('.')[0]) == basename][0] 
                        sample_this = [item for item in in_file if item.samples_name == sample_name][0] 
                        align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
                        vcf_tmp = vcf_dict_tmp[sample]
                        new_sample_name = [basename for basename in basenames if (os.path.basename(sample_other).split('.')[0]) == basename][0] 
                        break_point_folder_this = os.path.join(output_directory,sample_name,"breakpoints")
                        sample_this_new = [item for item in in_file if item.samples_name == new_sample_name][0] 
                        tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[0]].info, gt_string_tmp + add_data)
                        cnvs = cnv_calls.CNVFromVCF(all_vcfs[sample],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
                        for i, cnv in enumerate(cnvs.input_rows):
                            if cnv.event_id == vcf_dict_tmp[sample].ids:
                                align.align_reads(sample_this_new,temp_dir, "test.fa", skip=True)
                                cnv._event_id = '_'.join((cnv.event_id.split("_")[:len(cnv.event_id.split("_"))-2]))
                                cnvs.extract_windowed_bam_reads_other_samples(i,sample_this_new.bam_file)
                                add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":0" 
                                qual += float(cnvs.input_rows[i].GQ)
                                gt_string = vcf._create_gt_string(cnvs.input_rows[i])
                                tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[0]].info, gt_string +  add_data)
                    if (qual > 100):
                        tmp_vcf_row.set_filter("PASS")
                    else:
                        tmp_vcf_row.set_filter("lowqual")
                    tmp_vcf_row.set_qual(qual)
                    tmp_vcf_row.set_format(FORMAT)
                    vcf_output.write(str(tmp_vcf_row) +"\n")
                    flag_same_pos = False
                    flag_different_pos = False
                tmp_vcf_row = copy.deepcopy(all_vcfs[vcf_f]._vcf_dict[ids])
                vcf_dict_tmp = {} 
                for key_v in all_vcfs.keys(): 
                    vcf_dict_tmp[key_v] =""
                vcf_dict_tmp[vcf_f] = (all_vcfs[vcf_f]._vcf_dict[ids])
                qual = float(tmp_vcf_row.qual)
                create_row = False
            prev_file  = vcf_f
            prev_id  = id_reads
            if not only_cnv:
                ids_old = ids
        flag = 0
        if flag_same_pos:
            flag = 1 
        if flag_different_pos:
            flag = 2
        if flag_same_pos and flag_different_pos:
            flag = 3
        all_vcfs[prev_file]._vcf_dict[ids_old].set_flag(flag) 
        if vcf_dict_tmp[all_vcfs.keys()[0]] != ""  and vcf_dict_tmp[all_vcfs.keys()[1]] != "": 
            sample = all_vcfs.keys()[0]
            sample_other = all_vcfs.keys()[1]
            sample_name = [basename for basename in basenames if (os.path.basename(sample).split('.')[0]) == basename][0] 
            sample_this = [item for item in in_file if item.samples_name == sample_name][0] 
            vcf_tmp = vcf_dict_tmp[sample]
            add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
            new_sample_name = [basename for basename in basenames if (os.path.basename(sample_other).split('.')[0]) == basename][0] 
            break_point_folder_this = os.path.join(output_directory,new_sample_name,"breakpoints")
            sample_this_new = [item for item in in_file if item.samples_name == new_sample_name][0] 
            align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
            cnvs = cnv_calls.CNVFromVCF(all_vcfs[sample_other],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
            for i, cnv in enumerate(cnvs.input_rows):
                if cnv.event_id == vcf_dict_tmp[sample_other].ids:
                    cnv._event_id = '_'.join((cnv.event_id.split("_")[:len(cnv.event_id.split("_"))-2]))
                    cnvs.extract_windowed_bam_reads_other_samples(i,sample_this.bam_file)
                    gt_string = vcf._create_gt_string(cnvs.input_rows[i])  
                    tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[0]].info, vcf_dict_tmp[sample].gtstring + add_data)
            sample = all_vcfs.keys()[1]
            sample_other = all_vcfs.keys()[0]
            sample_name = [basename for basename in basenames if (os.path.basename(sample).split('.')[0]) == basename][0] 
            sample_this = [item for item in in_file if item.samples_name == sample_name][0] 
            vcf_tmp = vcf_dict_tmp[sample]
            add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
            new_sample_name = [basename for basename in basenames if (os.path.basename(sample_other).split('.')[0]) == basename][0] 
            break_point_folder_this = os.path.join(output_directory,new_sample_name,"breakpoints")
            sample_this_new = [item for item in in_file if item.samples_name == new_sample_name][0] 
            align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
            cnvs = cnv_calls.CNVFromVCF(all_vcfs[sample_other],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
            for i, cnv in enumerate(cnvs.input_rows):
                if cnv.event_id == vcf_dict_tmp[sample_other].ids:
                    cnv._event_id = '_'.join((cnv.event_id.split("_")[:len(cnv.event_id.split("_"))-2]))
                    cnvs.extract_windowed_bam_reads_other_samples(i,sample_this.bam_file)
                    gt_string = vcf._create_gt_string(cnvs.input_rows[i])  
                    tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[0]].info, vcf_dict_tmp[sample].gtstring + add_data)
        elif vcf_dict_tmp[all_vcfs.keys()[1]] != "": 
            sample = all_vcfs.keys()[1]
            vcf_tmp = vcf_dict_tmp[sample]
            sample_other = all_vcfs.keys()[0]
            sample_name = [basename for basename in basenames if (os.path.basename(sample).split('.')[0]) == basename][0] 
            sample_this = [item for item in in_file if item.samples_name == sample_name][0] 
            new_sample_name = [basename for basename in basenames if (os.path.basename(sample_other).split('.')[0]) == basename][0] 
            break_point_folder_this = os.path.join(output_directory,sample_name,"breakpoints")
            sample_this_new = [item for item in in_file if item.samples_name == new_sample_name][0] 
            align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
            align.align_reads(sample_this_new,temp_dir, "test.fa", skip=True)
            cnvs = cnv_calls.CNVFromVCF(all_vcfs[sample],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
            gt_string_tmp = (vcf_dict_tmp[all_vcfs.keys()[1]].gtstring)
            for i, cnv in enumerate(cnvs.input_rows):
                if cnv.event_id == vcf_dict_tmp[sample].ids:
                    cnv._event_id = '_'.join((cnv.event_id.split("_")[:len(cnv.event_id.split("_"))-2]))
                    cnvs.extract_windowed_bam_reads_other_samples(i,sample_this_new.bam_file)
                    add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":0" 
                    gt_string = vcf._create_gt_string(cnvs.input_rows[i])
                    tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[1]].info, gt_string + add_data)
            vcf_tmp = vcf_dict_tmp[sample]
            add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
            tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[1]].info, gt_string_tmp + add_data)
        elif vcf_dict_tmp[all_vcfs.keys()[0]] != "": 
            sample = all_vcfs.keys()[0]
            sample_other = all_vcfs.keys()[1]
            vcf_tmp = vcf_dict_tmp[sample]
            gt_string_tmp = (vcf_dict_tmp[all_vcfs.keys()[0]].gtstring)
            add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":" + str(vcf_tmp.flag)
            sample_name = [basename for basename in basenames if (os.path.basename(sample).split('.')[0]) == basename][0] 
            sample_this = [item for item in in_file if item.samples_name == sample_name][0] 
            vcf_tmp = vcf_dict_tmp[sample]
            new_sample_name = [basename for basename in basenames if (os.path.basename(sample_other).split('.')[0]) == basename][0] 
            break_point_folder_this = os.path.join(output_directory,sample_name,"breakpoints")
            sample_this_new = [item for item in in_file if item.samples_name == new_sample_name][0] 
            align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
            tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[0]].info, gt_string_tmp + add_data)
            cnvs = cnv_calls.CNVFromVCF(all_vcfs[sample],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
            for i, cnv in enumerate(cnvs.input_rows):
                if cnv.event_id == vcf_dict_tmp[sample].ids:
                    cnv._event_id = '_'.join((cnv.event_id.split("_")[:len(cnv.event_id.split("_"))-2]))
                    cnvs.extract_windowed_bam_reads_other_samples(i,sample_this_new.bam_file)
                    add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids  + ":0" 
                    gt_string = vcf._create_gt_string(cnvs.input_rows[i])
                    tmp_vcf_row.add_sample(sample_name, vcf_dict_tmp[all_vcfs.keys()[0]].info, gt_string +  add_data)
        tmp_vcf_row.set_qual(qual)
        tmp_vcf_row.set_format(FORMAT)
        vcf_output.write(str(tmp_vcf_row) +"\n")


def _process_clustering(cluster_file,all_vcfs, vcf_out):
    vcf_output = open(vcf_out,"w")
    sample_list = (sorted(all_vcfs.keys()))
    header_tmp = all_vcfs[all_vcfs.keys()[0]]._vcf_header.split("\n")
    header_tmp = "\n".join(header_tmp[:(len(header_tmp)-2)])
    sample_list_simplified = [os.path.basename(x).split(".vcf")[0] for x in sample_list] 
    header_tmp += """\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT""" + "\t" + "\t".join(sample_list_simplified)
    # Need to add Sample Quality to our files
    # TODO fix format to update to using sample quality
    FORMAT="GT:S:NS:LS:LNS:RS:RNS:GQ:POS:EVENTID:FLAG:SQ"
    vcf_output.write(header_tmp+"\n")
    with open(cluster_file) as f:
        prev_file = ""
        prev_id = ""
        create_row = True 
        flag_same_pos = False
        flag_different_pos = False
        first_sample = True
        for line in f:
            line_s = (line.strip().split("\t"))
            ids = line_s[3]
            vcf_f = line_s[4] 
            id_reads = int(line_s[5])
            only_cnv = False 
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
                    only_cnv = True 
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
                    tmp_vcf_row = all_vcfs[prev_file]._vcf_dict[ids_old]
                    vcf_dict_tmp[prev_file] = tmp_vcf_row
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
                    vcf_dict_tmp[prev_file] = (all_vcfs[prev_file]._vcf_dict[ids_old])
                    all_vcfs[prev_file]._vcf_dict[ids_old].set_flag(flag) 
                    all_vcfs[prev_file]._vcf_dict[ids_old].flag
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
                    vcf_output.write(str(vcf_row) +"\n")
                    flag_same_pos = False
                    flag_different_pos = False
                vcf_row = copy.deepcopy(all_vcfs[vcf_f]._vcf_dict[ids])
                vcf_dict_tmp = {} 
                vcf_dict_tmp[vcf_f] = (all_vcfs[vcf_f]._vcf_dict[ids])
                qual = float(vcf_row.qual)
                create_row = False
            prev_file  = vcf_f
            prev_id  = id_reads
            if not only_cnv:
                ids_old = ids
        flag = 0
        if flag_same_pos:
            flag = 1 
        if flag_different_pos:
            flag = 2
        if flag_same_pos and flag_different_pos:
            flag = 3
        all_vcfs[prev_file]._vcf_dict[ids_old].set_flag(flag) 
        for sample in sample_list:
            if sample in vcf_dict_tmp.keys():
                vcf_tmp = vcf_dict_tmp[sample]
                add_data = ":" + vcf_tmp.pos +":" + vcf_tmp.ids 
                vcf_row.add_sample(sample, vcf_dict_tmp[sample].info, vcf_dict_tmp[sample].gtstring + add_data)
            else:
                vcf_row.add_sample(sample, "", "0")
        vcf_row.set_qual(qual)
        vcf_row.set_format(FORMAT)
        vcf_output.write(str(vcf_row) +"\n")

def merge_cnvs_indiv(vcfs, temp_dir, output_directory, in_file):
    """
        Merge CNV individual
    """
    all_vcfs = {}
    __vcf_sort__ ="vcf-sort {0} |  bgzip -c > {0}.gz && tabix -fp vcf {0}.gz" 
    #__vcf_sort__ ="""vcf-sort {0} | awk '{{ if ($0 ~ /^#/ || $6 > 100 ){{print $0}} }}' | bgzip -c > {0}.gz && tabix -fp vcf {0}.gz""" 
    basenames = [os.path.basename(x).split(".vcf")[0] for x in vcfs]
        # TODO fix format to update to using sample quality
    FORMAT="GT:S:NS:LS:LNS:RS:RNS:GQ:POS:EVENTID:FLAG:SQ"
    logging.info("Merging all CNV calls")
    vcf_outputs  = []
    for j, vcf_f in enumerate(vcfs):
        output_vcf = os.path.join(output_directory, "vcfs", basenames[j] + ".vcf")
        output_vcf_f = open(output_vcf, "w")
        vcf_outputs.append(output_vcf)
        vcf_sort_command = __vcf_sort__.format(vcf_f)
        subprocess.check_call(vcf_sort_command,shell=True)
        all_vcfs[vcf_f] = vcf.VCFSimple(vcf_f + ".gz")
        header_tmp = all_vcfs[vcf_f]._vcf_header.split("\n")
        header_tmp = "\n".join(header_tmp[:(len(header_tmp)-2)])
        header_tmp += """\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT""" + "\t" + "\t".join(basenames)
        # Need to add Sample Quality to our files
        sample_this = [item for item in in_file if item.samples_name == basenames[j]][0] 
        align.align_reads(sample_this,temp_dir, "test.fa", skip=True)
        vcf_new_this = os.path.join(output_directory, sample_this.samples_name+".vcf.gz")
        output_vcf_f.write(header_tmp + "\n")
        fasta.index_fasta(sample_this, temp_dir, skip=True)
        break_point_folder_this = os.path.join(output_directory,basenames[j], "breakpoints")
        # VCF ids are fucking stupid, lets try and fix them here.
        cnvs = cnv_calls.CNVFromVCF(all_vcfs[vcf_f],sample_this.bam_file, sample_this.fasta_one, sample_this.fasta_two, break_point_folder_this)
        # Hello 
        logging.info("Generating merged VCF for sample = {0}".format(sample_this.samples_name))
        for i in range(len(cnvs)):
            qual = 0.0
            tmp_vcf_row = copy.deepcopy(all_vcfs[vcf_f]._vcf_dict[cnvs.input_rows[i].event_id]) 
            tmp_vcf_row.set_format(FORMAT) 
            logging.info("Processing {0}".format(str(cnvs.input_rows[i].event_id)))
            for j, sample in enumerate(in_file):
                sample_pairs = (sample.paired_end_list[0][0])
                base_fasta1 = sample_pairs.split(".gz")[0]
                sample_pairs = (sample.paired_end_list[0][1])
                base_fasta2 = sample_pairs.split(".gz")[0]
                cnvs.input_rows[i]._fasta_one = os.path.join(temp_dir, os.path.basename(base_fasta1) + ".fasta")
                cnvs.input_rows[i]._fasta_two = os.path.join(temp_dir, os.path.basename(base_fasta2) + ".fasta") 
                # Have to align
                if sample.samples_name == sample_this.samples_name:
                    cnvs.extract_windowed_bam_reads(i)
                    gt_string = vcf._create_gt_string(cnvs.input_rows[i])  
                    tmp_vcf_row.add_sample(sample.samples_name, tmp_vcf_row.info, gt_string)
                else:
                    align.align_reads(sample,temp_dir, "test.fa", skip=True)
                    #input_cnvs = (os.path.join(input_directory, sample.samples_name + ".cnv"))
                    # Index and create fastas from fastq.
                    cnvs.extract_windowed_bam_reads_other_samples(i,sample.bam_file)
                    gt_string = vcf._create_gt_string(cnvs.input_rows[i])  
                    tmp_vcf_row.add_sample(sample.samples_name, tmp_vcf_row.info, gt_string)
                qual += float(cnvs.input_rows[i].GQ)
            tmp_vcf_row.set_qual(qual)
                # For each locus in the sample perform lookup mapping in every other cross using this
                # samples breakpoint.
            output_vcf_f.write(str(tmp_vcf_row)+ "\n")
        output_vcf_f.close()

def pairwise_cnvs(vcfs,temp_dir, output_directory, in_file, cluster_merge_slop=0):
    """
        Merge CNVs, first step is to cluster the putatively identical deletions and duplications using bedtools. Then make a master VCF file.  
    """

    # Quality is lierally the sum of the previous VCF files.

    vcf_pairs = []
    for i in range(len(vcfs)):
        for j in range(i+1,len(vcfs)):
            vcf_pairs.append((vcfs[i],vcfs[j]))
    for p1, p2 in vcf_pairs: 
        print(p1)
        print(p2)
        __bedtools_all__= """awk '{{if ($0 !~ /^#/ && $6 > 10 ){{ split($8,a,";"); split(a[2],b,"=");print $1,$2,$2+b[2],$3,FILENAME}}}}' {0}  |tr ' ' '\t' |  sort -k 1,1 -k 2,2g | bedtools cluster -i - -d {1} > tmp_clusters.txt""" 
        bedtools_cmd = __bedtools_all__.format(" ".join([p1,p2]), cluster_merge_slop)
        subprocess.check_call(bedtools_cmd, shell=True)
        __vcf_sort__ ="vcf-sort {0} |  bgzip -c > {0}.gz && tabix -p vcf {0}.gz" 
        vcf_sort_command = __vcf_sort__.format(p1)
        subprocess.check_call(vcf_sort_command, shell=True)
        __vcf_sort__ ="vcf-sort {0} |  bgzip -c > {0}.gz && tabix -p vcf {0}.gz" 
        vcf_sort_command = __vcf_sort__.format(p2)
        subprocess.check_call(vcf_sort_command, shell=True)
        p1 = p1 + ".gz"
        p2 = p2 + ".gz"
        all_vcfs = {}
        all_vcfs[p1] = vcf.VCFSimple(p1)
        all_vcfs[p2] = vcf.VCFSimple(p2)
        try:
            os.mkdir(os.path.join(output_directory, "paired_vcfs"))
        except OSError:
            pass
        output_file = os.path.join(output_directory,"paired_vcfs", os.path.basename(p1) + os.path.basename(p2) + ".vcf")
        _process_clustering_pairs("tmp_clusters.txt",all_vcfs, output_file,in_file, output_directory, temp_dir)

def merge_cnvs_clusters(vcfs,temp_dir, output_directory, cluster_merge_slop=0):
    """
        Merge CNVs, first step is to cluster the putatively identical deletions and duplications using bedtools. Then make a master VCF file.  
    """

    # Quality is lierally the sum of the previous VCF files.

    basenames = [os.path.basename(x) for x in vcfs]
    __bedtools_duplication_string__ = """awk '{{if ($0 !~ /^#/ && $6 > 100){{ split($8,a,";"); split(a[2],b,"=");print $1,$2,$2+b[2],$3,FILENAME}}}}' {0} | grep duplication | tr ' ' '\t' |  sort -k 1,1 -k 2,2g | bedtools cluster -i - -d {1} > tmp_clusters_duplication.txt""" 
    __bedtools_deletion_string__ = """awk '{{if ($0 !~ /^#/ && $6 > 100){{ split($8,a,";"); split(a[2],b,"=");print $1,$2,$2+b[2],$3,FILENAME}}}}' {0} | grep deletion | tr ' ' '\t' |  sort -k 1,1 -k 2,2g | bedtools cluster -i - -d {1} > tmp_clusters_deletion.txt""" 
    __bedtools_all__= """awk '{{if ($0 !~ /^#/ && $6 > 100 ){{ split($8,a,";"); split(a[2],b,"=");print $1,$2,$2+b[2],$3,FILENAME}}}}' {0}  |tr ' ' '\t' |  sort -k 1,1 -k 2,2g | bedtools cluster -i - -d {1} > tmp_clusters.txt""" 
    bedtools_cmd = __bedtools_deletion_string__.format(" ".join(vcfs), cluster_merge_slop)
    subprocess.check_call(bedtools_cmd, shell=True)
    bedtools_cmd = __bedtools_duplication_string__.format(" ".join(vcfs), cluster_merge_slop)
    subprocess.check_call(bedtools_cmd, shell=True)
    bedtools_cmd = __bedtools_all__.format(" ".join(vcfs), cluster_merge_slop)
    subprocess.check_call(bedtools_cmd, shell=True)
    all_vcfs = {}
    __vcf_sort__ ="vcf-sort {0} |  bgzip -c > {0}.gz && tabix -fp vcf {0}.gz" 
    for vcf_f in vcfs:
        vcf_sort_command = __vcf_sort__.format(vcf_f)
        subprocess.check_call(vcf_sort_command,shell=True)
        all_vcfs[vcf_f] = vcf.VCFSimple(vcf_f + ".gz")
    # Ok now we have all the VCFs in this format.
    # Fix final files
    try:
        os.mkdir(os.path.join(output_directory,"vcfs"))
    except OSError:
        pass
    _process_clustering("tmp_clusters_duplication.txt",all_vcfs,os.path.join(output_directory,"vcfs","duplications.vcf"))
    _process_clustering("tmp_clusters_deletion.txt",all_vcfs, os.path.join(output_directory, "vcfs","deletions.vcf"))
    _process_clustering("tmp_clusters.txt",all_vcfs, os.path.join(output_directory,"vcfs", "all.vcf"))
