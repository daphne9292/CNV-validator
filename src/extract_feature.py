#coding:utf-8
import pandas as pd
import pysam

import re
import math
import os
import multiprocessing
import subprocess
import numpy as np
def parse_cigar(cigar_str):
    """
    Parse a CIGAR string into tuples of (operation, length).

    Parameters:
    cigar_str (str): The CIGAR string to parse.
    """
    # 使用正则表达式拆分CIGAR字符串
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    cigar_operations = cigar_pattern.findall(cigar_str)

    # 打印拆分后的CIGAR操作
    # for length, operation in cigar_operations:
    #     print(f"Length: {length}, Operation: {operation}")
    return cigar_operations

def wc_count(file_name):
    import subprocess
    buffer= 1024*1024
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])

def get_samtools_depth(chrom,start,end,bam):
    linux_command = 'samtools depth -aa -r ' + chrom + ':' +str(start) + '-' + str(end) + ' '+ bam
    num = 0
    sum_depth = 0
    ave_depth = 0
    depth=[]
    with os.popen(linux_command) as command:
        for line in command:
            info = line.split()
            depth.append(int(info[2]))

    step = int(len(depth)/20)
    for i in range(0,len(depth),step):
        num = num + 1
        sum_depth += depth[i]
    if num > 0:
        ave_depth = sum_depth / num
    return ave_depth


def fetch_reads_nineth(bam_file, chrom, start, end,type,mean_insert_size = 500,std_insert_size = 50,window_size = 1000):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    expand_size = mean_insert_size
    abnormal_size = 1000
    if (end - start ) < 1500:
        abnormal_size = mean_insert_size + 2 * std_insert_size
    left = 0
    right = 0
    #training
    # if type == 'DUP':
    #     left = start
    #     right = end
    # elif type =='DEL':
    #     left = max(0, start - expand_size)
    #     right = end + expand_size

    total_read_num_left_bp_interval = 0
    total_read_num_right_bp_interval = 0
    valid_RP_num_left_bp_interval = 0
    valid_RP_num_right_bp_interval = 0
    valid_SR_num_left_bp_interval = 0
    valid_SR_num_right_bp_interval = 0
    k = 0
    #在左右断点expand_size范围内找寻异常RP、SR信号
    # methods callset
    left = max(start - window_size,0)
    for read in bam.fetch(chrom, left, start + window_size):
    # for read in bam.fetch(chrom, left,left+expand_size):
        k += 1
        if k > 4000:
            break
        if read.mapping_quality < 20:
            continue
        if (not read.is_paired) or read.mate_is_unmapped:
            continue
        mate_read = bam.mate(read)
        if not mate_read:
            continue
        insert_size = read.template_length
        cigar_operations = parse_cigar(str(read.cigarstring))
        index_max_S = -1
        index_max_M = -1
        for i in range(len(cigar_operations)):
            operation = cigar_operations[i][1]
            length = int(cigar_operations[i][0])
            if (operation == 'S' or operation == 'H') and (
                    index_max_S < 0 or length > int(cigar_operations[index_max_S][0])):
                index_max_S = i
            if operation == 'M' and (index_max_M < 0 or length > int(cigar_operations[index_max_M][0])):
                index_max_M = i
        total_read_num_left_bp_interval +=1
        if ('S' in str(read.cigarstring) or 'H' in str(read.cigarstring)) and read.has_tag('SA'):
            chrom = read.reference_name
            sa_tag = read.get_tag('SA')
            multi_sa_alignment = sa_tag.split(';')[0:-1]
            for item in multi_sa_alignment:
                detail = item.split(',')
                supplementary_mapping_quality = int(detail[4])
                if supplementary_mapping_quality >= 30:
                    valid_SR_num_left_bp_interval += 1
        if ((abs(insert_size) > abnormal_size or insert_size == 0) and (not mate_read.is_reverse) and (type == 'DUP')) or \
                ((type == 'DEL') and abs(insert_size) > abnormal_size):
            # print(read.query_name)
            mate_read_index_max_S = -1
            mate_read_index_max_M = -1
            mate_cigar_operations = parse_cigar(str(mate_read.cigarstring))
            for i in range(len(mate_cigar_operations)):
                operation = mate_cigar_operations[i][1]
                length = int(mate_cigar_operations[i][0])
                if (operation == 'S' or operation == 'H') and (mate_read_index_max_S < 0 or length > int(mate_cigar_operations[mate_read_index_max_S][0])):
                    mate_read_index_max_S = i
                if operation == 'M' and (mate_read_index_max_M < 0 or length > int(mate_cigar_operations[mate_read_index_max_M][0])):
                    mate_read_index_max_M = i

            if int(cigar_operations[index_max_M][0])/read.reference_length > 0.3 and \
                int(mate_cigar_operations[mate_read_index_max_M][0])/mate_read.reference_length > 0.3:
                valid_RP_num_left_bp_interval += 1

    # methods callset
    k = 0
    for read in bam.fetch(chrom, end - window_size, end + window_size):
    # for read in bam.fetch(chrom, right - expand_size,right):
        k += 1
        if k > 4000:
            break
        if read.mapping_quality < 20:
            continue
        if (not read.is_paired) or read.mate_is_unmapped:
            continue
        mate_read = bam.mate(read)
        if not mate_read:
            continue
        insert_size = read.template_length
        cigar_operations = parse_cigar(str(read.cigarstring))
        index_max_S = -1
        index_max_M = -1
        for i in range(len(cigar_operations)):
            operation = cigar_operations[i][1]
            length = int(cigar_operations[i][0])
            if (operation == 'S' or operation == 'H') and (
                    index_max_S < 0 or length > int(cigar_operations[index_max_S][0])):
                index_max_S = i
            if operation == 'M' and (index_max_M < 0 or length > int(cigar_operations[index_max_M][0])):
                index_max_M = i
        total_read_num_right_bp_interval += 1
        if ('S' in str(read.cigarstring) or 'H' in str(read.cigarstring)) and read.has_tag('SA'):
            chrom = read.reference_name

            sa_tag = read.get_tag('SA')
            multi_sa_alignment = sa_tag.split(';')[0:-1]
            for item in multi_sa_alignment:
                detail = item.split(',')
                supplementary_mapping_quality = int(detail[4])
                if supplementary_mapping_quality >= 30:
                    valid_SR_num_right_bp_interval += 1

        if ((abs(insert_size) > abnormal_size or insert_size == 0) and (not read.is_reverse) and type == 'DUP') or \
                (abs(insert_size) > abnormal_size and type == 'DEL'):
            # print(read.query_name)
            mate_read_index_max_S = -1
            mate_read_index_max_M = -1
            mate_cigar_operations = parse_cigar(str(mate_read.cigarstring))
            for i in range(len(mate_cigar_operations)):
                operation = mate_cigar_operations[i][1]
                length = int(mate_cigar_operations[i][0])
                if (operation == 'S' or operation == 'H') and (mate_read_index_max_S < 0 or length > int(mate_cigar_operations[mate_read_index_max_S][0])):
                    mate_read_index_max_S = i
                if operation == 'M' and (mate_read_index_max_M < 0 or length > int(mate_cigar_operations[mate_read_index_max_M][0])):
                    mate_read_index_max_M = i

            if int(cigar_operations[index_max_M][0]) / read.reference_length > 0.3 and \
                    int(mate_cigar_operations[mate_read_index_max_M][0]) / mate_read.reference_length > 0.3:
                valid_RP_num_right_bp_interval += 1

    # Close the BAM file
    bam.close()
    return total_read_num_left_bp_interval,total_read_num_right_bp_interval,valid_RP_num_left_bp_interval,valid_RP_num_right_bp_interval,\
        valid_SR_num_left_bp_interval,valid_SR_num_right_bp_interval

def fetch_snp_nineth(fasta,chrom, start, end,snp_vcf,normal_bam_dir,tumor_bam_dir):
    # fasta = "/data/DATA/Reference/human/GRCh38.d1.vd1/genome/GRCh38.d1.vd1.fa"
    read_fa = pysam.FastaFile(fasta)
    chrom_dict = {'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9,
                  'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17,
                  'chr18': 18, 'chr19': 19, 'chr20': 20, 'chr21': 21, 'chr22': 22, 'chrX': 23, 'chrY': 24}
    chrom2index = chrom_dict[chrom] - 1
    chrom_length = read_fa.lengths[chrom2index]
    #根据区间尺寸左右外扩expand_size
    svlen = end - start
    if svlen < 10000:
        expand_size = 3000
    elif svlen < 100000:
        expand_size = min(int(svlen * 0.3), 8000)
    else:
        expand_size = 8000
    left = max(1,start - expand_size)
    right = min(chrom_length,end + expand_size)

    # 文件太大，根据区间坐标猜测一个位置实现文件快速读取
    total_line = wc_count(snp_vcf)
    # line_offset = 2
    line_offset = math.ceil(left / chrom_length * total_line)
    if line_offset > 1 and total_line > 1 and line_offset <= total_line:

        res= subprocess.run(['sed','-n','{}p'.format(line_offset),snp_vcf],stdout=subprocess.PIPE,text=True,check=True)
        line = res.stdout.strip()
        detail = line.split('\t',2)

        if int(detail[1]) < left:
            for j in range(line_offset +50,total_line,50):
                res = subprocess.run(['sed', '-n', '{}p'.format(j), snp_vcf], stdout=subprocess.PIPE, text=True,
                                     check=True)
                line = res.stdout.strip()
                detail = line.split('\t', 2)
                if int(detail[1]) < left:
                    line_offset += 50
                else:
                    break
        # print(line_offset)
        res = subprocess.run(['sed', '-n', '{}p'.format(line_offset), snp_vcf], stdout=subprocess.PIPE, text=True,
                             check=True)
        line = res.stdout.strip()
        detail = line.split('\t', 2)
        while int(detail[1]) > left:
            line_offset = max(line_offset - 50,1)
            res = subprocess.run(['sed', '-n', '{}p'.format(line_offset), snp_vcf], stdout=subprocess.PIPE, text=True,
                                 check=True)
            line = res.stdout.strip()
            detail = line.split('\t', 2)
            if line_offset == 1:
                break

    ave_tumor_baf_leftside_of_interval = 0
    ave_normal_baf_leftside_of_interval = 0
    ave_tumor_baf_of_interval = 0
    ave_normal_baf_of_interval = 0
    ave_tumor_baf_rightside_of_interval = 0
    ave_normal_baf_rightside_of_interval = 0
    num4 = 0
    num5 = 0
    num6 = 0
    if line_offset == 1:
        line_offset = 2
    for i in range(line_offset,total_line):
        res = subprocess.run(['sed', '-n', '{}p'.format(i), snp_vcf], stdout=subprocess.PIPE, text=True,
                             check=True)

        line = res.stdout.strip()
        detail = line.split('\t', 2)
        if int(detail[1]) > right:
            break
        elif int(detail[1]) < left:
            continue
        else:
            info = line.split()
            normal_baf = float(info[6].strip("%"))/100
            tumor_baf =float(info[10].strip("%"))/100
            if int(detail[1]) < start:
                if abs(normal_baf - 0.5) < 0.15:
                    ave_tumor_baf_leftside_of_interval += tumor_baf
                    ave_normal_baf_leftside_of_interval += normal_baf
                    num4 += 1
            elif int(detail[1]) <= end:
                if abs(normal_baf - 0.5) < 0.15:
                    ave_tumor_baf_of_interval += tumor_baf
                    ave_normal_baf_of_interval += normal_baf
                    num5 += 1
            elif int(detail[1]) <= right:
                if abs(normal_baf - 0.5) < 0.15:
                    ave_tumor_baf_rightside_of_interval += tumor_baf
                    ave_normal_baf_rightside_of_interval += normal_baf
                    num6 += 1


    ave_t_depth_leftside_of_interval = get_samtools_depth(chrom,left,start - 1,tumor_bam_dir)
    ave_n_depth_leftside_of_interval = get_samtools_depth(chrom, left, start - 1, normal_bam_dir)
    ave_t_depth_of_interval = get_samtools_depth(chrom,start,end,tumor_bam_dir)
    ave_n_depth_of_interval = get_samtools_depth(chrom,start,end,normal_bam_dir)
    ave_t_depth_rightside_of_interval = get_samtools_depth(chrom,end+1,right,tumor_bam_dir)
    ave_n_depth_rightside_of_interval = get_samtools_depth(chrom,end+1,right,normal_bam_dir)

    t_depth_ratio_of_leftside_with_interval = ave_t_depth_of_interval / ave_t_depth_leftside_of_interval if ave_t_depth_leftside_of_interval > 0 else 1
    t_depth_ratio_of_rightside_with_interval = ave_t_depth_of_interval / ave_t_depth_rightside_of_interval if ave_t_depth_rightside_of_interval > 0 else 1
    n_depth_ratio_of_leftside_with_interval = ave_n_depth_of_interval / ave_n_depth_leftside_of_interval if ave_n_depth_leftside_of_interval > 0 else 1
    n_depth_ratio_of_rightside_with_interval = ave_n_depth_of_interval / ave_n_depth_rightside_of_interval if ave_n_depth_rightside_of_interval > 0 else 1

    if num4 > 0 :
        ave_tumor_baf_leftside_of_interval = ave_tumor_baf_leftside_of_interval/num4
        ave_normal_baf_leftside_of_interval = ave_normal_baf_leftside_of_interval/num4
    if num5 > 0:
        ave_tumor_baf_of_interval = ave_tumor_baf_of_interval / num5
        ave_normal_baf_of_interval=ave_normal_baf_of_interval/num5
    if num6 > 0:
        ave_tumor_baf_rightside_of_interval = ave_tumor_baf_rightside_of_interval / num6
        ave_normal_baf_rightside_of_interval = ave_normal_baf_rightside_of_interval/num6

    t_BAF_diff_of_leftside_with_interval = abs(ave_tumor_baf_leftside_of_interval - ave_tumor_baf_of_interval)
    t_BAF_diff_of_rightside_with_interval = abs(ave_tumor_baf_rightside_of_interval - ave_tumor_baf_of_interval)
    n_BAF_diff_of_leftside_with_interval = abs(ave_normal_baf_leftside_of_interval - ave_normal_baf_of_interval)
    n_BAF_diff_of_rightside_with_interval = abs(ave_normal_baf_rightside_of_interval - ave_normal_baf_of_interval)
    return t_depth_ratio_of_leftside_with_interval,t_depth_ratio_of_rightside_with_interval,\
              n_depth_ratio_of_leftside_with_interval,n_depth_ratio_of_rightside_with_interval,\
              t_BAF_diff_of_leftside_with_interval,t_BAF_diff_of_rightside_with_interval,\
        n_BAF_diff_of_leftside_with_interval,n_BAF_diff_of_rightside_with_interval

def chunker_RP_SR_nineth(list,RP_SR_out_dir,n_bam,t_bam,mean_insert_size,std_insert_size,sample_name='s'):
    for i in range(len(list)):
        chrom = list[i][0]
        start = int(list[i][1])
        end = int(list[i][2])
        type = list[i][3]
        n_total_read_num_left_bp_interval, n_total_read_num_right_bp_interval, n_valid_RP_num_left_bp_interval, n_valid_RP_num_right_bp_interval, \
            n_valid_SR_num_left_bp_interval, n_valid_SR_num_right_bp_interval = fetch_reads_nineth(n_bam, chrom, start, end,type,mean_insert_size,std_insert_size)
        t_total_read_num_left_bp_interval, t_total_read_num_right_bp_interval, t_valid_RP_num_left_bp_interval, t_valid_RP_num_right_bp_interval, \
            t_valid_SR_num_left_bp_interval, t_valid_SR_num_right_bp_interval = fetch_reads_nineth(t_bam, chrom, start, end,type,mean_insert_size,std_insert_size)
        out_dir = RP_SR_out_dir + sample_name + '_'+chrom+'_'+str(start)+'_'+str(end)+'_'+list[i][3]+'.proof.txt'

        with open(out_dir, "w") as fout:
            print(n_total_read_num_left_bp_interval, n_total_read_num_right_bp_interval, n_valid_RP_num_left_bp_interval, n_valid_RP_num_right_bp_interval, \
            n_valid_SR_num_left_bp_interval, n_valid_SR_num_right_bp_interval,t_total_read_num_left_bp_interval, t_total_read_num_right_bp_interval,\
                  t_valid_RP_num_left_bp_interval, t_valid_RP_num_right_bp_interval, t_valid_SR_num_left_bp_interval, t_valid_SR_num_right_bp_interval, file=fout)
def chunker_depth_BAF_nineth(list,depth_BAF_out_dir,normal_bam_dir,tumor_bam_dir,fasta,snp_vcf_prefix,sample_name='s'):
    for i in range(len(list)):
        chrom = list[i][0]
        start = int(list[i][1])
        end = int(list[i][2])
        type = list[i][3]
        snp_vcf = snp_vcf_prefix + chrom + ".snv.vcf.snp"
        out_dir = depth_BAF_out_dir + sample_name+'_'+chrom + '_' + str(start) + '_' + str(end) + '_' + type + '.proof.txt'

        t_depth_ratio_of_leftside_with_interval, t_depth_ratio_of_rightside_with_interval, \
        n_depth_ratio_of_leftside_with_interval, n_depth_ratio_of_rightside_with_interval, \
        t_BAF_diff_of_leftside_with_interval, t_BAF_diff_of_rightside_with_interval, \
        n_BAF_diff_of_leftside_with_interval, n_BAF_diff_of_rightside_with_interval = fetch_snp_nineth(fasta,chrom, start,end, snp_vcf,normal_bam_dir,tumor_bam_dir)
        with open(out_dir, "w") as fout:
            print(t_depth_ratio_of_leftside_with_interval, t_depth_ratio_of_rightside_with_interval, \
                  n_depth_ratio_of_leftside_with_interval, n_depth_ratio_of_rightside_with_interval, \
                  t_BAF_diff_of_leftside_with_interval, t_BAF_diff_of_rightside_with_interval, \
                  n_BAF_diff_of_leftside_with_interval, n_BAF_diff_of_rightside_with_interval, file=fout)
def PR_SR_nineth(benchmark_bed,RP_SR_dir,normal_bam,tumor_bam,mean_insert_size,std_insert_size,thread_num=1,sample_name='s'):
    bed_list = []
    # bed需要四列，最后一列是类型‘DUP’或‘DEL'
    with open(benchmark_bed, 'r') as fin:
        for line in fin:
            info = line.strip().split()
            sta = int(info[1])
            end = int(info[2])
            svlen = end - sta
            if len(info) == 4:
                if info[3] == 'LOH':
                    continue
                else:
                    type = info[3]
                proof = RP_SR_dir + sample_name + '_' + info[0] + '_' + info[1] + '_' + info[2] + '_' + type + '.proof.txt'
                if svlen > 1000 and not os.path.exists(proof):
                    bed_list.append([info[0], info[1], info[2], type])

    print(len(bed_list))

    total_row = len(bed_list)
    num_chunks = thread_num
    chunk_size = total_row // num_chunks + 1
    jobs = []

    for i in range(num_chunks):
        list_chunk = [row for row in bed_list[i * chunk_size: min((i + 1) * chunk_size, total_row)]]
        p = multiprocessing.Process(target=chunker_RP_SR_nineth, args=(list_chunk, RP_SR_dir,normal_bam,tumor_bam,mean_insert_size,std_insert_size,sample_name))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

def depth_BAF_nineth(benchmark_bed,depth_BAF_dir,normal_bam_dir,tumor_bam_dir,fasta,snp_vcf_prefix,thread_num=1,sample_name='s'):
    bed_list =[]
    with open(benchmark_bed,'r') as fin:
        for line in fin:
            info = line.strip().split()
            sta = int(info[1])
            end = int(info[2])
            svlen = end - sta
            if len(info) == 4:
                if info[3] == 'LOH':
                    continue
                else:
                    type = info[3]

                proof = depth_BAF_dir + sample_name + '_' + info[0] + '_' + info[1] + '_' + info[2] + '_' + type + '.proof.txt'
                if svlen > 1000 and not os.path.exists(proof):
                    bed_list.append([info[0], info[1], info[2], type])

    print(len(bed_list))

    # print(bed_list[0])
    total_row=len(bed_list)
    num_chunks = thread_num
    chunk_size = total_row // num_chunks + 1
    jobs = []

    for i in range(num_chunks):
        list_chunk = [row for row in bed_list[i * chunk_size: min((i + 1) * chunk_size,total_row)]]
        p = multiprocessing.Process(target=chunker_depth_BAF_nineth, args=(list_chunk,depth_BAF_dir,normal_bam_dir,tumor_bam_dir,fasta,snp_vcf_prefix,sample_name))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()


    thread = 48
    fasta = "/data/DATA/Reference/human/GRCh38.d1.vd1/genome/GRCh38.d1.vd1.fa"
    negative_bed = "/data/home/gyzhe/project/validateCNVBeakpoint/data/ChineseQuartet/analysis/methods.call.cnv/cnvkit/LCL6normal/proof/totalcnv.bed"
    RP_SR_dir = "/data/home/gyzhe/project/validateCNVBeakpoint/data/ChineseQuartet/analysis/methods.call.cnv/cnvkit/LCL6normal/proof/RP_and_SR_proof/"
    normal_bam = "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL6/ChineseQuartet.LCL6.GRCh38.ILM.bwa.markdup.bam"
    bam_file = '/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL5/ChineseQuartet.LCL5.GRCh38.ILM.bwa.markdup.bam'
    mean_insert_size = 427
    std_insert_size = 136
    s_name = 'LCL5.cnvkitcall'
    snp_vcf_prefix = "/data/home/gyzhe/project/validateCNVBeakpoint/data/ChineseQuartet/analysis/varscan/ChineseQuartet.GRCh38.ILM.varscan."
    PR_SR_nineth(negative_bed, RP_SR_dir, normal_bam, bam_file,mean_insert_size,std_insert_size,thread,s_name)
    depth_BAF_dir = "/data/home/gyzhe/project/validateCNVBeakpoint/data/ChineseQuartet/analysis/methods.call.cnv/cnvkit/LCL6normal/proof/depth_and_BAF_proof/"
    depth_BAF_nineth(negative_bed, depth_BAF_dir, normal_bam, bam_file,fasta,snp_vcf_prefix,thread,s_name)
def statisticproof(RP_SR_proof_dir,depth_BAF_proof_dir, out_statistic_csv):
    RP_SR_depth_BAF_statistic_df = pd.DataFrame()
    id = []
    report_len = []
    sv_type = []
    t_depth_ratio_of_leftside_with_interval=[]
    t_depth_ratio_of_rightside_with_interval=[]
    n_depth_ratio_of_leftside_with_interval=[]
    n_depth_ratio_of_rightside_with_interval=[]
    t_BAF_diff_of_leftside_with_interval=[]
    t_BAF_diff_of_rightside_with_interval=[]
    n_BAF_diff_of_leftside_with_interval=[]
    n_BAF_diff_of_rightside_with_interval=[]
    valid_t_RP_percent_left_bp_interval = []
    valid_t_SR_percent_left_bp_interval = []
    valid_t_RP_percent_right_bp_interval = []
    valid_t_SR_percent_right_bp_interval = []
    valid_n_RP_percent_left_bp_interval = []
    valid_n_SR_percent_left_bp_interval = []
    valid_n_RP_percent_right_bp_interval = []
    valid_n_SR_percent_right_bp_interval = []

    RP_SR_dir = RP_SR_proof_dir
    depth_BAF_dir = depth_BAF_proof_dir
    files = os.listdir(RP_SR_dir)

    for file in files:
        head = file.split('.')[0]
        filename = head.split('_')

        start = int(filename[1])
        end = int(filename[2])
        svlen = end - start
        type = filename[3]
        RP_SR_file_path = RP_SR_dir + file
        if os.path.isfile(RP_SR_file_path):
            with open(RP_SR_file_path, 'r') as fin:
                for line in fin:
                    info = line.split()
                    n_total_read_num_left_bp_interval = int(info[0])
                    n_total_read_num_right_bp_interval = int(info[1])
                    n_valid_RP_num_left_bp_interval = int(info[2])
                    n_valid_RP_num_right_bp_interval = int(info[3])
                    n_valid_SR_num_left_bp_interval = int(info[4])
                    n_valid_SR_num_right_bp_interval = int(info[5])
                    t_total_read_num_left_bp_interval = int(info[6])
                    t_total_read_num_right_bp_interval = int(info[7])
                    t_valid_RP_num_left_bp_interval = int(info[8])
                    t_valid_RP_num_right_bp_interval = int(info[9])
                    t_valid_SR_num_left_bp_interval = int(info[10])
                    t_valid_SR_num_right_bp_interval = int(info[11])
            # 将RP\SR证据转换为一些绝对数值和异常reads百分比的特征
            valid_n_RP_percent_left_bp_interval.append(n_valid_RP_num_left_bp_interval/n_total_read_num_left_bp_interval \
                                                         if n_total_read_num_left_bp_interval > 0 else np.nan)
            valid_n_SR_percent_left_bp_interval.append(n_valid_SR_num_left_bp_interval/n_total_read_num_left_bp_interval \
                                                         if n_total_read_num_left_bp_interval > 0 else np.nan)
            valid_n_RP_percent_right_bp_interval.append(n_valid_RP_num_right_bp_interval/n_total_read_num_right_bp_interval \
                                                      if n_total_read_num_right_bp_interval > 0 else np.nan)
            valid_n_SR_percent_right_bp_interval.append(n_valid_SR_num_right_bp_interval/n_total_read_num_right_bp_interval \
                                                      if n_total_read_num_right_bp_interval > 0 else np.nan)
            valid_t_RP_percent_left_bp_interval.append(t_valid_RP_num_left_bp_interval/t_total_read_num_left_bp_interval \
                                                         if t_total_read_num_left_bp_interval > 0 else np.nan)
            valid_t_SR_percent_left_bp_interval.append(t_valid_SR_num_left_bp_interval/t_total_read_num_left_bp_interval \
                                                         if t_total_read_num_left_bp_interval > 0 else np.nan)
            valid_t_RP_percent_right_bp_interval.append(t_valid_RP_num_right_bp_interval/t_total_read_num_right_bp_interval \
                                                      if t_total_read_num_right_bp_interval > 0 else np.nan)
            valid_t_SR_percent_right_bp_interval.append(t_valid_SR_num_right_bp_interval/t_total_read_num_right_bp_interval \
                                                      if t_total_read_num_right_bp_interval > 0 else np.nan)

            id.append(file.split('.', 1)[0])
            report_len.append(svlen)
            sv_type.append(type)

        depth_BAF_file_path = depth_BAF_dir + file
        if os.path.isfile(depth_BAF_file_path):
            with open(depth_BAF_file_path, 'r') as fin:
                for line in fin:
                    info = line.split()
                    t_depth_ratio_of_leftside_with_interval.append(float(info[0]))
                    t_depth_ratio_of_rightside_with_interval.append(float(info[1]))
                    n_depth_ratio_of_leftside_with_interval.append(float(info[2]))
                    n_depth_ratio_of_rightside_with_interval.append(float(info[3]))
                    t_BAF_diff_of_leftside_with_interval.append(float(info[4]))
                    t_BAF_diff_of_rightside_with_interval.append(float(info[5]))
                    n_BAF_diff_of_leftside_with_interval.append(float(info[6]))
                    n_BAF_diff_of_rightside_with_interval.append(float(info[7]))

    RP_SR_depth_BAF_statistic_df['id'] = id
    RP_SR_depth_BAF_statistic_df['report_len'] = report_len
    RP_SR_depth_BAF_statistic_df['svtype'] = sv_type
    RP_SR_depth_BAF_statistic_df['t_depth_ratio_of_leftside_with_interval'] = t_depth_ratio_of_leftside_with_interval
    RP_SR_depth_BAF_statistic_df['t_depth_ratio_of_rightside_with_interval'] = t_depth_ratio_of_rightside_with_interval
    RP_SR_depth_BAF_statistic_df['n_depth_ratio_of_leftside_with_interval'] = n_depth_ratio_of_leftside_with_interval
    RP_SR_depth_BAF_statistic_df['n_depth_ratio_of_rightside_with_interval'] = n_depth_ratio_of_rightside_with_interval
    RP_SR_depth_BAF_statistic_df['t_BAF_diff_of_leftside_with_interval'] = t_BAF_diff_of_leftside_with_interval
    RP_SR_depth_BAF_statistic_df['t_BAF_diff_of_rightside_with_interval'] = t_BAF_diff_of_rightside_with_interval
    RP_SR_depth_BAF_statistic_df['n_BAF_diff_of_leftside_with_interval'] = n_BAF_diff_of_leftside_with_interval
    RP_SR_depth_BAF_statistic_df['n_BAF_diff_of_rightside_with_interval'] = n_BAF_diff_of_rightside_with_interval

    RP_SR_depth_BAF_statistic_df['valid_t_RP_percent_left_bp_interval'] = valid_t_RP_percent_left_bp_interval
    RP_SR_depth_BAF_statistic_df['valid_t_RP_percent_right_bp_interval'] = valid_t_RP_percent_right_bp_interval
    RP_SR_depth_BAF_statistic_df['valid_n_RP_percent_left_bp_interval'] = valid_n_RP_percent_left_bp_interval
    RP_SR_depth_BAF_statistic_df['valid_n_RP_percent_right_bp_interval'] = valid_n_RP_percent_right_bp_interval
    RP_SR_depth_BAF_statistic_df['valid_t_SR_percent_left_bp_interval'] = valid_t_SR_percent_left_bp_interval
    RP_SR_depth_BAF_statistic_df['valid_t_SR_percent_right_bp_interval'] = valid_t_SR_percent_right_bp_interval
    RP_SR_depth_BAF_statistic_df['valid_n_SR_percent_left_bp_interval'] = valid_n_SR_percent_left_bp_interval
    RP_SR_depth_BAF_statistic_df['valid_n_SR_percent_right_bp_interval'] = valid_n_SR_percent_right_bp_interval
    print(RP_SR_depth_BAF_statistic_df.shape)
    RP_SR_depth_BAF_statistic_df.to_csv(out_statistic_csv, index=False)

def generate_proof(tumor,normal,ref,bed,sample_name,RP_dir,RD_dir,mean_insert_size,std_insert_size,snp_vcf_prefix,thread,out_statistic_csv):
    PR_SR_nineth(bed, RP_dir, normal, tumor,mean_insert_size,std_insert_size,sample_name,thread)
    depth_BAF_nineth(bed, RD_dir, normal, tumor,ref,snp_vcf_prefix,sample_name,thread)
    statisticproof(RP_dir, RD_dir, out_statistic_csv)



