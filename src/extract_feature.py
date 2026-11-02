#coding:utf-8
import pandas as pd
import pysam
import numpy as np
import re
import math
import os
import multiprocessing
import random
import subprocess
import gzip
chrom_dict = {'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9,'chr10': 10, \
                  'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17,'chr18': 18, \
                  'chr19': 19, 'chr20': 20, 'chr21': 21, 'chr22': 22, 'chrX': 23, 'chrY': 24}

def parse_cigar(cigar_str):
    """
    Parse a CIGAR string into tuples of (operation, length).

    Parameters:
    cigar_str (str): The CIGAR string to parse.
    """
    # 使用正则表达式拆分CIGAR字符串
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    cigar_operations = cigar_pattern.findall(cigar_str)

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

    step =  max(1, int(len(depth)/20))
    for i in range(0,len(depth),step):
        num = num + 1
        sum_depth += depth[i]
    if num > 0:
        ave_depth = sum_depth / num
    return ave_depth

def cluster_pos(chr_pos_list, dist_thresh=500):
    """
    chr_pos_list: [(chr:str, pos:int), ...]
    """

    if not chr_pos_list:
        return (0, 0)
    from collections import defaultdict

    chr_dict = defaultdict(list)
    for c,p in chr_pos_list:
        chr_dict[c].append(p)

    all_clusters = []
    for _, pos_arr in chr_dict.items():
        pos_sorted = sorted(pos_arr)
        cur_cluster = [pos_sorted[0]]
        for p in pos_sorted[1:]:
            if p - cur_cluster[-1] <= dist_thresh:
                cur_cluster.append(p)
            else:
                all_clusters.append(cur_cluster)
                cur_cluster = [p]
        all_clusters.append(cur_cluster)

    # total_num = sum(len(c) for c in all_clusters)
    cluster_count = len(all_clusters)
    max_cluster_size = max(len(c) for c in all_clusters) if all_clusters else 0

    return cluster_count, max_cluster_size

def cal(abnormal_rp,sa_pos,total_read_of_window,the_other_bp,neighbor_size,chrom):
    rp_cnt_l = len(abnormal_rp)
    sa_cnt_l = len(sa_pos)

    f1 = rp_cnt_l / total_read_of_window if total_read_of_window else 0.0
    f2 = sa_cnt_l / total_read_of_window if total_read_of_window else 0.0

    rp_pos_l = [(c,p) for c,p,s,ins in abnormal_rp ]
    f3, f4 = cluster_pos(rp_pos_l)
    f5 = f4/total_read_of_window if total_read_of_window >0 else 0
    opp_s = the_other_bp - neighbor_size
    opp_e = the_other_bp + neighbor_size
    in_opp_l = sum(1 for c,p,s,ins in abnormal_rp if c==chrom and opp_s <= p <= opp_e)
    f6 = in_opp_l / rp_cnt_l if rp_cnt_l else 0.0 # proportion of reads which MATE pos is in the other bp neighbor

    same_l = sum(1 for c,p,s,ins in abnormal_rp if s)# proportion of reads which direction is same with MATE
    f7 = same_l / rp_cnt_l if rp_cnt_l else 0.0

    zero_ins_l = sum(1 for c,p,s,ins in abnormal_rp if ins == 0)
    f8 = zero_ins_l / rp_cnt_l if rp_cnt_l else 0.0 # proportion of reads which insert size is 0

    f9, f10 = cluster_pos(sa_pos)
    f11 = f10/total_read_of_window if total_read_of_window >0 else 0
    sa_in_opp_l = sum(1 for c,p in sa_pos if c==chrom and opp_s <= p <= opp_e)
    f12 = sa_in_opp_l / sa_cnt_l if sa_cnt_l else 0.0
    left = [total_read_of_window,rp_cnt_l,f1,f3,f4,f5,f6,f7,f8,sa_cnt_l,f2,f9,f10,f11,f12]
    return left

def fetch_reads_tenth(bam_file, chrom, start, end, type,mean_insert_size =413,std_insert_size = 106,expand_size = 500):
    bam = pysam.AlignmentFile(bam_file, "rb")
    abnormal_size = mean_insert_size + 2 * std_insert_size if (end - start) < 1500 else 1000

    left_win_s = max(0, start - expand_size)
    left_win_e = start + expand_size
    right_win_s = max(0, end - expand_size)
    right_win_e = end + expand_size

    opp_left_break = end
    opp_right_break = start
    opp_win = expand_size

    left_total_read = 0
    left_abnormal_rp = []  # (chr, pos, is_same_strand, ins_val)
    left_sa_pos = []  # (chr, pos)
    total_lbp_read = set()

    k = 0
    # mate() seeks on bam; keep the window iterator on a separate handle.
    for read in bam.fetch(chrom, left_win_s, left_win_e, multiple_iterators=True):
        k += 1
        if k > 4000: break
        if read.mapping_quality < 20: continue
        if not read.is_paired or read.mate_is_unmapped: continue
        try:
            mate = bam.mate(read)
        except ValueError as exc:
            # Region-filtered BAMs can retain paired flags but omit the mate.
            if str(exc) != "mate not found":
                raise
            continue
        if not mate: continue
        if read.query_name not in total_lbp_read:
            total_lbp_read.add(read.query_name)
            left_total_read += 1
        else:
            continue

        # SR SA
        if ("S" in read.cigarstring or "H" in read.cigarstring) and read.has_tag("SA"):
            sa = read.get_tag("SA")
            multi_sa_alignment = sa.split(';')[0:-1]
            if len(multi_sa_alignment) > 2:
                continue
            for seg in multi_sa_alignment:
                items = seg.split(",")
                supplementary_mapping_quality = int(items[4])
                if supplementary_mapping_quality >= 30:
                    left_sa_pos.append((items[0], int(items[1])))

        # CIGAR
        cig_r = parse_cigar(str(read.cigarstring))
        cig_m = parse_cigar(str(mate.cigarstring))
        maxM_r = max([int(x[0]) for x in cig_r if x[1] == "M"], default=0)
        maxM_m = max([int(x[0]) for x in cig_m if x[1] == "M"], default=0)
        cond_m = (maxM_r / read.query_length > 0.3) and (maxM_m / mate.query_length > 0.3)

        ins_abs = abs(read.template_length)
        if (ins_abs > abnormal_size or ins_abs == 0)  and cond_m:
            if type == 'DUP':
                if (not read.is_reverse):
                    continue
            elif type == 'DEL':
                if read.is_reverse:
                    continue

            r1_plus = not read.is_reverse
            r2_plus = not mate.is_reverse
            rp_same_dir = (r1_plus and r2_plus) or ((not r1_plus) and (not r2_plus))
            left_abnormal_rp.append((mate.reference_name, mate.reference_start, rp_same_dir, ins_abs))

    right_total_read = 0
    right_abnormal_rp = []
    right_sa_pos = []
    total_rbp_read =set()

    k = 0
    for read in bam.fetch(chrom, right_win_s, right_win_e, multiple_iterators=True):
        k += 1
        if k > 4000: break
        if read.mapping_quality < 20: continue
        if not read.is_paired or read.mate_is_unmapped: continue
        try:
            mate = bam.mate(read)
        except ValueError as exc:
            if str(exc) != "mate not found":
                raise
            continue
        if not mate: continue
        if read.query_name not in total_rbp_read:
            total_rbp_read.add(read.query_name)
            right_total_read += 1
        else:
            continue

        if ("S" in read.cigarstring or "H" in read.cigarstring) and read.has_tag("SA"):
            sa = read.get_tag("SA")
            multi_sa_alignment = sa.split(';')[0:-1]
            if len(multi_sa_alignment) > 2:
                continue
            for seg in multi_sa_alignment:
                items = seg.split(",")
                supplementary_mapping_quality = int(items[4])
                if supplementary_mapping_quality >= 30:
                    right_sa_pos.append((items[0], int(items[1])))

        cig_r = parse_cigar(str(read.cigarstring))
        cig_m = parse_cigar(str(mate.cigarstring))
        maxM_r = max([int(x[0]) for x in cig_r if x[1] == "M"], default=0)
        maxM_m = max([int(x[0]) for x in cig_m if x[1] == "M"], default=0)
        cond_m = (maxM_r / read.query_length > 0.3) and (maxM_m / mate.query_length > 0.3)

        ins_abs = abs(read.template_length)
        if (ins_abs > abnormal_size or ins_abs == 0) and cond_m:
            if type == 'DUP':
                if read.is_reverse:
                    continue
            elif type == 'DEL':
                if not read.is_reverse:
                    continue
            r1_plus = not read.is_reverse
            r2_plus = not mate.is_reverse
            rp_same_dir = (r1_plus and r2_plus) or ((not r1_plus) and (not r2_plus))
            right_abnormal_rp.append((mate.reference_name, mate.reference_start, rp_same_dir, ins_abs))
    bam.close()

    left_f = cal(left_abnormal_rp, left_sa_pos, left_total_read, opp_left_break, opp_win, chrom)
    right_f = cal(right_abnormal_rp,right_sa_pos, right_total_read, opp_right_break, opp_win, chrom)
    return left_f + right_f
def fetch_snp_baf(fasta,chrom, start, end,snp_vcf,normal_bam_dir,tumor_bam_dir):
    read_fa = pysam.FastaFile(fasta)
    chrom_length = read_fa.get_reference_length(chrom)
    svlen = end - start
    if svlen < 10000:
        expand_size = 3000
    elif svlen < 100000:
        expand_size = min(int(svlen * 0.3), 8000)
    else:
        expand_size = 8000
    left = max(1,start - expand_size)
    right = min(chrom_length,end + expand_size)
    total_line = wc_count(snp_vcf)
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

    ave_tumor_baf_dev_leftside_of_interval = 0
    ave_normal_baf_dev_leftside_of_interval = 0
    ave_tumor_baf_dev_of_interval = 0
    ave_normal_baf_dev_of_interval = 0
    ave_tumor_baf_dev_rightside_of_interval = 0
    ave_normal_baf_dev_rightside_of_interval = 0
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
                    ave_tumor_baf_dev_leftside_of_interval += abs(tumor_baf-0.5)
                    ave_normal_baf_dev_leftside_of_interval += abs(normal_baf-0.5)
                    num4 += 1
            elif int(detail[1]) <= end:
                if abs(normal_baf - 0.5) < 0.15:
                    ave_tumor_baf_dev_of_interval += abs(tumor_baf-0.5)
                    ave_normal_baf_dev_of_interval += abs(normal_baf-0.5)
                    num5 += 1
            elif int(detail[1]) <= right:
                if abs(normal_baf - 0.5) < 0.15:
                    ave_tumor_baf_dev_rightside_of_interval += abs(tumor_baf-0.5)
                    ave_normal_baf_dev_rightside_of_interval += abs(normal_baf-0.5)
                    num6 += 1

    ave_t_depth_leftside_of_interval = get_samtools_depth(chrom,left,start - 1,tumor_bam_dir)
    ave_n_depth_leftside_of_interval = get_samtools_depth(chrom, left, start - 1, normal_bam_dir)
    ave_t_depth_of_interval = get_samtools_depth(chrom,start,end,tumor_bam_dir)
    ave_n_depth_of_interval = get_samtools_depth(chrom,start,end,normal_bam_dir)
    ave_t_depth_rightside_of_interval = get_samtools_depth(chrom,end+1,right,tumor_bam_dir)
    ave_n_depth_rightside_of_interval = get_samtools_depth(chrom,end+1,right,normal_bam_dir)

    t_depth_ratio_of_leftside_with_interval = ave_t_depth_of_interval / ave_t_depth_leftside_of_interval if ave_t_depth_leftside_of_interval > 0 else np.nan
    t_depth_ratio_of_rightside_with_interval = ave_t_depth_of_interval / ave_t_depth_rightside_of_interval if ave_t_depth_rightside_of_interval > 0 else np.nan
    n_depth_ratio_of_leftside_with_interval = ave_n_depth_of_interval / ave_n_depth_leftside_of_interval if ave_n_depth_leftside_of_interval > 0 else np.nan
    n_depth_ratio_of_rightside_with_interval = ave_n_depth_of_interval / ave_n_depth_rightside_of_interval if ave_n_depth_rightside_of_interval > 0 else np.nan

    if num4 > 0 :
        ave_tumor_baf_dev_leftside_of_interval = ave_tumor_baf_dev_leftside_of_interval/num4
        ave_normal_baf_dev_leftside_of_interval = ave_normal_baf_dev_leftside_of_interval/num4
    else:
        ave_tumor_baf_dev_leftside_of_interval = np.nan
        ave_normal_baf_dev_leftside_of_interval = np.nan
    if num5 > 0:
        ave_tumor_baf_dev_of_interval = ave_tumor_baf_dev_of_interval / num5
        ave_normal_baf_dev_of_interval=ave_normal_baf_dev_of_interval/num5
    else:
        ave_tumor_baf_dev_of_interval = np.nan
        ave_normal_baf_dev_of_interval = np.nan
    if num6 > 0:
        ave_tumor_baf_dev_rightside_of_interval = ave_tumor_baf_dev_rightside_of_interval / num6
        ave_normal_baf_dev_rightside_of_interval = ave_normal_baf_dev_rightside_of_interval/num6
    else:
        ave_tumor_baf_dev_rightside_of_interval= np.nan
        ave_normal_baf_dev_rightside_of_interval = np.nan

    return t_depth_ratio_of_leftside_with_interval,t_depth_ratio_of_rightside_with_interval,\
              n_depth_ratio_of_leftside_with_interval,n_depth_ratio_of_rightside_with_interval,\
              ave_tumor_baf_dev_leftside_of_interval,ave_normal_baf_dev_leftside_of_interval,\
              ave_tumor_baf_dev_of_interval,ave_normal_baf_dev_of_interval,\
            ave_tumor_baf_dev_rightside_of_interval,ave_normal_baf_dev_rightside_of_interval

def chunker_RP_SR_nineth(list,RP_SR_out_dir,n_bam,t_bam,mean_insert_size,std_insert_size,sample_name='',expand_size=500):
    for i in range(len(list)):
        chrom = list[i][0]
        start = int(list[i][1])
        end = int(list[i][2])
        type = list[i][3]

        sum_read_of_lbp_w,ab_rp_num_lbp_w,nfl1,nfl2,nfl3,nfl4,nfl5,nfl6,nfl7,sa_num_lbp_w,nfl8,nfl9,nfl10,nfl11,nfl12,\
            sum_read_of_rbp_w,ab_rp_num_rbp_w,nfr1,nfr2,nfr3,nfr4,nfr5,nfr6,nfr7,sa_num_rbp_w,nfr8,nfr9,nfr10,nfr11,nfr12 = fetch_reads_tenth(n_bam, chrom, start, end,type,mean_insert_size ,std_insert_size ,expand_size)
        t_sum_read_of_lbp_w,t_ab_rp_num_lbp_w,tfl1, tfl2, tfl3, tfl4, tfl5, tfl6, tfl7, t_sa_num_lbp_w,tfl8, tfl9, tfl10, tfl11, tfl12,\
            t_sum_read_of_rbp_w,t_ab_rp_num_rbp_w,tfr1, tfr2, tfr3, tfr4, tfr5, tfr6, tfr7,t_sa_num_rbp_w, tfr8, tfr9, tfr10, tfr11, tfr12 = fetch_reads_tenth(t_bam, chrom, start, end,type,mean_insert_size,std_insert_size ,expand_size)
        out_dir = RP_SR_out_dir + sample_name + '_'+chrom+'_'+str(start)+'_'+str(end)+'_'+list[i][3]+'.proof.txt'
        if sample_name == '':
            out_dir = RP_SR_out_dir + chrom + '_' + str(start) + '_' + str(end) + '_' + list[i][3] + '.proof.txt'

        with open(out_dir, "w") as fout:
            print(sum_read_of_lbp_w,ab_rp_num_lbp_w,nfl1,nfl2,nfl3,nfl4,nfl5,nfl6,nfl7,sa_num_lbp_w,nfl8,nfl9,nfl10,nfl11,nfl12,\
                  sum_read_of_rbp_w,ab_rp_num_rbp_w,nfr1,nfr2,nfr3,nfr4,nfr5,nfr6,nfr7,sa_num_rbp_w,nfr8,nfr9,nfr10,nfr11,nfr12,\
                t_sum_read_of_lbp_w,t_ab_rp_num_lbp_w,tfl1,tfl2,tfl3,tfl4,tfl5,tfl6,tfl7,t_sa_num_lbp_w,tfl8,tfl9,tfl10,tfl11,tfl12,\
                  t_sum_read_of_rbp_w,t_ab_rp_num_rbp_w,tfr1,tfr2,tfr3,tfr4,tfr5,tfr6,tfr7,t_sa_num_rbp_w,tfr8,tfr9,tfr10,tfr11,tfr12, file=fout)

def chunker_depth_BAF_nineth(list,depth_BAF_out_dir,normal_bam_dir,tumor_bam_dir,fasta,snp_vcf_prefix,sample_name=''):
    for i in range(len(list)):
        chrom = list[i][0]
        start = int(list[i][1])
        end = int(list[i][2])
        type = list[i][3]
        snp_vcf = snp_vcf_prefix + chrom + ".snv.vcf.snp"
        out_dir = depth_BAF_out_dir + sample_name + '_'+ chrom + '_' + str(start) + '_' + str(end) + '_' + type + '.proof.txt'
        if sample_name == '':
            out_dir = depth_BAF_out_dir + chrom + '_' + str(start) + '_' + str(end) + '_' + type + '.proof.txt'
        t_depth_ratio_of_leftside_with_interval, t_depth_ratio_of_rightside_with_interval, \
            n_depth_ratio_of_leftside_with_interval, n_depth_ratio_of_rightside_with_interval, \
            ave_tumor_baf_dev_leftside_of_interval, ave_normal_baf_dev_leftside_of_interval, \
            ave_tumor_baf_dev_of_interval, ave_normal_baf_dev_of_interval, \
            ave_tumor_baf_dev_rightside_of_interval, ave_normal_baf_dev_rightside_of_interval = fetch_snp_baf(fasta, chrom, start, end, snp_vcf,normal_bam_dir,tumor_bam_dir)
        with open(out_dir, "w") as fout:
            print(t_depth_ratio_of_leftside_with_interval,t_depth_ratio_of_rightside_with_interval,\
              n_depth_ratio_of_leftside_with_interval,n_depth_ratio_of_rightside_with_interval,\
              ave_tumor_baf_dev_leftside_of_interval,ave_normal_baf_dev_leftside_of_interval,\
              ave_tumor_baf_dev_of_interval,ave_normal_baf_dev_of_interval,\
            ave_tumor_baf_dev_rightside_of_interval,ave_normal_baf_dev_rightside_of_interval,file=fout)

def PR_SR_nineth(benchmark_bed,RP_SR_dir,normal_bam,tumor_bam,mean_insert_size,std_insert_size,sample_name='',expand_size=500,thread_num=1):
    bed_list = []
    with open(benchmark_bed, 'r') as fin:
        for line in fin:
            info = line.strip().split()

            sta = int(info[1])
            end = int(info[2])
            svlen = end - sta
            if len(info) == 4:
                type = info[3]
                if type !='DUP' and type !='DEL':
                    continue
                proof = RP_SR_dir + sample_name + '_'+ info[0] + '_' + info[1] + '_' + info[2] + '_' + type + '.proof.txt'
                if sample_name == '':
                    proof =  RP_SR_dir +info[0] + '_' + info[1] + '_' + info[2] + '_' + type + '.proof.txt'
                if svlen > 1000 and not os.path.exists(proof):
                    bed_list.append([info[0], info[1], info[2], type])

    total_row = len(bed_list)
    chunk_size = total_row // thread_num + 1
    jobs = []

    for i in range(thread_num):
        list_chunk = [row for row in bed_list[i * chunk_size: min((i + 1) * chunk_size, total_row)]]
        p = multiprocessing.Process(target=chunker_RP_SR_nineth, args=(list_chunk, RP_SR_dir,normal_bam,tumor_bam,mean_insert_size,std_insert_size,sample_name,expand_size))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

def depth_BAF_nineth(benchmark_bed,depth_BAF_dir,normal_bam_dir,tumor_bam_dir,fasta,snp_vcf_prefix,sample_name,thread_num):
    bed_list =[]
    with open(benchmark_bed,'r') as fin:
        for line in fin:
            info = line.strip().split()
            sta = int(info[1])
            end = int(info[2])
            svlen = end - sta
            if len(info) == 4:
                type = info[3]
                if type !='DUP' and type !='DEL':
                    continue
                proof = depth_BAF_dir + sample_name + '_'+  info[0] + '_' + info[1] + '_' + info[2] + '_' + type + '.proof.txt'
                if sample_name == '':
                    proof =  depth_BAF_dir +info[0] + '_' + info[1] + '_' + info[2] + '_' + type + '.proof.txt'
                if svlen > 1000 and not os.path.exists(proof):
                    bed_list.append([info[0], info[1], info[2], type])

    total_row=len(bed_list)
    chunk_size = total_row // thread_num + 1
    jobs = []

    for i in range(thread_num):
        list_chunk = [row for row in bed_list[i * chunk_size: min((i + 1) * chunk_size,total_row)]]
        p = multiprocessing.Process(target=chunker_depth_BAF_nineth, args=(list_chunk,depth_BAF_dir,normal_bam_dir,tumor_bam_dir,fasta,snp_vcf_prefix,sample_name))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

def statisticproof(RP_SR_proof_dir,depth_BAF_proof_dir,out_statistic_csv):
    RP_SR_depth_BAF_statistic_df = pd.DataFrame()
    id = []
    report_len = []
    t_depth_ratio_of_leftside_with_interval = []
    t_depth_ratio_of_rightside_with_interval = []
    n_depth_ratio_of_leftside_with_interval = []
    n_depth_ratio_of_rightside_with_interval = []

    valid_t_RP_percent_left_bp_interval = []
    RP_mate_loci_cluster_num_t_lbp = []
    maxC_RPnum_divide_totalRP_t_lbp = []
    mate_in_rbp_neigb_divide_abRP_t_lbp =[]
    valid_t_SR_percent_left_bp_interval = []
    SR_SA_cluster_num_t_lbp =[]
    maxC_SRnum_divide_totalR_t_lbp = []
    SA_in_rbp_neigb_divide_SR_t_lbp = []

    valid_t_RP_percent_right_bp_interval = []
    RP_mate_loci_cluster_num_t_rbp = []
    maxC_RPnum_divide_totalRP_t_rbp = []
    mate_in_lbp_neigb_divide_abRP_t_rbp = []
    valid_t_SR_percent_right_bp_interval = []
    SR_SA_cluster_num_t_rbp = []
    maxC_SRnum_divide_totalR_t_rbp = []
    SA_in_lbp_neigb_divide_SR_t_rbp = []

    valid_n_RP_percent_left_bp_interval = []
    RP_mate_loci_cluster_num_n_lbp = []
    maxC_RPnum_divide_totalRP_n_lbp = []
    mate_in_rbp_neigb_divide_abRP_n_lbp = []
    SR_SA_cluster_num_n_lbp = []
    maxC_SRnum_divide_totalR_n_lbp = []
    SA_in_rbp_neigb_divide_SR_n_lbp = []
    valid_n_SR_percent_left_bp_interval = []

    valid_n_RP_percent_right_bp_interval = []
    RP_mate_loci_cluster_num_n_rbp = []
    maxC_RPnum_divide_totalRP_n_rbp = []
    mate_in_lbp_neigb_divide_abRP_n_rbp = []
    SR_SA_cluster_num_n_rbp = []
    maxC_SRnum_divide_totalR_n_rbp = []
    SA_in_lbp_neigb_divide_SR_n_rbp = []
    valid_n_SR_percent_right_bp_interval = []

    ave_tumor_baf_dev_leftside_of_interval = []
    ave_normal_baf_dev_leftside_of_interval = []
    ave_tumor_baf_dev_of_interval = []
    ave_normal_baf_dev_of_interval = []
    ave_tumor_baf_dev_rightside_of_interval = []
    ave_normal_baf_dev_rightside_of_interval = []
    sv_type = []

    RP_SR_dir = RP_SR_proof_dir
    depth_BAF_dir = depth_BAF_proof_dir
    files = os.listdir(depth_BAF_dir)

    for file in files:
        proof_id = file[:-len('.proof.txt')]
        filename = proof_id.rsplit('_', 4)
        sample = filename[0] if len(filename) == 5 else ''
        chrom, start, end, type = filename[-4:]
        start = int(start)
        end = int(end)
        svlen = end - start
        RP_SR_file_path = RP_SR_dir + file

        if os.path.isfile(RP_SR_file_path):
            with open(RP_SR_file_path, 'r') as fin:
                for line in fin:
                    info = line.split()
                    valid_n_RP_percent_left_bp_interval.append(float(info[2]))
                    RP_mate_loci_cluster_num_n_lbp.append(int(info[3]))
                    maxC_RPnum_divide_totalRP_n_lbp.append(float(info[5]))
                    mate_in_rbp_neigb_divide_abRP_n_lbp.append(float(info[6]))
                    SR_SA_cluster_num_n_lbp.append(int(info[11]))
                    maxC_SRnum_divide_totalR_n_lbp.append(float(info[13]))
                    SA_in_rbp_neigb_divide_SR_n_lbp.append(float(info[14]))
                    valid_n_SR_percent_left_bp_interval.append(float(info[10]))

                    valid_n_RP_percent_right_bp_interval.append(float(info[17]))
                    RP_mate_loci_cluster_num_n_rbp.append(int(info[18]))
                    maxC_RPnum_divide_totalRP_n_rbp.append(float(info[20]))
                    mate_in_lbp_neigb_divide_abRP_n_rbp.append(float(info[21]))
                    SR_SA_cluster_num_n_rbp.append(int(info[26]))
                    maxC_SRnum_divide_totalR_n_rbp.append(float(info[28]))
                    SA_in_lbp_neigb_divide_SR_n_rbp.append(float(info[29]))
                    valid_n_SR_percent_right_bp_interval.append(float(info[25]))

                    valid_t_RP_percent_left_bp_interval.append(float(info[32]))
                    RP_mate_loci_cluster_num_t_lbp.append(int(info[33]))
                    maxC_RPnum_divide_totalRP_t_lbp.append(float(info[35]))
                    mate_in_rbp_neigb_divide_abRP_t_lbp.append(float(info[36]))
                    SR_SA_cluster_num_t_lbp.append(int(info[41]))
                    maxC_SRnum_divide_totalR_t_lbp.append(float(info[43]))
                    SA_in_rbp_neigb_divide_SR_t_lbp.append(float(info[44]))
                    valid_t_SR_percent_left_bp_interval.append(float(info[40]))

                    valid_t_RP_percent_right_bp_interval.append(float(info[47]))
                    RP_mate_loci_cluster_num_t_rbp.append(int(info[48]))
                    maxC_RPnum_divide_totalRP_t_rbp.append(float(info[50]))
                    mate_in_lbp_neigb_divide_abRP_t_rbp.append(float(info[51]))
                    SR_SA_cluster_num_t_rbp.append(int(info[56]))
                    maxC_SRnum_divide_totalR_t_rbp.append(float(info[58]))
                    SA_in_lbp_neigb_divide_SR_t_rbp.append(float(info[59]))
                    valid_t_SR_percent_right_bp_interval.append(float(info[55]))

            id.append(proof_id)
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
                    ave_tumor_baf_dev_leftside_of_interval.append(float(info[4]))
                    ave_normal_baf_dev_leftside_of_interval.append(float(info[5]))
                    ave_tumor_baf_dev_of_interval.append(float(info[6]))
                    ave_normal_baf_dev_of_interval.append(float(info[7]))
                    ave_tumor_baf_dev_rightside_of_interval.append(float(info[8]))
                    ave_normal_baf_dev_rightside_of_interval.append(float(info[9]))

    RP_SR_depth_BAF_statistic_df['id'] = id
    RP_SR_depth_BAF_statistic_df['report_len'] = report_len
    RP_SR_depth_BAF_statistic_df['svtype'] = sv_type
    RP_SR_depth_BAF_statistic_df['ave_tumor_baf_dev_leftside_of_interval'] = ave_tumor_baf_dev_leftside_of_interval
    RP_SR_depth_BAF_statistic_df['ave_tumor_baf_dev_of_interval'] = ave_tumor_baf_dev_of_interval
    RP_SR_depth_BAF_statistic_df['ave_tumor_baf_dev_rightside_of_interval'] = ave_tumor_baf_dev_rightside_of_interval
    RP_SR_depth_BAF_statistic_df['t_depth_ratio_of_leftside_with_interval'] = t_depth_ratio_of_leftside_with_interval
    RP_SR_depth_BAF_statistic_df['t_depth_ratio_of_rightside_with_interval'] = t_depth_ratio_of_rightside_with_interval
    RP_SR_depth_BAF_statistic_df['valid_t_RP_percent_left_bp_interval'] = valid_t_RP_percent_left_bp_interval
    RP_SR_depth_BAF_statistic_df['RP_mate_loci_cluster_num_t_lbp'] = RP_mate_loci_cluster_num_t_lbp
    RP_SR_depth_BAF_statistic_df['maxC_RPnum_divide_totalRP_t_lbp'] = maxC_RPnum_divide_totalRP_t_lbp
    RP_SR_depth_BAF_statistic_df['mate_in_rbp_neigb_divide_abRP_t_lbp'] = mate_in_rbp_neigb_divide_abRP_t_lbp
    RP_SR_depth_BAF_statistic_df['valid_t_SR_percent_left_bp_interval'] = valid_t_SR_percent_left_bp_interval
    RP_SR_depth_BAF_statistic_df['SR_SA_cluster_num_t_lbp'] = SR_SA_cluster_num_t_lbp
    RP_SR_depth_BAF_statistic_df['maxC_SRnum_divide_totalR_t_lbp'] = maxC_SRnum_divide_totalR_t_lbp
    RP_SR_depth_BAF_statistic_df['SA_in_rbp_neigb_divide_SR_t_lbp'] = SA_in_rbp_neigb_divide_SR_t_lbp

    RP_SR_depth_BAF_statistic_df['valid_t_RP_percent_right_bp_interval'] = valid_t_RP_percent_right_bp_interval
    RP_SR_depth_BAF_statistic_df['RP_mate_loci_cluster_num_t_rbp'] = RP_mate_loci_cluster_num_t_rbp
    RP_SR_depth_BAF_statistic_df['maxC_RPnum_divide_totalRP_t_rbp'] = maxC_RPnum_divide_totalRP_t_rbp
    RP_SR_depth_BAF_statistic_df['mate_in_lbp_neigb_divide_abRP_t_rbp'] = mate_in_lbp_neigb_divide_abRP_t_rbp
    RP_SR_depth_BAF_statistic_df['valid_t_SR_percent_right_bp_interval'] = valid_t_SR_percent_right_bp_interval
    RP_SR_depth_BAF_statistic_df['SR_SA_cluster_num_t_lbp'] = SR_SA_cluster_num_t_lbp
    RP_SR_depth_BAF_statistic_df['maxC_SRnum_divide_totalR_t_rbp'] = maxC_SRnum_divide_totalR_t_rbp
    RP_SR_depth_BAF_statistic_df['SA_in_lbp_neigb_divide_SR_t_rbp'] = SA_in_lbp_neigb_divide_SR_t_rbp


    RP_SR_depth_BAF_statistic_df['ave_normal_baf_dev_leftside_of_interval'] = ave_normal_baf_dev_leftside_of_interval
    RP_SR_depth_BAF_statistic_df['ave_normal_baf_dev_of_interval'] = ave_normal_baf_dev_of_interval
    RP_SR_depth_BAF_statistic_df['ave_normal_baf_dev_rightside_of_interval'] = ave_normal_baf_dev_rightside_of_interval
    RP_SR_depth_BAF_statistic_df['n_depth_ratio_of_leftside_with_interval'] = n_depth_ratio_of_leftside_with_interval
    RP_SR_depth_BAF_statistic_df['n_depth_ratio_of_rightside_with_interval'] = n_depth_ratio_of_rightside_with_interval
    RP_SR_depth_BAF_statistic_df['valid_n_RP_percent_left_bp_interval'] = valid_n_RP_percent_left_bp_interval
    RP_SR_depth_BAF_statistic_df['RP_mate_loci_cluster_num_n_lbp'] = RP_mate_loci_cluster_num_n_lbp
    RP_SR_depth_BAF_statistic_df['maxC_RPnum_divide_totalRP_n_lbp'] = maxC_RPnum_divide_totalRP_n_lbp
    RP_SR_depth_BAF_statistic_df['mate_in_rbp_neigb_divide_abRP_n_lbp'] = mate_in_rbp_neigb_divide_abRP_n_lbp
    RP_SR_depth_BAF_statistic_df['valid_n_SR_percent_left_bp_interval'] = valid_n_SR_percent_left_bp_interval
    RP_SR_depth_BAF_statistic_df['SR_SA_cluster_num_n_lbp'] = SR_SA_cluster_num_n_lbp
    RP_SR_depth_BAF_statistic_df['maxC_SRnum_divide_totalR_n_lbp'] = maxC_SRnum_divide_totalR_n_lbp
    RP_SR_depth_BAF_statistic_df['SA_in_rbp_neigb_divide_SR_n_lbp']= SA_in_rbp_neigb_divide_SR_n_lbp

    RP_SR_depth_BAF_statistic_df['valid_n_RP_percent_right_bp_interval'] = valid_n_RP_percent_right_bp_interval
    RP_SR_depth_BAF_statistic_df['RP_mate_loci_cluster_num_n_rbp'] = RP_mate_loci_cluster_num_n_rbp
    RP_SR_depth_BAF_statistic_df['maxC_RPnum_divide_totalRP_n_rbp'] = maxC_RPnum_divide_totalRP_n_rbp
    RP_SR_depth_BAF_statistic_df['mate_in_lbp_neigb_divide_abRP_n_rbp'] = mate_in_lbp_neigb_divide_abRP_n_rbp
    RP_SR_depth_BAF_statistic_df['valid_n_SR_percent_right_bp_interval'] = valid_n_SR_percent_right_bp_interval
    RP_SR_depth_BAF_statistic_df['SR_SA_cluster_num_n_lbp'] = SR_SA_cluster_num_n_rbp
    RP_SR_depth_BAF_statistic_df['maxC_SRnum_divide_totalR_n_rbp'] = maxC_SRnum_divide_totalR_n_rbp
    RP_SR_depth_BAF_statistic_df['SA_in_lbp_neigb_divide_SR_n_rbp'] = SA_in_lbp_neigb_divide_SR_n_rbp
    RP_SR_depth_BAF_statistic_df.to_csv(out_statistic_csv, index=False)

def generate_proof(tumor,normal,ref,bed,sample_name,RP_dir,RD_dir,mean_insert_size,std_insert_size,expand_size,snp_vcf_prefix,thread,out_statistic_csv):
    PR_SR_nineth(bed, RP_dir, normal, tumor,mean_insert_size,std_insert_size,sample_name,expand_size,thread)
    depth_BAF_nineth(bed, RD_dir, normal, tumor, ref,snp_vcf_prefix,sample_name,thread)
    statisticproof(RP_dir, RD_dir, out_statistic_csv)
