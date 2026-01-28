import argparse
from ./src/extract_feature import generate_proof

def main():
    parser = argparse.ArgumentParser(description='Extract feature for training or evaluation')
    parser.add_argument('--tumor_bam', type=str, required=True, help='Path to tumor BAM')
    parser.add_argument('--normal_bam', type=str, required=True, help='Path to normal BAM')
    parser.add_argument('--fasta', type=str, required=True,default='./model.pth', help='Path to reference genome')
    parser.add_argument('--bed', type=str,required=True, default='./model.pth', help='Path to CNV bed,need four columns, "ch1\tstart\tend\tDUP\nchr1\tstart\tend\tDEL\n"')
    parser.add_argument('--sample_name', required=True,type=str, default='./model.pth', help='Name of sample')
    parser.add_argument('--RP_dir', type=str, required=True,default='./model.pth', help='Path to store read-pair split-read related information')
    parser.add_argument('--RD_dir', type=str, required=True,default='./model.pth', help='Path to store read-depth BAF related information')
    parser.add_argument('--mean_insert_size', type=int, required=True,default=500, help='Mean insert size of tumor or normal BAM')
    parser.add_argument('--std_insert_size', type=int, required=True,default=50, help='Standard deviation of insert size of tumor or normal BAM')
    parser.add_argument('--snp_vcf_prefix', type=str, required=True,default='./model.pth', help='Prefix of SNP vcf called by Varscan,"/xxx/xxx/HG002.ILM.varscan."')
    parser.add_argument('--thread_num', type=int, required=True,default=1, help='Num of thread using for extracting feature')
    parser.add_argument('--output', type=str, required=True,default='./model.pth', help='Output CSV dir')

    args = parser.parse_args()
    generate_proof(args.tumor_bam,args.normal_bam,args.fasta,args.bed,args.sample_name,args.RP_dir,args.RD_dir,args.mean_insert_size,args.std_insert_size,args.snp_vcf_prefix,args.thread,args.output)

if __name__ == '__main__':
    main()