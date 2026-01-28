configfile: "./config/config.yaml"  
samtools = config["samtools"]
varscan = config["varscan"]
tumor = config['tumor']
normal = config['normal']
snp_bed = config["snp_bed"]
dir_out = config["dir_out_varscan"]
dir_data =config["dir_out_samtools"]
path_ref = config["path_ref"]
sample_name = config['sample']
wildcard_constraints:
    chrom ="|".join(['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
                'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                'chr21', 'chr22'])

targets = []
for chr in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
                'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                'chr21', 'chr22']:
    targets.append(dir_out + f"{sample_name}.varscan.{chr}.snv.vcf.snp")


rule all:
    input:
        targets

rule samtools_n_mpileup:
    input:
        bam = normal
    output:
        dir_data + "{sample}.N.{chrom}.samtools.pileup"
    run:
        shell("{samtools} mpileup -r {wildcards.chrom} -f {path_ref} -q 1  {input.bam} -o {output}")
rule samtools_t_mpileup:
    input:
        bam = tumor
    output:
        dir_data + "{sample}.T.{chrom}.samtools.pileup"
    run:
        shell("{samtools} mpileup -r {wildcards.chrom} -f {path_ref} -q 1  {input.bam} -o {output}")

rule varscan_call_snp_indel:
    input:
        normal_pileup=dir_data + "{sample}.N.{chrom}.samtools.pileup",
        tumor_pileup=dir_data + "{sample}.T.{chrom}.samtools.pileup",
    output:
        dir_out + "{sample}.varscan.{chrom}.snv.vcf.snp",
    run:
        shell("{varscan} somatic {input.normal_pileup} {input.tumor_pileup} {dir_out}/{wildcards.sample}.varscan.{wildcards.chrom}.snv.vcf")
