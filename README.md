# CNV-validator

CNV-validator is a deep learning–based model trained on both simulated and real sequencing data.CNV-validator evaluates candidate CNV sets by assessing both the CNV intervals and the accuracy of their breakpoints, providing a robust framework for CNV validation.

## Requirements

### Operation systems

Linux-based systems required, including but are not limited to:

- MacOS, Big Sur
- Ubuntu (including Windows Subsystem)
- CentOS

### Software and packages

```python
GPU Driver Version  455.23.05 
CUDA Version        11.1
Python              3.9.0
Snakemake           7.25.0
pandas              1.5.1
sklearn             1.2.1
scipy               1.7.3
```


## Usage

---
CNV-validator is composed of three modules: feature extraction, model training, and evaluation of independent CNV callsets.
For a given CNV callset, features can be extracted from the paired tumor–normal BAM files and used as input to either a newly trained model or the provided pre-trained model for CNV validation.

### Feature extraction

Because the B-allele frequency (BAF) of heterozygous SNPs are used in feature extraction, SNP calling must be performed on both normal and tumor samples firstly. In this work, SNP calling was carried out using VarScan. The detailed commands are provided in varscan.smk, and users are required to modify the corresponding config.yaml file in the config directory before running the pipeline. The usage of the .smk file is as follows:

```snakemake
nohup snakemake -s ./varscan.smk --cluster "qsub -l nodes=1:ppn=1 " -j 1 >sublog 2>&1 &
```

#### Required Input/Ouput parameters

```python
--tumor_bam         Path to tumor BAM
--normal_bam        Path to normal BAM
--fasta             Path to reference genome
--bed               Path to CNV bed (Need four columns. eg: "ch1\tstart\tend\tDUP\nchr1\tstart\tend\tDEL\n")
--sample_name       Name of sample
--RP_dir            Path to store read-pair split-read related information
--RD_dir            Path to store read-depth BAF related information
--mean_insert_size  Mean insert size of tumor or normal BAM
--std_insert_size   Standard deviation of insert size of tumor or normal BAM
--snp_vcf_prefix    Prefix of SNP vcf called by Varscan (eg:"/xxx/xxx/HG002.ILM.varscan.")
--thread_num        Num of thread using for extracting feature
--output            Output CSV dir
```

The detailed commands used is listed below.

```python
python extract.py --tumor_bam --normal_bam --fasta --bed --sample_name --RP_dir --RD_dir --mean_insert_size --std_insert_size --snp_vcf_prefix --thread_num --output
```

### CNV callset validation

For a CNV callset produced by an external method, users can directly apply the validator command to assess its reliability using the trained model.

```python
feature extraction parameters
--tumor_bam             Path to tumor BAM
--normal_bam            Path to normal BAM
--fasta                 Path to reference genome
--bed                   Path to CNV bed (Need four columns. eg: "ch1\tstart\tend\tDUP\nchr1\tstart\tend\tDEL\n")
--sample_name           Name of sample
--RP_dir                Path to store read-pair split-read related information
--RD_dir                Path to store read-depth BAF related information
--mean_insert_size      Mean insert size of tumor or normal BAM
--std_insert_size       Standard deviation of insert size of tumor or normal BAM
--snp_vcf_prefix        Prefix of SNP vcf called by Varscan (eg:"/xxx/xxx/HG002.ILM.varscan.")
--thread_num            Num of thread using for extracting feature
--output_feature        Output extrated feature, CSV format.

validation parameters
--device                Devices using for model validation.eg:"cpu","gpu1".
--model                 Path to model
--output_model_predict  Output model validation result,CSV format.

```

The detailed commands used is listed below.

```python
python validator.py --tumor_bam --normal_bam --fasta --bed --sample_name --RP_dir --RD_dir --mean_insert_size --std_insert_size --snp_vcf_prefix --thread_num --output_feature --device --model --output_model_predict
```