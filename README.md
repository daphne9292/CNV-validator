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
Python              3.9.0
Snakemake           7.25.0
pandas              1.5.1
sklearn             1.2.1
scipy               1.7.3
```


## Usage

CNV-validator is composed of three modules: feature extraction, model training, and evaluation of independent CNV callsets.
For a given CNV callset, features can be extracted from the paired tumor–normal BAM files and used as input to either a newly trained model or the provided pre-trained model for CNV validation.

### Feature extraction

Because the B-allele frequency (BAF) deviation of heterozygous SNPs are used in feature extraction, SNP calling must be performed on both normal and tumor samples firstly. In this work, SNP calling was carried out using VarScan. The detailed commands are provided in varscan.smk, and users are required to modify the corresponding config.yaml file in the config directory before running the pipeline. The usage of the .smk file is as follows:

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
--expand_size       Expand size for RP feature extraction
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
We provide a pre-trained model in the model directory, which can be applied to the evaluation of external CNV callsets.

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
--expand_size           Expand size for RP feature extraction
--snp_vcf_prefix        Prefix of SNP vcf called by Varscan (eg:"/xxx/xxx/HG002.ILM.varscan.")
--thread_num            Num of thread using for extracting feature
--output_feature        Output extrated feature, CSV format.

validation parameters
--model                 Path to model
--output_model_predict  Output model validation result,CSV format.

```

The detailed commands used is listed below.

```python
python validator.py --tumor_bam --normal_bam --fasta --bed --sample_name --RP_dir --RD_dir --mean_insert_size --std_insert_size --snp_vcf_prefix --thread_num --output_feature --model --output_model_predict
```

#### Trained model

We trained a CatBoost classifier using both simulated and real datasets. The model classifies CNV candidates into six categories. Class 0 represents regions without detectable CNV signals. Class 1 represents CNVs with accurately positioned left and right breakpoints. Classes 2 and 3 represent potential CNVs with accurately positioned left and right breakpoints, respectively. Class 4 represents intervals formed by combining breakpoints from two distinct events. Class 5 represents regions with potential loss of heterozygosity (LOH), but with both breakpoints inaccurately positioned.

Users can select different prediction categories according to their analytical requirements. For applications requiring high-confidence CNVs with accurate intervals and breakpoints, only Class 1 predictions should be retained, although this may exclude some true CNV events. In contrast, when maximizing CNV recall is the priority, Classes 2 and 3 can also be considered to retain more potential CNV events.

A pre-trained model is provided in the "model" directory.You can use the model directly to validate your CNV callset.

#### Demo

The `demo` directory provides a small example dataset in `demo/data/` and a test script, `demo/demo.test.sh`, to demonstrate feature extraction and CNV callset validation using the pre-trained model. Before running the demo, check the input paths and parameters in the script and ensure that the required dependencies are installed. Run the following command from the project root:

```bash
bash demo/demo.test.sh
```

### CNV detection Software

 We list several widely used software for CNV detection.The callset of these tools were used for evaluation the performance of CNV-validator.Command-line parameters are listed below.

#### FACETS

Linux

```python
snp-pileup -A -g -q 10 -Q 10 -r 5,5 -P 100 SNP.benchmark.vcf.gz output_pileup.csv.gz normal.bam tumor.bam
```

R

```python
library(facets)
rcmat <- readSnpMatrix('output_pileup.csv.gz')
x <-preProcSample(rcmat, ndepth=10, het.thresh=0.25, snp.nbhd=100, cval=25)
y <- procSample(x, cval = 50, min.nhet = 5)
fit <- emcncf(y)
write.table(fit$cncf, file = "facets_segments.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

#### CNVkit

```python
cnvkit.py batch tumor.bam --normal normal.bam --method wgs --fasta ref.fasta --output-dir ./ -p 48 --access access.hg38.bed
```

#### Suquenza

python

```python
sequenza-utils gc_wiggle -w 50 -f ref.fasta -o hg38.gc50Base.wig.gz
sequenza-utils bam2seqz -n normal.bam -t tumor.bam --fasta ref.fasta -gc hg38.gc50Base.wig.gz -o out.seqz.gz -S ~/miniconda3/bin/samtools --parallel 22 -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 -T ~/miniconda3/bin/tabix
sequenza-utils seqz_binning -s chr1_22.seqz.gz -w 50 -o small.chr1_22.seqz.gz -T ~/miniconda3/bin/tabix
```

R

```python
library('sequenza')
test<-sequenza.extract('small.chr1_22.seqz.gz')
cp<-sequenza.fit(test)
sequenza.results(sequenza.extract = test,cp.table=cp,sample.id = 'tumor',out.dir='./')

```
