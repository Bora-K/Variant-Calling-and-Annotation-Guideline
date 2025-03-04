# Variant Calling and Annotation Example

Due to technical difficulties (Macbook M1), I couldn't recreate a previous iteration of this repo. However, since I want to keep a version of this on my GitHub, I will instead write a guide for Variant Calling and Annotation. This guide will not be as deep as a 
university course but will provide a great outline on what we are doing & why we are doing it. 

---
<!--TODO: Add lists under every code block that explains each argument-->
## Step 0: Required Software / Packages

###### QC & Trimming:
- FastQC
- MultiQC
- Trimmomatic

###### Alignment:
- BWA-MEM
- samtools

###### Duplicate Removal:
- Picard

###### Variant Calling:
- GATK

###### Annotation:
- VEP

## Step 1: Acquiring Data

For this tutorial case we will use a whole genome sequencing (WGS) dataset. However, this method can be applied to whole exome sequencing (WES), targeted sequencing, etc. 

We will be using publicly available data from ["Cytidine deaminase role on replicative stress in pancreatic cancer cells"](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982645/#:~:text=Cytidine%20deaminase%20reduces%20replication%20stress,that%20could%20enhance%20treatment%20response.) study. The runs are available under [SRR27639115](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR27639115&display=metadata) accession code.

To get the runs in the fastq format we will use the following commands from the **STA Toolkit**:

- prefetch: For getting the sra file
- fasterq-dump: For extracting the raw reads


```bash
prefetch SRR27639115
fasterq-dump --split-files SRR27639115
gzip SRR27639115_1.fastq SRR27639115_2.fastq
```
After this step, we will have two fastq files: one for reverse reads and one for forward reads. This is because these reads are **paired-end**. 

## Step 2: Quality Control (QC)

We will be using **FastQC** in this step. This will determine the sequencing quality of our reads. **MultiQC** will allow us to combine the QC reports.

```bash
fastqc SRR27639115_1.fastq.gz SRR27639115_2.fastq.gz
multiqc .
```

Then we will use **Trimmomatic** to get rid of the low-quality reads and adapters.

```bash
trimmomatic PE -phred33 SRR27639115_1.fastq.gz SRR27639115_2.fastq.gz output_forward_paired.fastq.gz output_forward_unpaired.fastq.gz output_reverse_paired.fastq.gz output_reverse_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:2 MINLEN:36
```

To verify the success of the trimming, you can check the post-trimming quality using:

```bash
fastqc output_foward_paired.fastq.gz output_reverse_paired.fastq.gz
```

## Step 3: Alignment

To align our reads we need a reference genome. Since our organism is *Homo sapiens*, we can easily use the publicly available genome.

```bash
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

We will use the **Burrows-Wheeler Aligner (BWA)**, the most commonly used aligner. First, we need to index the reference genome. 

```bash
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Then we can align the reads using BWA-MEM, an algoirthm of BWA designed for short reads. We will also use **samtools** to convert the SAM file to a BAM file.

```bash
bwa mem -t 8 Homo_sapiens.GRCh38.dna.primary_assembly.fa output_foward_paired.fastq.gz\
output_reverse_paired.fastq.gz | samtools view -bS -o sample.bam
```

Finally, we will sort and index the BAM file.

```bash
samtools sort sample.bam -o sample_sorted.bam
samtools index sample_sorted.bam
```

We can remove the duplicates here using **Picard**.
```bash
picard MarkDuplicates I=sample_sorted.bam O=sample_dedup.bam M=duplication_metrics.txt REMOVE_DUPLICATES=true
```

## Step 4: Variant Calling

**GATK** works in 3 steps. First, the **BaseRecalibrator** tool analyzes the base quality scores and creates a recalibration table. This table contains information about the systematic errors observed in the data.

```bash
gatk BaseRecalibrator -I sample_dedup.bam -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --known-sites dbsnp.vcf -O recal_data.table
```

The **ApplyBQSR** tool uses the recalibration table to adjust the base quality scores in the BAM file.

```bash
gatk ApplyBQSR -R Homo_sapiens.GRCh38.dna.primary_assembly.fa -I sample_dedup.bam \
  --bqsr-recal-file recal_data.table -O sample_recal.bam
```

The **HaplotypeCaller** tool in GATK is used to call variants (both SNPs and indels) in the recalibrated BAM file.

```bash
gatk HaplotypeCaller -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -I sample_recal.bam -O sample_raw_variants.vcf -ERC GVCF
```

For somatic variant calling **Mutect2** is used, which is specifically designed for detecting somatic mutations. This is a grea tool to use in cancer samples. 

```bash
gatk Mutect2 -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -I sample_recal.bam -tumor MiaPaCa2 -O somatic_variants.vcf
```

## Step 5: Annotation

Now that we have the raw variants we need to annotate them. We will use **Variant Effect Predictor (VEP)** for this purpose.

```bash
vep --cache --dir /path/to/vep/cache --assembly GRCh38 \
  --input_file sample_raw_variants.vcf --output_file sample_annotated.vcf --format vcf --everything
```
<!--TODO: We can add alternative annotation tools here-->

## Step 6: Filtering Variants

We can filter variants based on quality score and depth. We can keep the high-confidence variants. 

```bash
bcftools filter -i 'QUAL>30 && DP>10' sample_annotated.vcf > sample_filtered.vcf
```

This is it for variant calling and annotation. 

<!--TODO: Add a downstream analysis as well-->