# Variant Calling and Annotation Example

Due to technical difficulties (Macbook M1), I couldn't recreate a previous iteration of this repo. However, since I want to keep a version of this on my GitHub, I will instead write a guide for Variant Calling and Annotation. This guide will not be as deep as a 
university course but will provide a great outline on what we are doing & why we are doing it. 

---

## Step 0: Required Software / Packages

###### QC & Trimming:
- FastQC
- MultiQC
- Trimmomatic

###### Alignment:
- BWA-MEM
- samtools

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

### Step 3: Alignment

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
bwa mem -t 8 Homo_sapiens.GRCh38.dna.primary_assembly.fa output_foward_paired.fastq.gz \
 | samtools view -bS -o sample.bam
```

Finally, we will sort and index the BAM file.

```bash
samtools sort sample.bam -o sample_sorted.bam
samtools index sample_sorted.bam
```

<!-- ### Step 4: Variant Calling

Here we will use  -->