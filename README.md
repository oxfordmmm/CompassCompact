# CompassCompact
A compact version of Oxford Compass (Complete Pathogen Analytical Software Solution)


## Overview

A nextflow-docker paired pipeline for processing pathogen bacterial sequencing data generated using Illumina sequencing platform. It runs following:

    1. Take paired fastq files or bam files in a file folder
    2. Map reads to reference genome using Stampy (main_stampy.nf) or BWA (main_bwa.nf)
    3. SNP calling using samtools and bcftools, and 
    4. Annotate VCF using masked reference and create a consensus sequence fasta file

### Input
    input files directory
    fastq or bam files pattern 
    reference genome

### Output
    *.basecallstats.txt
    *.consensus.fasta.gz
    *.basecall_indel.vcf.gz
    *.basecall.vcf.gz

## Variant Call with Stampy

### 1. Get Docker Image from Docker Hub
```bash
    docker pull oxfordmmm/compasscompact:{version}
```
### 2. Run nextflow with docker images
```bash
    nextflow run main_stampy.nf --help
    nextflow run main_stampy.nf --test -profile test_docker
    nextflow run main_stampy.nf \
    --input_dir tests/data/input_dir/ \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_2.fasta \
    --fastq true \
    --pattern "*_{1,2}.fastq.gz" \
    -profile test_docker

    nextflow run main_stampy.nf \
    --input_dir tests/data/input_dir/ \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_2.fasta \
    --fastq false \
    --pattern "*.bam" \
    -profile test_docker
```
#### Mandatory arguments:
```bash
    --input_dir               DIR       path of fastq files, or bam files
    --fastq                   Boolean   Input files are fastq format
    --pattern                 String    fastq file name pattern, such as "*_{1,2}.fastq.gz"
    --output_dir              DIR       path for transformed fastq files
    --ref                     FILE      reference genome
```
#### Optional arguments:
```bash
    --threads                 INT       number of threads to run, default 4
```
#### Run Validation Tests with Stampy

    1. Copy a pair input fastq files or a bam file to tests/data/test_input
    2. Copy genomo reference fasta as tests/data/reference/NC_000962_3.fasta
    3. Copy expected basecall output fasta to tests/data/test_output/expected_output
```bash
    python3 tests/test_stampy.py bam (under CompassCompact, test bam input)
    python3 tests/test_stampy.py fastq (under CompassCompact, test fastq input)
```
    The test will run the stampy nextflow pipeline and compare the output fasta file with expected fasta file.

## Variant Call with BWA

### 1. Get Docker Image from Docker Hub
```bash
    docker pull oxfordmmm/compasscompact:{version}
```
### 2. Generate reference mask array file (will be used in step 3)
```bash
    nextflow run mask_ref.nf --help
    nextflow run mask_ref.nf --test -profile test_docker
    nextflow run mask_ref.nf \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_2.fasta \
    --mask true \
    -profile test_docker
```
#### Mandatory arguments:
```bash
    --ref                    FILE         reference genome
```
#### Optional arguments:
```bash
    --mask                   Boolean      use self-blast to mask repeated region, default true
```
### 3. Run nextflow with docker images
```bash
    nextflow run main_bwa.nf --help
    nextflow run main_bwa.nf --test
    nextflow run main_bwa.nf \
    --input_dir tests/data/input_dir/ \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_3.fasta \
    --mask_file "tests/data/reference/NC_000962_3_repmask.array" \
    --fastq true \
    --pattern "*_{1,2}.fastq.gz" \
    -profile test_docker

    nextflow run main_bwa.nf --input_dir tests/data/input_dir/ \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_3.fasta \
    --mask_file = "tests/data/reference/NC_000962_3_repmask.array" \
    --fastq false \
    --pattern "*.bam" \
    -profile test_docker

```
#### Mandatory arguments:
```bash
    --input_dir               DIR       path of fastq files, or bam files
    --fastq                   Boolean   Input files are fastq format
    --pattern                 String    fastq file name pattern, such as "*_{1,2}.fastq.gz"
    --output_dir              DIR       path for transformed fastq files
    --ref                     FILE      reference genome
    --mask_file               FILE      reference mask array
```
#### Optional arguments:
```bash
    --threads                 INT       number of threads to run, default 4
```
#### Run Validation Tests with BWA

    1. Copy a pair input fastq files or a bam file to tests/data/test_input
    2. Copy genomo reference fasta as tests/data/reference/NC_000962_3.fasta
    3. Copy genomo reference mask array as tests/data/reference/NC_000962_3_repmask.array
    4. Copy expected basecall output fasta to tests/data/test_output/expected_output
```bash
    python3 tests/test_bwa.py bam (under CompassCompact, test bam input)
    python3 tests/test_bwa.py fastq (under CompassCompact, test fastq input)
```
    The test will run the bwa nextflow pipeline and compare the output fasta file with expected fasta file.

### Bioinformatic tools

    bwa-0.7.15	
    https://github.com/lh3/bwa

    GenomeAnalysisTK-3.7-0	
    https://github.com/broadinstitute/gatk/ 

    ncbi-blast-2.2.23+	
    https://blast.ncbi.nlm.nih.gov/Blast.cgi

    picard-tools-1.123	
    https://github.com/broadinstitute/picard/

    stampy-1.0.23	
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3106326/

    samtools-1.4.1	
    https://github.com/samtools/samtools/releases

    vcftools_0.1.9	
    https://vcftools.github.io/downloads.html

    All bioinformatic tools are configured in  `docker/compass/lib/compass.cfg`

    To debug with different version of tools or parameters, change `nextflow.config` the volume host path 
    from    `/home/docker/Code/CompassCompact/docker/compass` 
    to wherever the compass code directory is (typically, where you clone the reponsitory to + `/docker/compass`).

