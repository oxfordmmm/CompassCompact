#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    =========================================
                   MMM Compass Compact 
    =========================================
    Usage:

    nextflow run main_bwa.nf --help
    nextflow run main_bwa.nf --test -profile test_docker
    nextflow run main_bwa.nf \
    --input_dir tests/data/input_dir/ \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_3.fasta \
    --mask_file tests/data/reference/NC_000962_3_repmask.array \
    --fastq true \
    --pattern "*_{1,2}.fastq.gz" \
    -profile test_docker

    nextflow run main_bwa.nf \
    --input_dir tests/data/input_dir/ \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_3.fasta \
    --mask_file tests/data/reference/NC_000962_3_repmask.array \
    --fastq false \
    --pattern "*.bam" \
    -profile test_docker

    Mandatory arguments:
    --input_dir               DIR       path of fastq files, or bam files
    --fastq                   Boolean   Input files are fastq format
    --pattern                 String    fastq file name pattern, such as "*_{1,2}.fastq.gz"
    --output_dir              DIR       path for transformed fastq files
    --ref                     FILE      reference genome
    --mask_file               FILE      reference mask array

    Optional arguments:
    --threads                 INT       number of threads to run, default 4

    """.stripIndent()
}

params.help = false;
params.test = false;

if (params.test){
    params.fastq = true
    params.pattern = "*_{1,2}.fastq.gz"
    params.input_dir = "tests/data/input_dir/"
    params.output_dir = "tests/data/output_dir"
    params.threads = 8
    params.ref = "tests/data/reference/NC_000962_3.fasta"
    params.mask_file = "tests/data/reference/NC_000962_3_repmask.array"
}else{
    params.fastq = false
    params.pattern = ""
    params.input_dir = ""
    params.output_dir = ""
    params.threads = 8
    params.ref = ""
    params.mask_file = ""
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

ref = file(params.ref)
ref_folder = ref.getParent()
ref_name = ref.getBaseName()

masked_ref = file(params.mask_file)
masked_ref_folder = masked_ref.getParent()
masked_ref_name = masked_ref.getBaseName()

threads = params.threads
data_path = params.input_dir + params.pattern

if (params.fastq == false){
    bam_files_channel = Channel.fromPath(data_path)

    process bam2fastq {
        scratch true
        echo true
        //publishDir "${params.output_dir}/bam2fastq", mode: "copy"

        tag {bam_file.getBaseName()}

        input:
        file bam_file from bam_files_channel

        output:
        set val("${bam_file.getBaseName()}"), file("${bam_file.getBaseName()}_1.fastq"), file("${bam_file.getBaseName()}_2.fastq") into bwa_read_pairs

        """
        ${SAMTOOL1x}/samtools bam2fq \
        -1 ${bam_file.getBaseName()}_1.fastq \
        -2 ${bam_file.getBaseName()}_2.fastq \
        ${bam_file.getBaseName()}.bam
        """
    }
}else{
    Channel
    .fromFilePairs( data_path , flat: true) //Find reads from a given path
    .ifEmpty{ error "Cannot find any reads matching: ${params.pattern}"}
    .set { bwa_read_pairs }
}

process bwa_index {
    memory '1 GB'

    echo true
    scratch true

    publishDir "${params.output_dir}/indexed_ref", mode: "copy"

    tag {ref}

    input:
        file ref

    output:
        file "*" into bwa_index

    script:
    """
    $BWA/bwa index ${ref}
    """
}
// Map reads to reference genome with BWA MEM
process bwa{
    memory '4 GB'

    echo true
    scratch true

    //publishDir "${params.output_dir}/bwa_out", mode: "copy", pattern: "*.sam"

    tag {dataset_id}

    input:
    set dataset_id, file(forward), file(reverse) from bwa_read_pairs
    file index from bwa_index
    file ref

    output:
    set dataset_id, file("${dataset_id}_alignment.sam") into bwa_map

    """
    $BWA/bwa mem -R '@RG\tID:${dataset_id}\tSM:null\tLB:null\tCN:null' \
    -t ${threads} \
    ${ref} ${forward} ${reverse} \
    > ${dataset_id}_alignment.sam
    """
}

// Merging SAMfiles from different readgroup mappings
process bwa_merge {
    memory '13 GB'

    echo true
    scratch true

    publishDir "${params.output_dir}/bwa_merge", mode: "copy", pattern: "${dataset_id}*"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}_alignment.sam") from bwa_map

    output:
    set dataset_id, file("${dataset_id}_alignment.bam"), file("${dataset_id}_seqstats.txt"), file("${dataset_id}_flagstats.txt") into bwa_merge

    """
    python $COMPASS_ROOT/nf_bwa_merge.py \
    -b ${dataset_id}_alignment.sam \
    -o ${dataset_id}_alignment.bam \
    -execute \
    -ss ${dataset_id}_seqstats.txt \
    -fs ${dataset_id}_flagstats.txt
    """
}

process mpileup {
    memory '12 GB'

    echo true
    scratch true

    //publishDir "${params.output_dir}/mpileup", mode: "copy", pattern: "${dataset_id}*"

    tag {dataset_id}

    input:
    set dataset_id,  file("${dataset_id}_alignment.bam") from bwa_merge
    file ref

    output:
    set dataset_id, file("${dataset_id}.out.vcf"), file("${dataset_id}.pileup.vcf") into snpcalling

    """
    python $COMPASS_ROOT/nf_mpileup.py \
    -o 40 -e 20 -H 100 -m 2 -F 0.002 -D -S -M0 -q 30 \
    -Q25 -E -c -g -K -L -t0.01 -i -1 -p0.5 -P full \
    -B ${dataset_id}_alignment.bam \
    -R ${ref} \
    -out ${dataset_id}.out.vcf \
    -outpileup ${dataset_id}.pileup.vcf
    """
}


// Add extra information for VCF file
process annotvcf {
    memory '2 GB'

    echo true
    scratch true

    publishDir "${params.output_dir}/annotvcf", mode: "copy", pattern: "${dataset_id}*"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}.out.vcf"), file("${dataset_id}.pileup.vcf") from snpcalling
    file masked_ref
    //file "${ref_name}_repmask.array" from mask_ref


    output:
    set dataset_id, file("${dataset_id}.annotvcf.vcf") into annotcvf

    """
    python $COMPASS_ROOT/nf_annotvcf.py \
    -vcf ${dataset_id}.out.vcf \
    -mpileup ${dataset_id}.pileup.vcf \
    -refmask ${masked_ref} \
    -o ${dataset_id}.annotvcf.vcf \
    -basecall -selfblastR -hqdepthinfo -lgcdepthinfo
    """
}


process basecall {
    memory '1 GB'

    echo true
    scratch true

    publishDir "${params.output_dir}/basecall", mode: "move", pattern: "${dataset_id}*"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}.annotvcf.vcf") from annotcvf
    output:
    set dataset_id, file("*") into basecall

    """
    python $COMPASS_ROOT/nf_basecall.py \
    -A 25 -e DISABLED -E DISABLED -g DISABLED -G DISABLED -K0.90 -J DISABLED \
    -invcf ${dataset_id}.annotvcf.vcf \
    -outvcf ${dataset_id}.basecall.vcf.gz \
    -outvcfIndel ${dataset_id}.basecall_Indel.vcf.gz \
    -outfasta ${dataset_id}.consensus.fasta.gz \
    -outstats ${dataset_id}.basecallstats.txt \
    -u ${dataset_id} \
    -refid ${ref_name} \
    -Q30 -q30 -m30 -n5 -S25 -I25 -z -B1 -p -N -f0.35
    """
}
