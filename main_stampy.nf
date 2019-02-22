#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    =========================================
                   MMM Compass Next 
    =========================================
    Usage:

    nextflow run main_stampy.nf --help
    nextflow run main_stampy.nf --test
    nextflow run main_stampy.nf --input_dir tests/data/input_dir/fastqs/ \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_2.fasta \
    --fastq true
    --pattern "*_{1,2}P.fastq.gz"

    nextflow run main_stampy.nf --input_dir tests/data/input_dir/bams/ \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_2.fasta \
    --fastq false \
    --pattern "*.bam"

    Mandatory arguments:
    --input_dir               DIR       path of fastq files, or bam files
    --fastq                   Boolean   Input files are fastq format
    --pattern                 String    fastq file name pattern, such as "*_{1,2}P.fastq.gz"
    --output_dir              DIR       path for transformed fastq files
    --ref                     FILE      reference genome

    Optional arguments:
    --threads                 INT       number of threads to run, default 4

    """.stripIndent()
}

params.help = false;
params.test = false;

if (params.test){
    params.fastq = true
    params.pattern = "*_{1,2}.fastq.gz"
    params.input_dir = "tests/data/input_dir/fastqs/"
    params.output_dir = "tests/data/output_dir"
    params.threads = 8
    params.ref = "tests/data/reference/NC_000962_2.fasta"
}else{
    params.fastq = false
    params.pattern = ""
    params.input_dir = ""
    params.output_dir = ""
    params.threads = 8
    params.ref = ""
}

ref = file(params.ref)
ref_name = ref.getBaseName()

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

threads = params.threads
data_path = params.input_dir + params.pattern

if (params.fastq){
    Channel
    .fromFilePairs( data_path , flat: true) //Find reads from a given path
    .ifEmpty{ error "Cannot find any reads matching: ${params.pattern}"}
    .set { read_pairs }
    // Convert fastq to unmapped bam
    process fastq2bam {
        
        scratch true
        echo true
        tag {dataset_id}

        // publishDir "${params.output_dir}/bams", mode: "copy"

        input:
        set dataset_id, file(forward), file(reverse) from read_pairs

        output:
        file("${dataset_id}.bam") into bam_files

        """
        java -jar ${PICARD}/FastqToSam.jar \\
        F1=${forward} \\
        F2=${reverse} \\
        OUTPUT=${dataset_id}.bam.tmp \\
        READ_GROUP_NAME=${dataset_id} \\
        SAMPLE_NAME=${dataset_id} \\
        LIBRARY_NAME=unknown \\
        PLATFORM=Illumina \\
        SEQUENCING_CENTER=unknown \\
        RUN_DATE=null

        samtools view -H ${dataset_id}.bam.tmp  |\
        sed -e '2i\\@SQ\\tSN\\:unmapped\\tLN\\:0' |\
        samtools reheader - ${dataset_id}.bam.tmp > ${dataset_id}.bam

        rm ${dataset_id}.bam.tmp
        """
        }
    }
else{
     Channel
    .fromPath(data_path)
    .ifEmpty{ error "Cannot find any bam files"} //If not found, return error
    .into{ bam_files }
}

process GenerateMaskReference {
    memory '4 GB'

    scratch true
    echo true
    tag {ref}

    publishDir "${params.output_dir}/ref", mode: "copy" , pattern: "${ref.getBaseName()}*"

    input:
    file ref

    output:
    file("${ref_name}_repmask.array") into masked_ref
    set file(ref), file("${ref_name}.stidx"), file("${ref_name}.sthash") into (stampy_index, mpileup_index)
    """
    python ${COMPASS_ROOT}/nf_ref_index.py -r ${ref}
    """
}

// Map reads to reference genome with BWA MEM
process Stampy {
    memory '10 GB'

    scratch true
    echo true

    tag {"${bam.getBaseName()}"}

    //publishDir "${params.output_dir}/stampy", mode: "copy", pattern: "${bam.getBaseName()}*"

    input:
    file(bam) from bam_files
    set file(ref), file("${ref_name}.stidx"), file("${ref_name}.sthash") from stampy_index

    output:
    set val("${bam.getBaseName()}"), file("${bam.getBaseName()}_alignment_stampy.bam"), file("${bam.getBaseName()}_seqstats.txt"), file("${bam.getBaseName()}_flagstats.txt") into stampy_map

    """
    python $COMPASS_ROOT/nf_stampy.py \
    -b $bam \
    -r $ref_name \
    -o  ${bam.getBaseName()}_alignment_stampy.bam \
    -ss ${bam.getBaseName()}_seqstats.txt \
    -fs ${bam.getBaseName()}_flagstats.txt
    """
}

// SNP calling
process Mpileup {
    // publishDir "${params.out_dir}/Mpileup2nd", mode: "copy", pattern: "${dataset_id}*"

    memory '12 GB'
    scratch true
    echo true
    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}_alignment_stampy.bam") from stampy_map
    set file(ref), file("${ref_name}.stidx"), file("${ref_name}.sthash") from mpileup_index

    output:
    set dataset_id, file("${dataset_id}.out_stampy.vcf"), file("${dataset_id}.pileup_stampy.vcf") into snpcalling

    """
    python $COMPASS_ROOT/nf_mpileup.py \
    -o 40 -e 20 -H 100 -m 2 -F 0.002 -D -S -M0 -q 30 \
    -Q25 -E -c -g -K -L -t0.01 -i -1 -p0.5 -P full \
    -B ${dataset_id}_alignment_stampy.bam \
    -R $ref \
    -out ${dataset_id}.out_stampy.vcf \
    -outpileup ${dataset_id}.pileup_stampy.vcf

    if [-d ../tmp]; then
        rm -rf ../tmp
    fi
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
