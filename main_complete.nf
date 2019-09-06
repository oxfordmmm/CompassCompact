indir = params.indir
readpat = params.readpat
output_dir = params.output_dir
kraken2_db = file(params.kraken2_db).toAbsolutePath()
ref = params.ref
input_filetype = params.input_filetype
refmap = params.refmap


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

//
// If inputs are in bam format, convert to fastq.gz
//
if (input_filetype == "bam") {
    channel_bam_files = Channel.fromPath(indir + readpat)

    process bam2fastq {
	label "clockwork"

	tag { bam_file.getBaseName() }

	input:
	file bam_file from channel_bam_files

	output:
	set val("${bam_file.getBaseName()}"), file("${bam_file.getBaseName()}_1.fastq.gz"), file("${bam_file.getBaseName()}_2.fastq.gz") into channel_fastqs

	"""
        samtools bam2fq -1 ${bam_file.getBaseName()}_1.fastq \
                        -2 ${bam_file.getBaseName()}_2.fastq \
                           ${bam_file}
        gzip ${bam_file.getBaseName()}_1.fastq || true
        gzip ${bam_file.getBaseName()}_2.fastq || true
        """
    }
}

//
// If files are in uncompressed fastq format, compress them
//
if (input_filetype == "fastq") {
    Channel.fromFilePairs(indir + readpat).ifEmpty("No files found").set{ channel_fastqs1 }

    //
    // used to turn the input format from list to whatever?
    //
    process compress_fastq {
	label "clockwork"

	tag { dataset_id }

	input:
	set dataset_id, reads from channel_fastqs1

	output:
	set dataset_id, file("${dataset_id}.fastq.gz"), file("${dataset_id}.fastq.gz") into channel_fastqs

	"""
        gzip -c ${reads[0]} > ${dataset_id}.fastq.gz || true
        gzip -c ${reads[1]} > ${dataset_id}.fastq.gz || true
        """
    }
}

//
// If input file are in fastq.gz format, don't do anything
//
if (input_filetype == "fastq.gz") {
    Channel.fromFilePairs(indir + readpat).map{ [it[0],it[1][0],it[1][1]] }.set{ channel_fastqs }
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

process kraken2_get_genus {
    label "fatos"

    publishDir "${output_dir}/${dataset_id}/speciation", mode: "copy"

    memory '10 GB'

    tag { dataset_id }

    input:
    set dataset_id, read1, read2 from channel_fastqs

    output:
    set dataset_id, read1, read2, stdout into input_contamremoval
    set dataset_id, read1, read2 into input_classification
    file("${dataset_id}_kraken2.tab")

    script:
    """
    /usr/bin/kraken2 --db ${kraken2_db} --output - --report "${dataset_id}_kraken2.tab" --paired "${read1}" "${read2}" 1>/dev/null

    KRAKEN=`python3 /usr/bin/get_genus.py "${dataset_id}_kraken2.tab"`
    curl -s --header "Content-Type: application/json" --request POST --data "{\\"pipeline_name\\": \\"${params.pipeline_name}\\", \\"run_uuid\\":\\"${params.run_uuid}\\", \\"sample_name\\": \\"${dataset_id}\\", \\"tag_type\\": \\"kraken2\\", \\"tag_name\\": \\"\$KRAKEN\\" }" http://${params.head_node_ip}:12000/add_sample_tag > /dev/null

    printf `python3 /usr/bin/get_genus.py "${dataset_id}_kraken2.tab"`
    """
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

process classification {
    label "fatos"

    tag { run_id }

    publishDir "${output_dir}/${run_id}/classification", mode: 'copy'

    input:
    set val(run_id), read1, read2 from input_classification

    output:
    set val(run_id), file("${run_id}.classification_non_human_read_list.txt") into non_human_list

    script:
    kraken2_summary              = "${run_id}.species_classification.txt"
    kraken2_read_classification  = "${run_id}.read_classification.txt"
    kraken2_human_read_list      = "${run_id}.classification_human_read_list.txt"
    kraken2_non_human_read_list  = "${run_id}.classification_non_human_read_list.txt"
    """
    cat <<- EOF >> filter.py
    import sys

    filename = sys.argv[1]

    tbl1 = open(filename).readlines()
    tbl1 = [x.split('\t') for x in tbl1]

    #
    # add missing indexes to sample reads 1-9 (bug in kraken2)
    # for whatever reason, in the kraken2 output the first 10 lines
    # don't have the .N ending. So if the 11th line has .12 we
    # add them back
    #
    if tbl1[11][1][:-3] == '.12':
        for i, row in enumerate(tbl1[0:9]):
            row[1] = f'{row[1]}.{i+1}'

    #
    # write non-human read ids
    #
    human_id = '9606'

    for row in tbl1:
        if row[2] != human_id:
            print(row[1])
    EOF

    kraken2 --threads ${task.cpus} --db ${kraken2_db} --report ${kraken2_summary} --output ${kraken2_read_classification} --paired ${read1} ${read2}

    echo "==== kraken2 ====" > ${kraken2_human_read_list}
    cat ${kraken2_summary} | grep 9606 >> ${kraken2_human_read_list}
    echo "==== human reads ====" >> ${kraken2_human_read_list}
    awk '\$3==\"9606\" { print \$2 }' ${kraken2_read_classification} >> ${kraken2_human_read_list}
    python3 filter.py ${kraken2_read_classification} > ${kraken2_non_human_read_list}
    """
}

non_human_list.join(input_contamremoval).set { input_contamremoval }

process contam_removal {
    label "fatos"

    tag { run_id }

    input:
    set val(run_id), file(nonhm), read1, read2, kraken2family from input_contamremoval

    output:
    set run_id, file("${run_id}.clean.1.fq.gz"), file("${run_id}.clean.2.fq.gz"), kraken2family into channel_after_remove_contam

    script:
    """
    cat <<- EOF >> fixheaders.py
    import re
    import sys

    p = re.compile("^(@[^\\s]+)\\/([0-9]+)\$")

    prevLine = None

    for line in sys.stdin:
        m = p.match(line)
        if m and prevLine != "+\\n":
            sys.stdout.write(m.group(1) + '\\n')
        else:
            sys.stdout.write(line)
        prevLine = line
    EOF

    zcat ${read1} | python3 fixheaders.py | gzip > ${run_id}_1.fix
    zcat ${read2} | python3 fixheaders.py | gzip > ${run_id}_2.fix

    seqtk subseq ${run_id}_1.fix ${nonhm} | gzip > "${run_id}.clean.1.fq.gz"
    seqtk subseq ${run_id}_2.fix ${nonhm} | gzip > "${run_id}.clean.2.fq.gz"
    """
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

process speciation_mykrobe {
    label "clockwork"

    memory '12 GB'

    publishDir "${output_dir}/${dataset_id}/speciation/", mode: "copy", pattern: "mykrobe_output.json"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}.clean.1.fq.gz"), file("${dataset_id}.clean.2.fq.gz"), kraken2family from channel_after_remove_contam
    when:
    kraken2family == "Mycobacteriaceae"


    output:
    file("mykrobe_output.json")
    set dataset_id, file("${dataset_id}.clean.1.fq.gz"), file("${dataset_id}.clean.2.fq.gz"), file('mykrobe_output.json') into channel_pick_reference

    script:
    """
    #
    # convert pe fastq to unmapped bam
    #
    java -Xmx8G -jar /bioinf-tools/picard.jar FastqToSam FASTQ=${dataset_id}.clean.1.fq.gz FASTQ2=${dataset_id}.clean.2.fq.gz OUTPUT=${dataset_id}.bam SAMPLE_NAME=${dataset_id}

    #
    # run mykrobe
    #
    mykrobe predict tb_sample_id tb -1 ${dataset_id}.bam --format json --output mykrobe_output.json
    """
}

process pick_reference {
    //
    // given mykrobe species, find the reference directory and pass it on to other channels
    //

    label "fatos"

    publishDir "${output_dir}/${dataset_id}/speciation/", mode: "copy", pattern: "reference_info.txt"

    memory '1 GB'

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}.clean.1.fq.gz"), file("${dataset_id}.clean.2.fq.gz"), file('mykrobe_output.json') from channel_pick_reference

    output:
    set dataset_id, file("${dataset_id}.clean.1.fq.gz"), file("${dataset_id}.clean.2.fq.gz"), stdout into channel_fastqs_clean_samtools
    // set dataset_id, file("${dataset_id}.clean.1.fq.gz"), file("${dataset_id}.clean.2.fq.gz"), stdout into channel_fastqs_clean_fastqc
    set dataset_id, file("${dataset_id}.clean.1.fq.gz"), file("${dataset_id}.clean.2.fq.gz"), stdout  into channel_fastqs_clean_trim
    set dataset_id, file("${dataset_id}.clean.1.fq.gz"), file("${dataset_id}.clean.2.fq.gz"), stdout  into channel_fastqs_for_bwa
    file("reference_info.txt")

    script:
    """
    cat <<- END_OF_FILE >> pick_reference.py
    import sys, pathlib, json

    mykrobe = json.loads(sys.stdin.read())
    input_ref = sys.argv[1]
    refmap = json.loads(sys.argv[2])
    dataset_id = sys.argv[3]

    if dataset_id in refmap:
        input_ref = refmap[dataset_id]

    # -----
    species_map = [("Mycobacterium_tuberculosis", "NC_000962.3"),
                   ("Mycobacterium_africanum", "NC_000962.3"),
                   ("Mycobacterium_bovis", "NC_002945.4"),
                   ("Mycobacterium_abscessus", "NC_010397.1"),
                   ("Mycobacterium_intracellulare", "NC_016946.1"),
                   ("Mycobacterium_avium", "NC_002944.2"),
                   ("Mycobacterium_chelonae", "NZ_CP007220.1"),
                   ("Mycobacterium_kansasii", "NC_022663.1"),
                   ("Mycobacterium_fortuitum", "NZ_CP011269.1"),
                   ("Mycobacterium_chimaera", "NZ_CP012885.2")]

    clockwork_references = { 'Mycobacterium_tuberculosis',
                             'Mycobacterium_africanum',
                             'Mycobacterium_bovis' }

    species = list(mykrobe['tb_sample_id']['phylogenetics']['species'].keys())
    top_species = species[0]
    # -----

    ref = 'failed'
    for k,v in species_map:
        if k == top_species:
            ref = v

    ref_base = '/data/references/clockwork/qc_vc/'
    ref_prefix = 'Reference.'

    if input_ref == 'AUTO':
        ref_dir = pathlib.Path(ref_base) / (ref_prefix + ref)
    else:
        ref_dir = pathlib.Path(input_ref)
        ref = str(ref_dir).split('Reference.')[1]

    if ref_dir.is_dir():
        sys.stdout.write(ref)
    else:
        sys.stdout.write('failed')

    if ref in clockwork_references:
        sys.stderr.write('clockwork')
    else:
        sys.stderr.write('compass')

    def write_reference_info(top_species, dataset_id, pick_taxid, ref_dir):
        with open('reference_info.txt', 'w') as f:
            out = { "mykrobe_species": top_species,
                    "dataset_id":      dataset_id,
                    "pick_taxid":      ref,
                    "reference_dir":   str(ref_dir) + '/', }
            f.write(json.dumps(out))

    write_reference_info(top_species, dataset_id, ref, ref_dir)
    END_OF_FILE

    cat mykrobe_output.json | python3 pick_reference.py '${ref}' '${refmap}' '${dataset_id}'

    SPECIES=`cat reference_info.txt | jq -r .mykrobe_species`
    curl -s --header "Content-Type: application/json" --request POST --data "{\\"pipeline_name\\": \\"${params.pipeline_name}\\", \\"run_uuid\\":\\"${params.run_uuid}\\", \\"sample_name\\": \\"${dataset_id}\\", \\"tag_type\\": \\"mykrobe\\", \\"tag_name\\": \\"\$SPECIES\\" }" http://${params.head_node_ip}:12000/add_sample_tag > /dev/null
    """
}



process bwa {
    label "compass"
    memory '4 GB'

    echo true
    scratch true

    tag { dataset_id }

    input:
    set dataset_id, file(forward), file(reverse), ref_dir from channel_fastqs_for_bwa
    when:
    ref_dir != "failed" 

    output:
    set dataset_id, file("${dataset_id}_alignment.sam") into bwa_map

    """
    $BWA/bwa mem -R '@RG\tID:${dataset_id}\tSM:null\tLB:null\tCN:null' -t ${task.cpus} ${ref_dir}/ref.fa ${forward} ${reverse} > ${dataset_id}_alignment.sam
    """
}

process bwa_merge {
    label "compass"
    memory '13 GB'

    echo true
    scratch true

    tag { dataset_id }

    input:
    set dataset_id, file("${dataset_id}_alignment.sam") from bwa_map

    output:
    set dataset_id, file("${dataset_id}_alignment.bam"), file("${dataset_id}_seqstats.txt"), file("${dataset_id}_flagstats.txt") into bwa_merge

    """
    python $COMPASS_ROOT/nf_bwa_merge.py -b ${dataset_id}_alignment.sam -o ${dataset_id}_alignment.bam -execute -ss ${dataset_id}_seqstats.txt -fs ${dataset_id}_flagstats.txt
    """
}

process mpileup {
    memory '12 GB'

    echo true
    scratch true

    tag { dataset_id }

    input:
    set dataset_id,  file("${dataset_id}_alignment.bam") from bwa_merge
    file ref

    output:
    set dataset_id, file("${dataset_id}.out.vcf"), file("${dataset_id}.pileup.vcf") into snpcalling

    """
    python $COMPASS_ROOT/nf_mpileup.py -o 40 -e 20 -H 100 -m 2 -F 0.002 -D -S -M0 -q 30 -Q25 -E -c -g -K -L -t0.01 -i -1 -p0.5 -P full -B ${dataset_id}_alignment.bam -R ${ref_dir}/ref.fa -out ${dataset_id}.out.vcf -outpileup ${dataset_id}.pileup.vcf
    """
}

process annotvcf {
    label "compass"
    memory '2 GB'

    echo true
    scratch true

    tag { dataset_id }

    input:
    set dataset_id, file("${dataset_id}.out.vcf"), file("${dataset_id}.pileup.vcf") from snpcalling


    output:
    set dataset_id, file("${dataset_id}.annotvcf.vcf") into annotcvf

    """
    python $COMPASS_ROOT/nf_annotvcf.py -vcf ${dataset_id}.out.vcf -mpileup ${dataset_id}.pileup.vcf -refmask /data/references/compass/mask/empty_file_repmask.array -o ${dataset_id}.annotvcf.vcf -basecall -selfblastR -hqdepthinfo -lgcdepthinfo
    """
}

process basecall {
    label "compass"
    memory '1 GB'

    echo true
    scratch true

    publishDir "${params.output_dir}/${dataset_id}", mode: "move"

    tag { dataset_id }

    input:
    set dataset_id, file("${dataset_id}.annotvcf.vcf") from annotcvf

    output:
    file("${dataset_id}.consensus.fasta.gz")
    file("${dataset_id}.basecall.vcf.gz")

    """
    python $COMPASS_ROOT/nf_basecall.py -A 25 -e DISABLED -E DISABLED -g DISABLED -G DISABLED -K0.90 -J DISABLED -invcf ${dataset_id}.annotvcf.vcf -outvcf ${dataset_id}.basecall.vcf.gz -outvcfIndel ${dataset_id}.basecall_Indel.vcf.gz -outfasta ${dataset_id}.consensus.fasta.gz -outstats ${dataset_id}.basecallstats.txt -u ${dataset_id} -refid ${ref_name} -Q30 -q30 -m30 -n5 -S25 -I25 -z -B1 -p -N -f0.35
    """
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////