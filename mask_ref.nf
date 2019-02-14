#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    =========================================
              MMM Compass Compact 
            Generate Masked Reference
    =========================================
    Usage:

    nextflow run mask_ref.nf --help
    nextflow run mask_ref.nf --test -profile test_docker
    nextflow run mask_ref.nf \
    --output_dir tests/data/output_dir \
    --ref tests/data/reference/NC_000962_2.fasta \
    --mask true \
    -profile test_docker

    Mandatory arguments:
    --ref                    FILE         reference genome

    Optional arguments:
    --mask                   Boolean      use self-blast to mask repeated region, default true
    """.stripIndent()
}

params.help = false
params.test = false
params.mask = true


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

if (params.test){
    params.ref = "tests/data/reference/NC_000962_3.fasta"
    params.output_dir = "tests/data/output_dir"
}else{
    params.ref = ""
}

ref = file(params.ref)
ref_folder = ref.getParent()
ref_name = ref.getBaseName()

if (params.mask == true)
{
    process generatemaskarray {
        echo true
        scratch true

        publishDir "${params.output_dir}/ref_mask/${ref.getBaseName()}", mode: "copy" , pattern: "${ref.getBaseName()}*"

        tag {ref}

        input:
        file ref
        
        output:
        file("${ref_name}_repmask.array") into mask_ref
        
        """
        python ${COMPASS_ROOT}/nf_ref_index.py -r ${ref}
        """
    }
}
else{
    process generatenomaskarray {
        echo true
        scratch true

        publishDir "${params.output_dir}/ref_mask/", mode: "copy" 

        tag {ref}

        input:
        file ref
        
        output:
        file("${ref_name}_repmask.array") into mask_ref
        
        """
        python3 ${COMPASS_ROOT}/nf_ref_nomask.py -r ${ref}
        """
    }
}