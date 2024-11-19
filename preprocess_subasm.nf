#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

nextflow.enable.dsl = 2

/*
========================================================
        MPRA-NextFlow Barcode Counting Pipeline
========================================================
#### Documentation
https://github.com/adelaidetovar/MPRA-NextFlow
#### Last updated
2024-09-11
#### Authors
Jacob Kitzman <kitzmanj@umich.edu>
Adelaide Tovar <tovar@umich.edu>
*/

// def help message

def helpMsg () {
    log.info'''
    ========================================================
            MPRA-NextFlow Barcode Counting Pipeline
    ========================================================

    Usage:

    nextflow run preprocess_subasm.nf -params-file sa_input.json

    Mandatory arguments:
      --fastadir                                Fasta directory - in quotes
      --bc-len
      --seq-adapt-a
      --seq-adapt-b
      --seq-adapt-c

    Recommended arguments:
        --outdir                    Output directory for results (default: preprocess_subasm)

    Extra arguments:
        --h, --help                 Print help message
    '''.stripIndent()
}

if (params.containsKey('h') || params.containsKey('help')){
    helpMsg()
    exit 0
}

process first_clip {

    publishDir "${params.outdir}/clip", overwrite: true
    cpus 8
    memory '32 GB'
    time '1h'
    queue 'standard'
    tag "$libname"

    input:
    tuple val(libname), path(fwd_fq), path(rev_fq)

    output:
    tuple val(libname), path("${libname}.clip.r*.gz"), path("${libname}.clip.log"), emit: clip_out

    script:
    """
    cutadapt \
        --adapter ggcgggccgCCGA...NNNNNNNNNNNNNNNCCGCGATCGCCTTACCGGTC \
        --adapter GACCGGTAAGGCGATCGCGGNNNNNNNNNNNNNNN...TCGGCGGCCCGCC \
        --discard-untrimmed \
        --action=lowercase \
        -e 0.1 -O 44 \
        -j 10 \
        --minimum-length 10 \
        -q 16 \
        -o ${libname}.clip.r2.gz \
        -p ${libname}.clip.r1.gz \
        ${rev_fq} ${fwd_fq}\
            > ${libname}.clip.log
    """
}

process extract_bc {
    
    publishDir "${params.results}/bc", mode: 'rellink'
    cpus 8
    memory '10 GB'
    time '12h'
    queue 'standard'
    tag "$libname"

    input:
    tuple val(libname), path(clip_out)

    output:
    tuple val(libname), path("${libname}.bc.r*.gz"), emit: extract_out

    script:
    """
    ${IONICE} extract_bc_and_tag_fqs \
        --in_bc_fq ${params.outdir}/clip/${libname}.clip.r2.gz \
        --in_paired_fq ${params.outdir}/clip/${libname}.clip.r1.gz \
        --out_bc_fq ${libname}.bc.r2.gz \
        --out_paired_fq ${libname}.bc.r1.gz \
        --min_bc_len ${params.min_bc_len} \
        --max_bc_len ${params.max_bc_len}
    """
}

process second_clip {

    publishDir "${params.outdir}/clip", overwrite: true
    cpus 8
    memory '48 GB'
    time '1h'
    queue 'standard'
    tag "$libname"

    input:
    tuple val(libname), path(extract_out)

    output:
    tuple val(libname), path("${libname}.clip2.r*.gz"), path("${libname}.clip2.log"), emit: clip2_out

    script:
    """
    cutadapt \
        -G 'tctagaggttcgagaccgcgatcgccttaccggtctcgcagtg;min_overlap=40' \
        --times 2 \
        --discard-untrimmed \
        --action=trim \
        -e 0.1 \
        -j 10 \
        --minimum-length 10 \
        -q 10 \
        -o ${libname}.clip2.r1.gz \
        -p ${libname}.clip2.r2.gz \
        ${libname}.bc.r1.gz ${libname}.bc.r2.gz \
            > ${libname}.clip2.log
    """
}

process bcs_to_fq {
    
    publishDir "${params.results}/clip2", mode: 'rellink'
    container 'library://'
    errorStrategy 'retry'
    maxRetries 1
    time '2h'
    maxForks 10
    memory '40GB'
    cpus '4'

    input:
    tuple val(libname), path(clip2_out)

    output:
    tuple val(libname), path("${libname}.plasbc.temp.gz"), emit: plasbc_temp

    shell:
    """
    bcs_from_fqhdr_to_fq \
        --in_fq ${libname}.clip2.r1.gz \
        --out_bc_fq ${libname}.plasbc.temp.gz
    """
}

process fix_header {
    
    publishDir "${params.results}/clip2", mode: 'rellink'
    container 'library://'
    errorStrategy 'retry'
    maxRetries 1
    time ''
    maxForks 10
    memory ''

    input:
    tuple val(libname), path(plasbc_temp)

    output:
    tuple val(libname), path(${libname}.plasbc.gz)

    """
    bcs_from_fqhdr_to_fq \
      --in_fq ${subasm}.oligo.fq.gz \
      --out_bc_fq ${subasm}.bc.fq.gz
    """
}

process pre_starcode {
    
    publishDir "${params.results}/pre_starcode", mode: 'rellink'
    container 'library://'
    errorStrategy 'retry'
    maxRetries 1
    time ''
    maxForks 10
    memory ''

    input:
        path "${subasm}.bc.fq.gz" from bcfq

    output:
        path "${subasm}.bc_precluster.txt" into prebcclust

    """
    starcode -t 32 --print-clusters -d 3 -i <( pigz -d -c -p 40 ${subasm}.bc.fq.gz | \
        sed -ne '2~4p' ) -o ${subasm}.bc_precluster.txt
    """
}

process pre_histo {
    
    publishDir "${params.results}/pre_starcode", mode: 'rellink'
    container 'library://'
    errorStrategy 'retry'
    maxRetries 1
    time ''
    maxForks 10
    memory ''

    input:
        path "${subasm}.bc_precluster.txt" from prebcclust

    output:
        path "${subasm}.bc_hist.png"

    shell:
        """
        python /home/tovar/toolscripts/plot_tag_count_histos.py \
            --linCountHisto ${subasm}.bc_precluster.txt \
            --lNames ${subasm} \
            --annotMaxY  \
            --outPlot ${subasm}.bc_hist.png
        """
}

workflow {
    libraries = params.libraries.keySet()
    subasm_fq1_fq2 = []

    for (subasm in libraries) {
        for (read in subasm_to_readgroups(subasm)) {
        fastqs = subasm_and_readgroup_to_fastqs(subasm, read)
        first_read = fastqs['1']
        second_read = fastqs['2']
        subasm_fq1_fq2 << [subasm, read, path(first_insert), path(second_insert)]
        }
    }

    subasm_fq1_fq2 = Channel.from(subasm_fq1_fq2)

    merged = merge(subasm_fq1_fq2)
    filtered = filter_reads(merged)
    pre_bc = tag_readnames_w_bc(filtered) | remove_flanks | extract_bc
    pre_clustered = pre_starcode(pre_bc)
    pre_plot = pre_histo(pre_clustered)

    post_clustered = post_starcode(pre_bc)
    clusters = load_clusters(post_clustered)
    aligned = align()
    
}