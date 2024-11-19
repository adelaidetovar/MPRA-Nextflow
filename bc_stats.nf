#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

nextflow.enable.dsl = 2

/*
====================================================
        MPRA-NextFlow Barcode Stats Pipeline
====================================================
#### Documentation
https://github.com/adelaidetovar/MPRA-NextFlow
#### Last updated
2024-09-05
#### Authors
Jacob Kitzman <kitzmanj@umich.edu>
Adelaide Tovar <tovar@umich.edu>
*/

// def help message

def helpMsg () {
    log.info'''
    ====================================================
            MPRA-NextFlow Barcode Stats Pipeline
    ====================================================

    Usage:

    nextflow run bc_stats.nf -params-file input.json

    Recommended arguments:
        --outdir                    Output directory for results (default: bc_count)
        --single_end                Use this flag if you have single-end reads (default: false)

    Other options:
        --ignore_umi                Use this flag if you do not want to pull UMIs and perform deduplication (default: false)
        --checksampbysamp        Use this flag if you want to calculate pairwise stats across samples and generate plots (default: false)

    Extra arguments:
        --h, --help                 Print help message
    '''.stripIndent()
}

if (params.containsKey('h') || params.containsKey('help')) {
    helpMsg()
    exit 0
}

// check params

// set output directory
if (params.containsKey('outdir')) {
    params.outdir = params['outdir']
} else {
    params.outdir = 'bc_count'
}

// ignore umi for deduplication
if (params.containsKey('ignore_umi')) {
    params.ignore_umi = true
} else {
    params.ignore_umi = false
}

if (params.containsKey('checksampbysamp')) {
    params.checksampbysamp = true
} else {
    params.checksampbysamp = false
}

/*
step 1: calculate read statistics across filtering and clustering steps
*/

process read_stats_bc {
    publishDir "${params.outdir}/postprocess", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 1
    memory '20 GB'
    time '30m'

    input:
    path outdir_ch

    output:
    path('read_stats.txt')

    script:
    """
    read_stats \
        --in_clip ${outdir_ch}/clip \
        --in-bc ${outdir_ch}/bc \
        --out_tab "read_stats.txt"
    """
}

process read_stats_umibc {
    publishDir "${params.outdir}/postprocess", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 1
    memory '20GB'
    time '30m'

    input:
    path outdir_ch

    output:
    path('read_stats.txt')

    script:
    """
    read_stats \
        --in_clip ${outdir_ch}/clip \
        --in_umi ${outdir_ch}/umi \
        --in_umibc ${outdir_ch}/parsed \
        --in_clust ${outdir_ch}/bc \
        --out_tab "read_stats.txt"
    """
}

/*
step 2: barcode stats
*/

process bc_stats {
    publishDir "${params.outdir}/postprocess", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 1
    memory '20GB'
    time '30m'

    input:
    path outdir_ch

    output:
    path('bc_stats.txt')

    script:
    """
    bc_stats \
        --in_dir ${outdir_ch}/bc \
        --out_tab "bc_stats.txt"
    """
}

/*
step 3: filter clustered barcodes
*/

process filter_bcgroups {
    publishDir "${params.outdir}/postprocess/filt_bc", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 4
    memory '24GB'
    tag "$libname"

    input:
    val(libname)

    output:
    tuple val(libname), path("${libname}.bc_cluster.filt.txt")

    script:
    """
    filter_starcode_bcgroups \
        --in_histo ${params.outdir}/bc/${libname}.bc_cluster.txt \
        --out_histo ${libname}.bc_cluster.filt.txt \
        --out_summary ${libname}.summ.txt
    """
}

/*
step 4: make waterfall plots with clustered barcodes
*/

process waterfall_plot {
    publishDir "${params.outdir}/postprocess/plots", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 2
    memory '24GB'
    tag "$libname"

    input:
    val(libname)

    output:
    tuple val(libname), path("${libname}.waterfall.png")

    script:
    """
    plot_tag_count_histos \
        --linCountHisto ${params.outdir}/bc/${libname}.bc_cluster.txt \
        --lNames ${libname} \
        --annotMaxY \
        --outPlot ${libname}.waterfall.png
    """
}

/*
step 5: combine clustered barcodes and counts across all samples into a joint matrix
*/

process joint_mtx {
    publishDir "${params.outdir}/postprocess/joint_mtx", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 4
    memory { 20.GB + (20.GB * task.attempt) }
    maxRetries 3
    time '2h'
    errorStrategy 'retry'

    input:
    path outdir_ch

    output:
    path("sampletab.txt"), emit: sample_tab
    path("joint_matrix.txt"), emit: joint_matrix

    script:
    """
    echo -e "sampleName\thistoPath" > sampletab.txt
    for filepath in \$(ls ${outdir_ch}/bc/*.bc_cluster.txt); do
        libname=\$(basename "\$filepath" .bc_cluster.txt)
        echo -e "\${libname}\t\$filepath" >> sampletab.txt
    done

    ${IONICE} join_starcode_histos \
        --inSampleKey sampletab.txt \
        --colSampname "sampleName" \
        --colHistoPath "histoPath" \
        --outHisto joint_matrix.txt
    """
}

/*
step 6: matrix math to generate barcode statistics across samples
*/

process sampbysamp_mtx {
    publishDir "${params.outdir}/postprocess/joint_mtx", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 4
    memory { 20.GB + (20.GB * task.attempt) }
    maxRetries 3
    time '6h'
    errorStrategy 'retry'

    input:
    path(sample_tab)
    path(joint_matrix)

    output:
    path('sampbysamp.*')

    script:
    """
    ${IONICE} jointhisto_to_sampbysamp_mtx \
        --inKey $sample_tab \
        --colLibname "sampleName" \
        --inJointHisto $joint_matrix \
        --outBase sampbysamp
    """
}

/*
step 7: plot stats generated in step 10/12
*/

process sampbysamp_plot {
    publishDir "${params.outdir}/postprocess/plots", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 4
    memory { 20.GB + (20.GB * task.attempt) }
    maxRetries 3
    time '6h'
    errorStrategy 'retry'

    input:
    path("${params.outdir}/postprocess/joint_mtx/sampbysamp.*")

    output:
    path('sampbysamp.*.png')

    script:
    """
    ${IONICE} samp2samp_mtx_cluster_plots \
        --inSamp2sampMtxBase ${params.outdir}/postprocess/joint_mtx/sampbysamp \
        --outPlotBase sampbysamp
    """
}

workflow {
    libnames = params.libnames.keySet()
    libname = Channel.from(libnames)

    outdir_ch = Channel.fromPath("${params.outdir}")

    if (params.ignore_umi) {
        read_stats_bc(outdir_ch)
} else {
        read_stats_umibc(outdir_ch)
    }

    bc_stats(outdir_ch)
    filter_bcgroups(libname)
    waterfall_plot(libname)

    if (params.checksampbysamp) {
        joint_mtx(outdir_ch) | sampbysamp_mtx | sampbysamp_plot
    }
}

// generate run summary report

workflow.onComplete {
    summary = """
    ===================================================
            MPRA-NextFlow Barcode Stats Summary
    ===================================================
    Time completed: ${workflow.complete}
    Duration:       ${workflow.duration}
    Success:        ${workflow.success}
    Exit status:    ${workflow.exitStatus}
    """

    println summary
}

