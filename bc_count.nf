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
2024-09-05
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

    nextflow run bc_count.nf -params-file input.json

    Recommended arguments:
        --outdir                    Output directory for results (default: bc_count)
        --single_end                Use this flag if you have single-end reads (default: false)

    Other options:
        --umi-len                   UMI length (default: 15)
        --bc-len                    MPRA barcode length (default: 20)
        --umi-dist                  Levenshtein distance used for clustering UMIs (default: 2)
        --bc-dist                   Levenshtein distance used for clustering (deduplicated) barcodes (default: 2)
        --ignore-umi                Use this flag if you do not want to pull UMIs and perform deduplication (default: false)
        --checksampbysamp           Use this flag if you want to calculate pairwise stats across samples and generate plots (default: false)

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

// set umi length for umi_tools/starcode-umi
if (params.containsKey('umi-len')) {
    params.umi_len = params['umi-len']
} else {
    params.umi_len = 15
}

// set bc length for starcode-umi/starcode
if (params.containsKey('bc-len')) {
    params.bc_len = params['bc-len']
} else {
    params.bc_len = 20
}

// set umi levenshtein distance for starcode-umi
if (params.containsKey('umi-dist')) {
    params.umi_dist = params['umi-dist']
} else {
    params.umi_dist = 2
}

// set bc levenshtein distance for starcode
if (params.containsKey('bc-dist')) {
    params.bc_dist = params['bc-dist']
} else {
    params.bc_dist = 2
}

// single end or paired end reads
if (params.containsKey('single-end')) {
    params.single_end = true
} else {
    params.single_end = false
}

// ignore umi for deduplication
if (params.containsKey('ignore-umi')) {
    params.ignore_umi = true
} else {
    params.ignore_umi = false
}

if (params.containsKey('checksampbysamp')) {
    params.checksampbysamp = true
} else {
    params.checksampbysamp = false
}

// define helper functions

def libname_to_fwd (libname) {
    return(params.libnames[libname].fwd.values().toList().collect())
}

def libname_to_rev (libname) {
    return(params.libnames[libname].rev.values().toList().collect())
}

def libname_to_umi_len (libname) {
    return(params.libnames[libname].umi_len)
}

/*
step 1: if necessary, concatenate fastq files from multiple lanes
*/

process cat_fq {
    publishDir "${params.outdir}/cat_fq", overwrite: true
    cpus 2
    memory '4 GB'
    time '15m'
    queue 'standard'
    tag "$libname"

    input:
    tuple val(libname), path(fwd_files), path(rev_files)

    output:
    tuple val(libname), path("${libname}.cat.r*.gz"), emit: concat_fq

    script:
    """
    cat ${fwd_files.join(' ')} > ${libname}.cat.r1.gz
    cat ${rev_files.join(' ')} > ${libname}.cat.r2.gz
    """
}

process cat_fq_umi {
    publishDir "${params.outdir}/cat_fq", overwrite: true
    cpus 2
    memory '4 GB'
    time '15m'
    queue 'standard'
    tag "$libname"

    input:
    tuple val(libname), val(umi_len), path(fwd_files), path(rev_files)

    output:
    tuple val(libname), val(umi_len), path("${libname}.cat.r*.gz"), emit: concat_fq

    script:
    """
    cat ${fwd_files.join(' ')} > ${libname}.cat.r1.gz
    cat ${rev_files.join(' ')} > ${libname}.cat.r2.gz
    """
}

/*
step 2: filter for sequences containing expected adapters
*/

process bc_clip {
    publishDir "${params.outdir}/clip", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 16
    memory '80 GB'
    time '60m'
    queue 'standard'
    tag "$libname"

    input:
    tuple val(libname), val(umi_len), path(concat_fq)

    output:
    tuple val(libname), path("${libname}.clip.r*.gz"), emit: clip_out
    path("${libname}.clip.log"), emit: clip_log
    
    script:
    """
    ${IONICE} cutadapt \
        -g CCGGTACTGTTGGTAAAGAACCTCTAGA...TCGGCNGCCC \
        -g CCGGTACTGTTGGTAAAGAACGGAAGA...TCGGCNGCCC \
        -g CCGGTACTGTTGGTAAAGAACCACCAGA...TCGGCNGCCC \
        -a TCGGCGGCCCGCCACCATGG \
        -e 0.1 -O 25 \
        -j 8 \
        --minimum-length 18 \
        --maximum-length 22 \
        --discard-untrimmed \
        -q 10 \
        --pair-filter=first \
        -o ${libname}.clip.r1.gz \
        -p ${libname}.clip.r2.gz \
        ${params.outdir}/cat_fq/${libname}.cat.r1.gz \
        ${params.outdir}/cat_fq/${libname}.cat.r2.gz \
        > ${libname}.clip.log
    """
}

process bc_clip_umi {
    publishDir "${params.outdir}/clip", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 16
    memory '80 GB'
    time '60m'
    queue 'standard'
    tag "$libname"

    input:
    tuple val(libname), val(umi_len), path(concat_fq)

    output:
    tuple val(libname), val(umi_len), path("${libname}.clip.r*.gz"), emit: clip_out
    path("${libname}.clip.log"), emit: clip_log

    script:
    """
    ${IONICE} cutadapt \
        -a TCGGCGGCCCGCCACCATGG \
        -g CCGGTACTGTTGGTAAAGAACCTCTAGA...TCGGCNGCCC \
        -g CCGGTACTGTTGGTAAAGAACGGAAGA...TCGGCNGCCC \
        -g CCGGTACTGTTGGTAAAGAACCACCAGA...TCGGCNGCCC \
        -e 0.1 -O 25 \
        -j 8 \
        --minimum-length 18 \
        --maximum-length 22 \
        --discard-untrimmed \
        -q 10 \
        --pair-filter=first \
        -o ${libname}.clip.r1.gz \
        -p ${libname}.clip.r2.gz \
        ${params.outdir}/cat_fq/${libname}.cat.r1.gz \
        ${params.outdir}/cat_fq/${libname}.cat.r2.gz \
        > ${libname}.clip.log
        """
}

/*
step 3: extract barcode seq and move to file used as input for starcode
*/

process prep_bc {
    publishDir "${params.outdir}/parsed", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    tag "$libname"
    cpus 10
    memory '50GB'
    time '1h'

    input:
    tuple val(libname), path(clip_out)

    output:
    tuple val(libname), path("${libname}.bc.txt"), emit: sepbc

    script:
    """
    seqtk seq -a ${libname}.clip.r1.gz | \
        tee >(awk 'NR % 2 == 0' > ${libname}.bc.txt)
    """
}

process extract_umi {
    publishDir "${params.outdir}/umi", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 16
    memory { 60.GB + (30.GB * task.attempt) }
    time { 8.hours + (4.hours * task.attempt) }
    tag "$libname"
    errorStrategy 'retry'
    maxRetries 3
    queue 'standard'

    input:
    tuple val(libname), val(umi_len), path(clip_out)

    output:
    tuple val(libname), val(umi_len), path("${libname}.umi.r1.gz"), emit: fq_fwd_umi
    path("${libname}.umi.log"), emit: umi_log

    script:
    """
    ${IONICE} umi_tools extract \
        --stdin ${libname}.clip.r2.gz \
        --stdout ${libname}.umi.r2.gz \
        --extract-method regex \
        --bc-pattern='(?P<umi_1>.{${umi_len}})(?P<discard_1>GCTCCTCGCCCTTGCTCA|TCGCCCTTGCTCACCATG){s<=2}' \
        --read2-in ${libname}.clip.r1.gz \
        --read2-out=${libname}.umi.r1.gz \
        -L /dev/stdout | tee ${libname}.umi.log
    """
}

process prep_umibc {
    publishDir "${params.outdir}/parsed", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    tag "$libname"

    input:
    tuple val(libname), val(umi_len), path(fq_fwd_umi)

    output:
    tuple val(libname), val(umi_len), path("${libname}.umibc.txt"), emit: umibc

    script:
    """
    seqtk seq -a ${fq_fwd_umi} | tee \
        >(awk 'NR % 2 == 1' | grep -Po '(?<=_)[ATGCN]*' > ${libname}.umi.txt) \
        >(awk 'NR % 2 == 0' > ${libname}.bc.txt)

    if [[ -s ${ libname }.umi.txt && -s ${ { libname } }.bc.txt ]]; then
        paste -d "" ${libname}.umi.txt ${libname}.bc.txt > ${libname}.umibc.txt
    else
        echo "One or both files empty." >&2
        exit 1
    fi
    """
}

/*
step 4: cluster barcodes using starcode
*/

process starcode_bc {
    publishDir "${params.outdir}/bc", overwrite: true, mode: "copy"
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 10
    memory { 30.GB + 30.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    tag "$libname"

    input:
    tuple val(libname), path(sepbc)

    output:
    tuple val(libname), path("${libname}.bc_cluster.txt"), emit: bc_cluster

    script:
    """
    ${IONICE} cat ${sepbc} | awk '{print substr(\$1, 1,${params.bc_len})}' | \
    starcode \
        -t 10 \
        --print-clusters \
        -d ${params.bc_dist} \
        -i /dev/stdin \
        -o ${libname}.bc_cluster.txt
    """
}

process starcode_umi {
    publishDir "${params.outdir}/parsed", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240716'
    tag "$libname"
    cpus 10
    memory { 80.GB + (40.GB * task.attempt) }
    time '12h'
    errorStrategy 'retry'
    maxRetries 3
    queue { task.attempt > 2 ? 'largemem' : 'standard' }

    input:
    tuple val(libname), val(umi_len), path(umibc)

    output:
    tuple val(libname), val(umi_len), path("${libname}.starumi"), emit: umi_out

    script:
    """
    ${IONICE} starcode-umi \
        --starcode-path /usr/local/bin/starcode \
        --umi-len ${umi_len} \
        --umi-d ${params.umi_dist} \
        ${umibc} > ${libname}.starumi

    if [ ! -s ${libname}.starumi ]; then
        echo "starcode-umi failed. Exiting!" >&2
        exit 1
    fi
    """
}

process starcode_umibc {
    publishDir "${params.outdir}/bc", overwrite: true, mode: "copy"
    container 'library://tovar/general/mpra_utilscripts:20240716'
    tag "$libname"
    cpus 10
    memory { 80.GB + (25.GB * task.attempt) }
    time { 8.hours + (4.hours * task.attempt) }
    errorStrategy 'retry'
    maxRetries 3
    queue 'standard'

    input:
    tuple val(libname), val(umi_len), path(umi_out)

    output:
    tuple val(libname), path("${libname}.bc_cluster.txt"), emit: bc_cluster

    script:
    """
    ${IONICE} cat ${libname}.starumi | awk '{print substr(\$1, ${umi_len}+1, ${params.bc_len})}' | \
    starcode \
        -t 10 \
        --print-clusters \
        -d ${params.bc_dist} \
        -i /dev/stdin \
        -o ${libname}.bc_cluster.txt
    """
}

/*
step 5: calculate read statistics across filtering and clustering steps
*/

process read_stats_bc {
    publishDir "${params.outdir}/postprocess", overwrite: true, mode: "copy"
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 1
    memory '20 GB'
    time '30m'

    input:
    path(clip_log)
    tuple val(libname), path(bc_cluster)
    val(params.outdir)

    output:
    path('read_stats.txt')

    script:
    """
    read_stats \
        --in_clip ${params.outdir}/clip \
        --in-bc ${params.outdir}/bc \
        --out_tab "read_stats.txt"
    """
}

process read_stats_umibc {
    publishDir "${params.outdir}/postprocess", overwrite: true, mode: "copy"
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 1
    memory '20GB'
    time '30m'

    input:
    path(clip_log)
    path(umi_log)
    tuple val(libname), val(umi_len), path(umibc)
    tuple val(libname), path(bc_cluster)
    val(params.outdir)

    output:
    path('read_stats.txt')

    script:
    """
    read_stats \
        --in_clip ${params.outdir}/clip \
        --in_umi ${params.outdir}/umi \
        --in_umibc ${params.outdir}/parsed \
        --in_clust ${params.outdir}/bc \
        --out_tab "read_stats.txt"
    """
}

/*
step 6: barcode stats
*/

process bc_stats {
    publishDir "${params.outdir}/postprocess", overwrite: true, mode: "copy"
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 1
    memory '20GB'
    time '30m'

    input:
    tuple val(libname), path(bc_cluster)
    val(params.outdir)

    output:
    path('bc_stats.txt')

    script:
    """
    bc_stats \
        --in_dir ${params.outdir}/bc \
        --out_tab "bc_stats.txt"
    """
}

/*
step 7: filter clustered barcodes
*/

process filter_bcgroups {
    publishDir "${params.outdir}/postprocess/filt_bc", overwrite: true, mode: "copy"
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 4
    memory '24GB'
    tag "$libname"

    input:
    tuple val(libname), path(bc_cluster)
    val(params.outdir)

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
step 8: make waterfall plots with clustered barcodes
*/

process waterfall_plot {
    publishDir "${params.outdir}/postprocess/plots", overwrite: true, mode: "copy"
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 2
    memory '24GB'
    tag "$libname"

    input:
    tuple val(libname), path(bc_cluster)
    val(params.outdir)

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
step 9: combine clustered barcodes and counts across all samples into a joint matrix
*/

process joint_mtx {
    publishDir "${params.outdir}/postprocess/joint_mtx", overwrite: true, mode: "copy"
    container 'library://tovar/general/mpra_utilscripts:20240716'
    cpus 4
    memory { 20.GB + (20.GB * task.attempt) }
    maxRetries 3
    time '2h'
    errorStrategy 'retry'

    input:
    tuple val(libname), path(bc_cluster)
    val(params.outdir)

    output:
    path("sampletab.txt"), emit: sample_tab
    path("joint_matrix.txt"), emit: joint_matrix

    script:
    """
    echo -e "sampleName\thistoPath" > sampletab.txt
    for filepath in \$(ls ${params.outdir}/bc/*.bc_cluster.txt); do
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
step 10: matrix math to generate barcode statistics across samples
*/

process sampbysamp_mtx {
    publishDir "${params.outdir}/postprocess/joint_mtx", overwrite: true, mode: "copy"
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
step 11: plot stats generated in step 10
*/

process sampbysamp_plot {
    publishDir "${params.outdir}/postprocess/plots", overwrite: true, mode: "copy"
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

    // create sample channels
    if (params.ignore_umi) {
        fq_in = Channel.from(
        libnames.collectMany { libname ->
            fwd_files = libname_to_fwd(libname)
            rev_files = libname_to_rev(libname)

            [libname, fwd_files, rev_files]
        }
    )

        clip_in = cat_fq(fq_in)
    } else {
        fq_in = Channel.from(
        libnames.collectMany { libname ->
            fwd_files = libname_to_fwd(libname)
            rev_files = libname_to_rev(libname)
            umi_len = libname_to_umi_len(libname)

            [[libname, umi_len, fwd_files, rev_files]]
        }
    )

        clip_in = cat_fq_umi(fq_in)
    }

    if (params.ignore_umi) {
        (clip_out, clip_log) = bc_clip(clip_in)
        starcode_out = prep_bc(clip_out) | starcode_bc
        clip_ch = clip_log.collect()
        bc_out = starcode_out.collect()
        read_stats_bc(clip_log, bc_out, params.outdir)
    } else {
        (clip_out, clip_log) = bc_clip_umi(clip_in)
        (fq_fwd_umi, umi_log) = extract_umi(clip_out)
        starcode_out = prep_umibc(fq_fwd_umi) | starcode_umi | starcode_umibc
        clip_ch = clip_log.collect()
        umi_ch = umi_log.collect()
        parsed_out = prep_umibc.out.umibc.collect()
        bc_out = starcode_out.collect()
        read_stats_umibc(clip_ch, umi_ch, parsed_out, bc_out, params.outdir)
    }

    bc_stats(bc_out, params.outdir)
    filter_bcgroups(starcode_out, params.outdir)
    waterfall_plot(starcode_out, params.outdir)
    if (params.checksampbysamp) {
        joint_mtx(bc_out, params.outdir) | sampbysamp_mtx | sampbysamp_plot
    }

}

// generate run summary report

workflow.onComplete {
    summary = """
    ========================================================
            MPRA-NextFlow Barcode Counting Summary
    ========================================================
    Time completed: ${workflow.complete}
    Duration:       ${workflow.duration}
    Success:        ${workflow.success}
    Exit status:    ${workflow.exitStatus}
    """

    println summary
}
