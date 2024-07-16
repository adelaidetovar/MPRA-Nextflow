#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

nextflow.enable.dsl = 2

/*
========================================================
        MPRA-NextFlow Barcode Counting Pipeline
========================================================
#### Documentation
https://github.com/adelaidetovar/MPRA-NextFlow
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
        --short_reads               Use this flag if you used custom primers to generate short reads instead of PE 150bp reads (default: false)

    Other options:
        --umi_len                   UMI length (default: 15)
        --bc_len                    MPRA barcode length (default: 20)
        --umi_dist                  Levenshtein distance used for clustering UMIs (default: 2)
        --bc_dist                   Levenshtein distance used for clustering (deduplicated) barcodes (default: 2)
        --ignore_umi                Use this flag if you do not want to pull UMIs and perform deduplication
        --check_samp_by_samp        Use this flag if you want to calculate pairwise stats across samples and generate plots

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
if (params.containsKey('umi_len')) {
    params.umi_len = params['umi_len']
} else {
    params.umi_len = 15
}

// set bc length for starcode-umi/starcode
if (params.containsKey('bc_len')) {
    params.bc_len = params['bc_len']
} else {
    params.bc_len = 20
}

// set umi levenshtein distance for starcode-umi
if (params.containsKey('umi_dist')) {
    params.umi_dist = params['umi_dist']
} else {
    params.umi_dist = 2
}

// set bc levenshtein distance for starcode
if (params.containsKey('bc_dist')) {
    params.bc_dist = params['bc_dist']
} else {
    params.bc_dist = 2
}

// single end or paired end reads
if (params.containsKey('single_end')) {
    params.single_end = true
} else {
    params.single_end = false
}

// short custom reads or long PE 150 reads
if (params.containsKey('short_reads')) {
    params.short_reads = true
} else {
    params.short_reads = false
}

// ignore umi for deduplication
if (params.containsKey('ignore_umi')) {
    params.ignore_umi = true
} else {
    params.ignore_umi = false
}

if (params.containsKey('check_samp_by_samp')){
    params.check_samp_by_samp = true
} else {
    params.check_samp_by_samp = false
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
check whether to use umi
*/

if (params.ignore_umi) {

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
        tuple val(libname), path("${libname}.cat.r*.gz")

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
        container 'library://tovar/general/mpra_utilscripts:20240619'
        cpus 16
        memory '80 GB'
        time '60m'
        queue 'standard'
        tag "$libname"

        input:
        tuple val(libname), val(umi_len), path("${params.outdir}/cat_fq/${libname}.cat.r*.gz")

        output:
        tuple val(libname), path("${libname}.clip.r*.gz"), path("${libname}.clip.log"), emit: clip_out

        script:
        if (params.short_reads) {
            """
            ${IONICE} cutadapt \
                -a TCGGCGGCCCGCCACCATGG \
                -e 0.1 -O 25 \
                -j 16 \
                --minimum-length 16 \
                --maximum-length 24 \
                --discard-untrimmed \
                -q 10 \
                --pair-filter=first \
                -o ${libname}.clip.r1.gz \
                -p ${libname}.clip.r2.gz \
                ${params.outdir}/cat_fq/${libname}.cat.r1.gz \
                ${params.outdir}/cat_fq/${libname}.cat.r2.gz \
                > ${libname}.clip.log
            """
        } else {
            """
            ${IONICE} cutadapt \
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
    }

    /*
    step 3: extract barcode seq and move to file used as input for starcode
    */

    process prep_bc {
        publishDir "${params.outdir}/parsed", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
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

    /*
    step 4: cluster barcodes using starcode
    */

    process starcode_bc {
        publishDir "${params.outdir}/bc", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
        cpus 10
        memory { 30.GB + 30.GB * task.attempt }
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

    /*
    step 5: calculate read statistics across filtering and clustering steps
    */

    process read_stats_bc {
        publishDir "${params.outdir}/postprocess", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
        cpus 1
        memory '20 GB'
        time '30m'

        input:
        tuple val(libname), path(bc_cluster)

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

} else {

    /*
    step 1: if necessary, concatenate fastq files from multiple lanes
    */

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
        tuple val(libname), val(umi_len), path("${libname}.cat.r*.gz")

        script:
        """
        cat ${fwd_files.join(' ')} > ${libname}.cat.r1.gz
        cat ${rev_files.join(' ')} > ${libname}.cat.r2.gz
        """
    }

    /*
    step 2: filter for sequences containing expected adapters
    */

    process bc_clip_umi {
        publishDir "${params.outdir}/clip", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
        cpus 16
        memory '80 GB'
        time '60m'
        queue 'standard'
        tag "$libname"

        input:
        tuple val(libname), val(umi_len), path("${params.outdir}/cat_fq/${libname}.cat.r*.gz")

        output:
        tuple val(libname), val(umi_len), path("${libname}.clip.r*.gz"), path("${libname}.clip.log")

        script:
        if (params.short_reads) {
            """
            ${IONICE} cutadapt \
                -a TCGGCGGCCCGCCACCATGG \
                -e 0.1 -O 25 \
                -j 16 \
                --minimum-length 16 \
                --maximum-length 24 \
                --discard-untrimmed \
                -q 10 \
                --pair-filter=first \
                -o ${libname}.clip.r1.gz \
                -p ${libname}.clip.r2.gz \
                ${params.outdir}/cat_fq/${libname}.cat.r1.gz \
                ${params.outdir}/cat_fq/${libname}.cat.r2.gz \
                > ${libname}.clip.log
            """
            } else {
            """
            ${IONICE} cutadapt \
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
        }
    
    /*
    step 3: use umi_tools to extract umis from reads
    */

    process extract_umi {
        publishDir "${params.outdir}/umi", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
        cpus 16
        memory { 60.GB + (30.GB * task.attempt) }
        time '5h'
        tag "$libname"
        maxRetries 3
        queue 'standard'

        input:
        tuple val(libname), val(umi_len), path("${libname}.clip.r*.gz")

        output:
        tuple val(libname), val(umi_len), path("${libname}.umi.r1.gz"), path("${libname}.umi.log")

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

    /*
    step 4: extract umis and bcs and move to input file for starcode-umi
    */

    process prep_umibc {
        publishDir "${params.outdir}/parsed", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
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
    step 5: use starcode-umi to cluster umis for deduplication
    */

    process starcode_umi {
        publishDir "${params.outdir}/parsed", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
        tag "$libname"
        cpus 10
        memory { 40.GB + (40.GB * task.attempt) }
        maxRetries 3
        queue 'standard'

        input:
        tuple val(libname), val(umi_len), path(umibc)

        output:
        tuple val(libname), val(umi_len), path("${libname}.starumi")

        script:
        """
        ${IONICE} starcode-umi \
            --starcode-path /usr/local/bin/starcode \
            --umi-len ${umi_len} \
            --umi-d ${params.umi_dist} \
            ${umibc} > ${libname}.starumi
        """
    }

    /*
    step 6: deduplicate barcodes and perform barcode clustering with starcode
    */

    process starcode_umibc {
        publishDir "${params.outdir}/bc", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
        tag "$libname"
        cpus 10
        memory { 40.GB + (40.GB * task.attempt) }
        maxRetries 3
        queue 'standard'

        input:
        tuple val(libname), val(umi_len), path("${libname}.starumi")

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
    step 7: calculate read statistics across filtering and clustering steps
    */

    process read_stats_umibc {
        publishDir "${params.outdir}/postprocess", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
        cpus 1
        memory '20GB'
        time '30m'

        input:
        tuple val(libname), path(bc_cluster)

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

}

/*
step 6/8: filter clustered barcodes
*/

process bc_stats {
    publishDir "${params.outdir}/postprocess", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240619'
    cpus 1
    memory '20GB'
    time '30m'

    input:
    tuple val(libname), path(bc_cluster)

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
step 7/9: filter clustered barcodes
*/

process filter_bcgroups {
    publishDir "${params.outdir}/postprocess/filt_bc", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240619'
    cpus 4
    memory '24GB'
    tag "$libname"

    input:
    tuple val(libname), path(bc_cluster)

    output:
    tuple val(libname), path("${libname}.bc_cluster.filt.txt")

    script:
    """
    filter_starcode_bcgroups \
        --in_histo $bc_cluster \
        --out_histo ${libname}.bc_cluster.filt.txt \
        --out_summary ${libname}.summ.txt
    """
}

/*
step 8/10: make waterfall plots with clustered barcodes
*/

process waterfall_plot {
    publishDir "${params.outdir}/postprocess/plots", overwrite: true
    container 'library://tovar/general/mpra_utilscripts:20240619'
    cpus 1
    memory '10GB'
    tag "$libname"

    input:
    tuple val(libname), path(bc_cluster)

    output:
    tuple val(libname), path("${libname}.waterfall.png")

    script:
    """
    plot_tag_count_histos \
        --linCountHisto $bc_cluster \
        --lNames ${libname} \
        --annotMaxY \
        --outPlot ${libname}.waterfall.png
    """
}


if (params.check_samp_by_samp) {

    /*
    step 9/11: combine clustered barcodes and counts across all samples into a joint matrix
    */

    process joint_mtx {
        publishDir "${params.outdir}/postprocess/joint_mtx", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'
        cpus 4
        memory { 20.GB + (20.GB * task.attempt) }
        maxRetries 3
        time '2h'

        input:
        tuple val(libname), path(bc_cluster)

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
    step 10/12: matrix math to generate barcode statistics across samples
    */

    process sampbysamp_mtx {
        publishDir "${params.outdir}/postprocess/joint_mtx", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'

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
    step 11/13: plot stats generated in step 10/12
    */

    process sampbysamp_plot {
        publishDir "${params.outdir}/postprocess/plots", overwrite: true
        container 'library://tovar/general/mpra_utilscripts:20240619'

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

}

workflow {

libnames = params.libnames.keySet()

// create sample channels
if (params.ignore_umi){
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
    bc_out = bc_clip(clip_in) | prep_bc | starcode_bc
} else {
    bc_out = bc_clip_umi(clip_in) | extract_umi | prep_umibc | starcode_umi | starcode_umibc
}

bc_out.bc_cluster | bc_stats
bc_out.bc_cluster | filter_bcgroups
bc_out.bc_cluster | waterfall_plot

if (params.ignore_umi) {
    read_stats_bc(bc_out.bc_cluster)
} else {
    read_stats_umibc(bc_out.bc_cluster)
}

if (params.check_samp_by_samp) {
    joint_mtx(bc_out.bc_cluster) | sampbysamp_mtx | sampbysamp_plot
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
