#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.keyfile   = null
params.fastq_dir = null
params.genome    = null
params.vcf       = null
params.outdir    = null
params.workdir   = null
params.conda_env = 'genomicsV2'
params.only      = null   // v8: comma-separated Sample_IDs for a test run; null = all

// ── Validate required parameters ────────────────────────────
if (!params.keyfile)   error "Missing required param: --keyfile"
if (!params.fastq_dir) error "Missing required param: --fastq_dir"
if (!params.genome)    error "Missing required param: --genome"
if (!params.vcf)       error "Missing required param: --vcf"
if (!params.outdir)    error "Missing required param: --outdir"

// ── Parse key file and resolve FASTQ pairs into the channel ──
// Row index (1-based) is used as the Illumina S-number in FASTQ filenames
Channel
    .fromPath(params.keyfile)
    .splitCsv(header: true)
    .toList()
    .flatMap { rows ->
        rows.withIndex().collect { row, idx ->
            def sid = row.'Sample_ID'
            def s = idx + 1
            // V8: EVERY LANE. Up to v7 this matched _L001_ only, so lane 2 of a
            // two-lane NovaSeq X run was never aligned. Lanes are now globbed,
            // sorted, checked to pair up R1/R2, and concatenated into bwa.
            // BaseSpace puts FASTQs in a sample-named subdirectory; else flat.
            def pick = { rd ->
                def g = files("${params.fastq_dir}/${sid}*/${sid}_S${s}_L00*_${rd}_001.fastq.gz")
                if (!g) g = files("${params.fastq_dir}/${sid}_S${s}_L00*_${rd}_001.fastq.gz")
                g.sort { it.name }
            }
            def f1 = pick('R1'), f2 = pick('R2')
            if (!f1) error "FASTQ R1 not found for sample ${sid} (S${s}) in ${params.fastq_dir}"
            def lanes1 = f1.collect { (it.name =~ /_L00\d_/)[0] }
            def lanes2 = f2.collect { (it.name =~ /_L00\d_/)[0] }
            if (lanes1 != lanes2) error "R1/R2 lanes differ for ${sid}: ${lanes1} vs ${lanes2}"
            tuple(sid, f1, f2)
        }
    }
    // filter AFTER the S-number is taken from the row index, so a test run
    // on a few samples still finds the right files
    .filter { params.only ? (it[0] in params.only.tokenize(',')) : true }
    .set { samples_ch }

process ALIGN_SORT_COUNT {
    tag "${sample_id}"
    cpus 2
    memory '4 GB'
    time '4h'      // v8: both lanes, about twice v7's reads

    publishDir "${params.outdir}/bam_files",   mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/bam_files",   mode: 'copy', pattern: '*.lanes'
    publishDir "${params.outdir}/count_files",  mode: 'copy', pattern: '*.txt'

    input:
    tuple val(sample_id), path(r1), path(r2)      // lists: one file per lane
    path genome
    path genome_indices
    path vcf
    path vcf_index

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.txt"), emit: counts
    path "${sample_id}.lanes"

    script:
    """
    # record exactly which files went in (audit trail for the lane fix)
    ls -l ${r1} ${r2} > ${sample_id}.lanes

    # Align and sort; all lanes concatenated (gzip members concatenate cleanly)
    bwa mem -t ${task.cpus} \\
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' \\
        ${genome} <(cat ${r1}) <(cat ${r2}) | \\
    sambamba view --nthreads=${task.cpus} --sam-input --format=bam --with-header /dev/stdin | \\
    sambamba sort --nthreads=${task.cpus} --tmpdir=\${TMPDIR:-/tmp} --out=${sample_id}.bam /dev/stdin

    # Count allele-specific reads
    gatk ASEReadCounter \\
        --min-mapping-quality 20 \\
        -I ${sample_id}.bam \\
        -V ${vcf} \\
        -R ${genome} \\
        -O ${sample_id}.txt
    """
}

workflow {
    // Collect reference files as value channels
    genome_file    = file(params.genome)
    // Include both sacCer3.fasta.* (bwa indices, .fai) and sacCer3.dict (GATK)
    genome_dir     = file(params.genome).parent
    genome_base    = file(params.genome).baseName  // "sacCer3" from "sacCer3.fasta"
    genome_indices = Channel.fromPath(["${params.genome}.*", "${genome_dir}/${genome_base}.dict"]).collect()
    vcf_file       = file(params.vcf)
    vcf_index      = file("${params.vcf}.*")

    ALIGN_SORT_COUNT(samples_ch, genome_file, genome_indices, vcf_file, vcf_index)
}
