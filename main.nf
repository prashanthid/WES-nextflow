#!/usr/bin/env nextflow
/*
 * ============================================================================
 * Somatic WES pipeline: BQSR -> Mutect2 + Strelka2 (Manta-assisted) -> Funcotator
 * ============================================================================
 * Input: paired tumor/normal BAMs, BWA-MEM aligned, Picard-deduplicated
 *        (NOT yet base-recalibrated)
 *
 * Usage:
 *   nextflow run main.nf -params-file params.yaml -resume
 * ============================================================================
 */

// ---------------------------------------------------------------------------
// Parameters (override via -params-file or --flag on the CLI)
// ---------------------------------------------------------------------------
params.samplesheet              = null   // CSV: patient_id,tumor_bam,tumor_bai,normal_bam,normal_bai
params.outdir                   = 'results'

// Reference
params.fasta                    = null
params.fasta_fai                = null
params.fasta_dict               = null
params.intervals                = null   // WES bait/target interval_list (required)

// Manta/Strelka2 need a SEPARATE bgzip+tabix'd BED for --callRegions.
// Prepare with: sort -k1,1 -k2,2n targets.bed | bgzip -c > targets.bed.gz && tabix -p bed targets.bed.gz
params.call_regions_bed         = null   // optional; if unset, Manta/Strelka2 run --exome without --callRegions

// BQSR known sites (repeatable: dbSNP + Mills + 1000G indels etc.)
// Each VCF must have a matching .tbi index sitting next to it.
params.known_sites              = []

// Mutect2 resources
params.germline_resource        = null   // gnomAD af-only VCF (+ .tbi), required for contamination estimation
params.pon                      = null   // panel of normals VCF (+ .tbi), optional

// Funcotator
params.funcotator_sources       = null   // path to Funcotator data source directory
params.funcotator_ref_version   = 'hg38' // hg19 | hg38
params.funcotator_output_format = 'MAF'  // MAF | VCF

// ---------------------------------------------------------------------------
// PROCESS: BaseRecalibrator
// ---------------------------------------------------------------------------
process BASE_RECALIBRATOR {
    tag "${patient_id}_${sample_type}"
    publishDir { "${params.outdir}/bqsr/${patient_id}" }, mode: 'copy', pattern: "*.table"

    input:
    tuple val(patient_id), val(sample_type), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path known_sites_all
    path intervals

    output:
    tuple val(patient_id), val(sample_type), path(bam), path(bai), path("${patient_id}_${sample_type}.recal.table"), emit: recal_table

    script:
    def known_sites_args = known_sites_all
        .findAll { !(it.name.endsWith('.tbi') || it.name.endsWith('.idx')) }
        .collect { "--known-sites ${it}" }
        .join(' ')
    """
    gatk BaseRecalibrator \\
        -R ${fasta} \\
        -I ${bam} \\
        ${known_sites_args} \\
        -L ${intervals} \\
        -O ${patient_id}_${sample_type}.recal.table
    """
}

// ---------------------------------------------------------------------------
// PROCESS: ApplyBQSR
// ---------------------------------------------------------------------------
process APPLY_BQSR {
    tag "${patient_id}_${sample_type}"
    publishDir { "${params.outdir}/bqsr/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), val(sample_type), path(bam), path(bai), path(recal_table)
    path fasta
    path fai
    path dict
    path intervals

    output:
    tuple val(patient_id), val(sample_type), path("${patient_id}_${sample_type}.recal.bam"), path("${patient_id}_${sample_type}.recal.bai"), emit: recal_bam

    script:
    """
    gatk ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -L ${intervals} \\
        -O ${patient_id}_${sample_type}.recal.bam
    """
}

// ---------------------------------------------------------------------------
// PROCESS: Mutect2 (tumor vs normal)
// ---------------------------------------------------------------------------
process MUTECT2 {
    tag "${patient_id}"
    errorStrategy 'retry'
    maxRetries 2
    publishDir { "${params.outdir}/mutect2/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(tumor_sm), val(normal_sm)
    path fasta
    path fai
    path dict
    path intervals
    path germline_resource_files
    path pon_files

    output:
    tuple val(patient_id),
          path("${patient_id}.mutect2.unfiltered.vcf.gz"),
          path("${patient_id}.mutect2.unfiltered.vcf.gz.tbi"),
          path("${patient_id}.mutect2.unfiltered.vcf.gz.stats"),
          path("${patient_id}.f1r2.tar.gz"), emit: mutect2_raw

    script:
    def germline_vcf = germline_resource_files.find { it.name.endsWith('.vcf.gz') || it.name.endsWith('.vcf') }
    def germline_arg = germline_vcf ? "--germline-resource ${germline_vcf}" : ''
    def pon_vcf       = pon_files.find { it.name.endsWith('.vcf.gz') || it.name.endsWith('.vcf') }
    def pon_arg       = pon_vcf ? "--panel-of-normals ${pon_vcf}" : ''
    """
    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${tumor_bam} -tumor ${tumor_sm} \\
        -I ${normal_bam} -normal ${normal_sm} \\
        ${germline_arg} \\
        ${pon_arg} \\
        -L ${intervals} \\
        --f1r2-tar-gz ${patient_id}.f1r2.tar.gz \\
        -O ${patient_id}.mutect2.unfiltered.vcf.gz
    """
}

// ---------------------------------------------------------------------------
// PROCESS: LearnReadOrientationModel
// ---------------------------------------------------------------------------
process LEARN_READ_ORIENTATION_MODEL {
    tag "${patient_id}"
    publishDir { "${params.outdir}/mutect2/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), path(f1r2)

    output:
    tuple val(patient_id), path("${patient_id}.read-orientation-model.tar.gz"), emit: rom

    script:
    """
    gatk LearnReadOrientationModel \\
        -I ${f1r2} \\
        -O ${patient_id}.read-orientation-model.tar.gz
    """
}

// ---------------------------------------------------------------------------
// PROCESS: GetPileupSummaries (run per tumor/normal)
// ---------------------------------------------------------------------------
process GET_PILEUP_SUMMARIES {
    tag "${patient_id}_${sample_type}"
    publishDir { "${params.outdir}/mutect2/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), val(sample_type), path(bam), path(bai)
    path germline_resource_files
    path intervals

    output:
    tuple val(patient_id), val(sample_type), path("${patient_id}_${sample_type}.pileups.table"), emit: pileups

    script:
    def germline_vcf = germline_resource_files.find { it.name.endsWith('.vcf.gz') || it.name.endsWith('.vcf') }
    """
    gatk GetPileupSummaries \\
        -I ${bam} \\
        -V ${germline_vcf} \\
        -L ${intervals} \\
        -O ${patient_id}_${sample_type}.pileups.table
    """
}

// ---------------------------------------------------------------------------
// PROCESS: CalculateContamination
// ---------------------------------------------------------------------------
process CALCULATE_CONTAMINATION {
    tag "${patient_id}"
    publishDir { "${params.outdir}/mutect2/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), path(tumor_pileups), path(normal_pileups)

    output:
    tuple val(patient_id),
          path("${patient_id}.contamination.table"),
          path("${patient_id}.segments.table"), emit: contamination

    script:
    """
    gatk CalculateContamination \\
        -I ${tumor_pileups} \\
        -matched ${normal_pileups} \\
        --tumor-segmentation ${patient_id}.segments.table \\
        -O ${patient_id}.contamination.table
    """
}

// ---------------------------------------------------------------------------
// PROCESS: FilterMutectCalls
// ---------------------------------------------------------------------------
process FILTER_MUTECT_CALLS {
    tag "${patient_id}"
    publishDir { "${params.outdir}/mutect2/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), path(unfiltered_vcf), path(unfiltered_tbi), path(stats),
          path(rom), path(contamination_table), path(segments_table)
    path fasta
    path fai
    path dict

    output:
    tuple val(patient_id),
          path("${patient_id}.mutect2.filtered.vcf.gz"),
          path("${patient_id}.mutect2.filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    """
    gatk FilterMutectCalls \\
        -R ${fasta} \\
        -V ${unfiltered_vcf} \\
        --contamination-table ${contamination_table} \\
        --tumor-segmentation ${segments_table} \\
        --ob-priors ${rom} \\
        -O ${patient_id}.mutect2.filtered.vcf.gz
    """
}

// ---------------------------------------------------------------------------
// PROCESS: Manta (produces candidate small indels used to boost Strelka2 sensitivity)
// ---------------------------------------------------------------------------
process MANTA {
    tag "${patient_id}"
    publishDir { "${params.outdir}/manta/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fai
    path dict
    path call_regions_files

    output:
    tuple val(patient_id),
          path("results/variants/candidateSmallIndels.vcf.gz"),
          path("results/variants/candidateSmallIndels.vcf.gz.tbi"), emit: manta_indels

    script:
    def bed = call_regions_files.find { it.name.endsWith('.bed.gz') }
    def bed_arg = bed ? "--exome --callRegions ${bed}" : "--exome"
    """
    configManta.py \\
        --normalBam ${normal_bam} \\
        --tumorBam ${tumor_bam} \\
        --referenceFasta ${fasta} \\
        ${bed_arg} \\
        --runDir .

    ./runWorkflow.py -m local -j ${task.cpus}
    """
}

// ---------------------------------------------------------------------------
// PROCESS: Strelka2 (somatic workflow)
// ---------------------------------------------------------------------------
process STRELKA2 {
    tag "${patient_id}"
    publishDir { "${params.outdir}/strelka2/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai),
          path(manta_indels), path(manta_indels_tbi)
    path fasta
    path fai
    path dict
    path call_regions_files

    output:
    tuple val(patient_id),
          path("results/variants/somatic.snvs.vcf.gz"),
          path("results/variants/somatic.snvs.vcf.gz.tbi"),
          path("results/variants/somatic.indels.vcf.gz"),
          path("results/variants/somatic.indels.vcf.gz.tbi"), emit: strelka_vcfs

    script:
    def bed = call_regions_files.find { it.name.endsWith('.bed.gz') }
    def bed_arg = bed ? "--exome --callRegions ${bed}" : "--exome"
    """
    configureStrelkaSomaticWorkflow.py \\
        --normalBam ${normal_bam} \\
        --tumorBam ${tumor_bam} \\
        --referenceFasta ${fasta} \\
        --indelCandidates ${manta_indels} \\
        ${bed_arg} \\
        --runDir .

    ./runWorkflow.py -m local -j ${task.cpus}
    """
}

// ---------------------------------------------------------------------------
// PROCESS: merge Strelka2 SNV + indel VCFs into one sorted, bgzipped VCF
// ---------------------------------------------------------------------------
process MERGE_STRELKA_VCFS {
    tag "${patient_id}"
    publishDir { "${params.outdir}/strelka2/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), path(snvs), path(snvs_tbi), path(indels), path(indels_tbi)

    output:
    tuple val(patient_id),
          path("${patient_id}.strelka2.somatic.merged.vcf.gz"),
          path("${patient_id}.strelka2.somatic.merged.vcf.gz.tbi"), emit: merged_vcf

    script:
    """
    bcftools concat -a ${snvs} ${indels} -Oz -o tmp.vcf.gz
    bcftools sort tmp.vcf.gz -Oz -o ${patient_id}.strelka2.somatic.merged.vcf.gz
    tabix -p vcf ${patient_id}.strelka2.somatic.merged.vcf.gz
    rm tmp.vcf.gz
    """
}

// ---------------------------------------------------------------------------
// PROCESS: Consensus VCF (Mutect2 ∩ Strelka2, PASS-only variants called by both)
// ---------------------------------------------------------------------------
// Restricts each caller's PASS calls, then intersects with bcftools isec.
// -n+2 keeps only sites present in both inputs; -w1 writes those sites using
// the Mutect2 record (its INFO/FORMAT fields tend to be more useful downstream,
// e.g. AF from the tumor sample). Site-level intersection only — this does not
// reconcile differing REF/ALT representations at indel sites (e.g. left-alignment
// or multi-allelic splitting), so consider `bcftools norm -m-` on both inputs
// upstream if your calls include many complex/multi-allelic indels.
process CONSENSUS_VCF {
    tag "${patient_id}"
    publishDir { "${params.outdir}/consensus/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), path(mutect2_vcf), path(mutect2_tbi), path(strelka2_vcf), path(strelka2_tbi)

    output:
    tuple val(patient_id),
          path("${patient_id}.consensus.vcf.gz"),
          path("${patient_id}.consensus.vcf.gz.tbi"), emit: consensus_vcf
    path("${patient_id}.isec_summary.txt"), emit: isec_summary

    script:
    """
    bcftools view -f PASS ${mutect2_vcf}  -Oz -o mutect2.pass.vcf.gz
    tabix -p vcf mutect2.pass.vcf.gz

    bcftools view -f PASS ${strelka2_vcf} -Oz -o strelka2.pass.vcf.gz
    tabix -p vcf strelka2.pass.vcf.gz

    # -n+2 : sites present in >=2 of the inputs (i.e. both, since there are only 2)
    # -w1  : write the surviving sites using file #1 (Mutect2) records
    bcftools isec -n+2 -w1 -p ${patient_id}_isec -Oz mutect2.pass.vcf.gz strelka2.pass.vcf.gz

    mv ${patient_id}_isec/0000.vcf.gz     ${patient_id}.consensus.vcf.gz
    mv ${patient_id}_isec/0000.vcf.gz.tbi ${patient_id}.consensus.vcf.gz.tbi

    {
        echo "patient_id: ${patient_id}"
        echo "mutect2 PASS sites:  \$(zcat mutect2.pass.vcf.gz  | grep -vc '^#')"
        echo "strelka2 PASS sites: \$(zcat strelka2.pass.vcf.gz | grep -vc '^#')"
        echo "consensus sites:     \$(zcat ${patient_id}.consensus.vcf.gz | grep -vc '^#')"
    } > ${patient_id}.isec_summary.txt
    """
}

// ---------------------------------------------------------------------------
// PROCESS: Funcotator (generic - reused for Mutect2, Strelka2, and consensus outputs)
// ---------------------------------------------------------------------------
process FUNCOTATOR {
    tag "${patient_id}_${caller}"
    publishDir { "${params.outdir}/funcotator/${patient_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), val(caller), path(vcf), path(vcf_tbi)
    path fasta
    path fai
    path dict
    path intervals
    path funcotator_sources

    output:
    tuple val(patient_id), val(caller), path("${patient_id}.${caller}.funcotated.*"), emit: annotated

    script:
    def out_ext  = params.funcotator_output_format == 'MAF' ? 'maf' : 'vcf'
    def out_name = "${patient_id}.${caller}.funcotated.${out_ext}"
    """
    gatk Funcotator \\
        -R ${fasta} \\
        -V ${vcf} \\
        -L ${intervals} \\
        --data-sources-path ${funcotator_sources} \\
        --output-file-format ${params.funcotator_output_format} \\
        --ref-version ${params.funcotator_ref_version} \\
        -O ${out_name}
    """
}

// ---------------------------------------------------------------------------
// WORKFLOW
// ---------------------------------------------------------------------------
workflow {

    // ---- validate required params ----
    if( !params.samplesheet )        error "ERROR: --samplesheet is required"
    if( !params.fasta )              error "ERROR: --fasta is required"
    if( !params.intervals )          error "ERROR: --intervals is required"
    if( !params.germline_resource )  error "ERROR: --germline_resource is required"
    if( !params.funcotator_sources ) error "ERROR: --funcotator_sources is required"

    // ---- reference value channels (broadcast to every process invocation) ----
    def fasta_path = file(params.fasta)
    def fai_path   = params.fasta_fai  ? file(params.fasta_fai)  : file("${params.fasta}.fai")
    def dict_path  = params.fasta_dict ? file(params.fasta_dict) : file(params.fasta.toString().replaceAll(/\.fa(sta)?$/, '.dict'))

    ref_fasta_ch = Channel.value(fasta_path)
    ref_fai_ch   = Channel.value(fai_path)
    ref_dict_ch  = Channel.value(dict_path)
    intervals_ch = Channel.value(file(params.intervals))

    known_sites_ch = Channel.value(
        params.known_sites.collectMany { f ->
            def vcf = file(f)
            def tbi = file("${f}.tbi")
            def idx = file("${f}.idx")
            def idxFile = tbi.exists() ? tbi : idx
            [ vcf, idxFile ]
        }
    )

    def germline_vcf = file(params.germline_resource)
    def germline_tbi = file("${params.germline_resource}.tbi")
    def germline_idx = file("${params.germline_resource}.idx")
    germline_resource_ch = Channel.value(
        [ germline_vcf, (germline_tbi.exists() ? germline_tbi : germline_idx) ]
    )

    pon_ch = Channel.value(
        params.pon ? [ file(params.pon), file("${params.pon}.tbi") ] : []
    )

    call_regions_ch = Channel.value(
        params.call_regions_bed ? [ file(params.call_regions_bed), file("${params.call_regions_bed}.tbi") ] : []
    )

    funcotator_sources_ch = Channel.value(file(params.funcotator_sources))

    // ---- sample sheet ----
    sample_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.patient_id,
                file(row.tumor_bam),  file(row.tumor_bai),
                file(row.normal_bam), file(row.normal_bai),
                row.tumor_sample_name,
                row.normal_sample_name
            )
        }

    bqsr_input_ch = sample_ch.flatMap { patient_id, tbam, tbai, nbam, nbai, tumor_sm, normal_sm ->
        [
            tuple(patient_id, 'tumor',  tbam, tbai),
            tuple(patient_id, 'normal', nbam, nbai)
        ]
    }

    sample_names_ch = sample_ch.map { patient_id, tbam, tbai, nbam, nbai, tumor_sm, normal_sm ->
        tuple(patient_id, tumor_sm, normal_sm)
    }

    // 1) BQSR on every tumor/normal BAM independently
    BASE_RECALIBRATOR(bqsr_input_ch, ref_fasta_ch, ref_fai_ch, ref_dict_ch, known_sites_ch, intervals_ch)
    APPLY_BQSR(BASE_RECALIBRATOR.out.recal_table, ref_fasta_ch, ref_fai_ch, ref_dict_ch, intervals_ch)

    // Split recalibrated BAMs back into tumor / normal channels and re-pair by patient
    recal_split = APPLY_BQSR.out.recal_bam.branch {
        tumor:  it[1] == 'tumor'
        normal: it[1] == 'normal'
    }

    tumor_recal_ch  = recal_split.tumor.map  { pid, type, bam, bai -> tuple(pid, bam, bai) }
    normal_recal_ch = recal_split.normal.map { pid, type, bam, bai -> tuple(pid, bam, bai) }

    paired_recal_ch = tumor_recal_ch
        .join(normal_recal_ch)   // joins on patient_id
        .map { pid, tbam, tbai, nbam, nbai -> tuple(pid, tbam, tbai, nbam, nbai) }

    // Sample names for -tumor/-normal flags: taken directly from the samplesheet's
    // tumor_sample_name / normal_sample_name columns (must exactly match each BAM's @RG SM: tag)
    mutect2_input_ch = paired_recal_ch
        .join(sample_names_ch)
        .map { pid, tbam, tbai, nbam, nbai, tumor_sm, normal_sm ->
            tuple(pid, tbam, tbai, nbam, nbai, tumor_sm, normal_sm)
        }

    // 2) Mutect2 branch
    MUTECT2(mutect2_input_ch, ref_fasta_ch, ref_fai_ch, ref_dict_ch, intervals_ch, germline_resource_ch, pon_ch)
    LEARN_READ_ORIENTATION_MODEL(MUTECT2.out.mutect2_raw.map { pid, vcf, tbi, stats, f1r2 -> tuple(pid, f1r2) })

    GET_PILEUP_SUMMARIES(APPLY_BQSR.out.recal_bam, germline_resource_ch, intervals_ch)

    pileup_split = GET_PILEUP_SUMMARIES.out.pileups.branch {
        tumor:  it[1] == 'tumor'
        normal: it[1] == 'normal'
    }

    tumor_pileup_ch  = pileup_split.tumor.map  { pid, type, table -> tuple(pid, table) }
    normal_pileup_ch = pileup_split.normal.map { pid, type, table -> tuple(pid, table) }

    contam_input_ch = tumor_pileup_ch.join(normal_pileup_ch)
    CALCULATE_CONTAMINATION(contam_input_ch)

    filter_input_ch = MUTECT2.out.mutect2_raw
        .map { pid, vcf, tbi, stats, f1r2 -> tuple(pid, vcf, tbi, stats) }
        .join(LEARN_READ_ORIENTATION_MODEL.out.rom)
        .join(CALCULATE_CONTAMINATION.out.contamination)

    FILTER_MUTECT_CALLS(filter_input_ch, ref_fasta_ch, ref_fai_ch, ref_dict_ch)

    // 3) Manta + Strelka2 branch
    MANTA(paired_recal_ch, ref_fasta_ch, ref_fai_ch, ref_dict_ch, call_regions_ch)

    strelka_input_ch = paired_recal_ch.join(MANTA.out.manta_indels)
    STRELKA2(strelka_input_ch, ref_fasta_ch, ref_fai_ch, ref_dict_ch, call_regions_ch)
    MERGE_STRELKA_VCFS(STRELKA2.out.strelka_vcfs)

    // 4) Consensus: intersect PASS-only Mutect2 + Strelka2 calls per patient
    consensus_input_ch = FILTER_MUTECT_CALLS.out.filtered_vcf
        .join(MERGE_STRELKA_VCFS.out.merged_vcf)   // joins on patient_id
        .map { pid, mvcf, mtbi, svcf, stbi -> tuple(pid, mvcf, mtbi, svcf, stbi) }

    CONSENSUS_VCF(consensus_input_ch)

    // 5) Funcotator on filtered Mutect2, merged Strelka2, AND the consensus VCF
    funcotator_input_ch = FILTER_MUTECT_CALLS.out.filtered_vcf
        .map { pid, vcf, tbi -> tuple(pid, 'mutect2', vcf, tbi) }
        .mix(
            MERGE_STRELKA_VCFS.out.merged_vcf.map { pid, vcf, tbi -> tuple(pid, 'strelka2', vcf, tbi) }
        )
        .mix(
            CONSENSUS_VCF.out.consensus_vcf.map { pid, vcf, tbi -> tuple(pid, 'consensus', vcf, tbi) }
        )

    FUNCOTATOR(funcotator_input_ch, ref_fasta_ch, ref_fai_ch, ref_dict_ch, intervals_ch, funcotator_sources_ch)
}
