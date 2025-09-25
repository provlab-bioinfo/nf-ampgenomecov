/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ampgenomecov_pipeline'
include { COLLAPSE_PRIMERS } from '../modules/local/collapse_primers'
include { MOSDEPTH as MOSDEPTH_GENOME } from '../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_AMPLICON } from '../modules/nf-core/mosdepth/main'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME } from '../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON } from '../modules/local/plot_mosdepth_regions'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AMPGENOMECOV {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (!params.ref || !file(params.ref).exists()) {
        exit(1, "ERROR: Reference file '${params.ref}' not found!")
    }

    ch_ref = Channel.fromPath(params.ref)
        .map { ref ->
            def meta2 = ref.baseName
            [meta2, ref]
        }
    if (!params.primer_bed || !file(params.primer_bed).exists()) {
        exit(1, "ERROR: primer_bed file '${params.primer_bed}' not found!")
    }
    ch_primer_bed = Channel.fromPath(params.primer_bed)


    // Separate tuples with and without bai
    ch_bam_bai = ch_samplesheet.filter { meta, bam, bai -> bai != null }

    //make bai if it is not existing
    ch_only_bam = ch_samplesheet
        .filter { meta, bam, bai -> bai == null }
        .map { meta, bam, bai ->
            [meta, bam]
        }

    SAMTOOLS_INDEX(ch_only_bam)

    // Merge both streams
    ch_bam_bai = ch_bam_bai.mix(ch_only_bam.join(SAMTOOLS_INDEX.out.bai))

    //prepare intput to MOSDEPTH_GENOME
    ch_bam_bai.combine(ch_ref)
        .multiMap { meta, bam, bai, meta2, ref ->
            bam_bai_bed: [meta, bam, bai, []]
            ref: [meta2, ref]
        }
        .set {
            ch_input_genome_depth
        }

    MOSDEPTH_GENOME(ch_input_genome_depth.bam_bai_bed, ch_input_genome_depth.ref)
    ch_mosdepth_multiqc = MOSDEPTH_GENOME.out.global_txt
    ch_multiqc_files = ch_multiqc_files.mix(ch_mosdepth_multiqc.collect { it[1] })
    ch_versions = ch_versions.mix(MOSDEPTH_GENOME.out.versions.first().ifEmpty(null))

    PLOT_MOSDEPTH_REGIONS_GENOME(
        MOSDEPTH_GENOME.out.regions_bed.collect { it[1] }
    )
    ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME.out.versions)

    // for amplicon
    COLLAPSE_PRIMERS(ch_primer_bed,params.primer_left_suffix,params.primer_right_suffix)
    ch_primer_collapsed_bed = COLLAPSE_PRIMERS.out.bed
    ch_versions = ch_versions.mix(COLLAPSE_PRIMERS.out.versions)

    ch_bam_bai_bed = ch_bam_bai.combine(ch_primer_collapsed_bed)
    ch_bam_bai_bed
        .combine(ch_ref)
        .multiMap { it ->
            bam_bai_bed: [it[0], it[1], it[2], it[3]]
            ref: [it[4], it[5]]
        }
        .set {
            ch_input_amp_depth
        }

    MOSDEPTH_AMPLICON(ch_input_amp_depth.bam_bai_bed, ch_input_amp_depth.ref)
    
    ch_versions = ch_versions.mix(MOSDEPTH_AMPLICON.out.versions.first().ifEmpty(null))
    PLOT_MOSDEPTH_REGIONS_AMPLICON(
        MOSDEPTH_AMPLICON.out.regions_bed.collect { it[1] }
    )
    ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AMPLICON.out.versions)
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'ampgenomecov_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions = ch_versions // channel: [ path(versions.yml) ]
}
