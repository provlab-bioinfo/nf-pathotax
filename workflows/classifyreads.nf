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
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pathotax_pipeline'

include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'

/* include {
    KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_ILLUMINA;
} from '../modules/nf-core/krakentools/combinekreports/main' */

include { SYLPH_PROFILE } from '../modules/local/sylph/profile/main'

include { SYLPHTAX_TAXPROF } from '../modules/local/sylphtax/taxprof/main'
include { SINGLEM_PIPE } from '../modules/local/singlem/pipe/main'
include { BRACKEN_BRACKEN } from '../modules/nf-core/bracken/bracken/main' 
include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../modules/nf-core/bracken/combinebrackenoutputs/main' 
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Function to select and validate database files
def selectDatabaseFiles(db_map, selected_databases) {
    //Get selected database keys
    def selected_keys
    if (selected_databases == 'all') {
        selected_keys = db_map.keySet()
    }
    else if (selected_databases instanceof List) {
        selected_keys = selected_databases
    }
    else {
        selected_keys = selected_databases.split(',')*.trim()
    }
    // Validate keys (optional but useful)
    def unknown_keys = selected_keys - db_map.keySet()
    if (unknown_keys) {
        error("Invalid database keys: ${unknown_keys.join(', ')}. Allowed: ${db_map.keySet().join(', ')}")
    }

    // Get the file paths
    def selected_dbs = selected_keys.collect { db_map[it] }

    // Convert to Nextflow file objects
    def db_paths = selected_dbs.collect { file(it) }

    return db_paths
}

workflow CLASSIFYREADS {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC(
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //classify
    if (!params.skip_kraken2) {
        KRAKEN2_KRAKEN2(
            ch_samplesheet,
            params.kraken2_db,
            false,
            true,
        )
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect { it[1] })

        BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, params.kraken2_db)
        
        ch_to_combine_bracken_report = BRACKEN_BRACKEN.out.reports
            .map{
                meta, report -> report
            }
            .collect()
            .map{
                reports -> tuple([id:"bracken_combined_report"], reports)
            }
        BRACKEN_COMBINEBRACKENOUTPUTS(ch_to_combine_bracken_report)
        ch_versions = ch_versions.mix(BRACKEN_COMBINEBRACKENOUTPUTS.out.versions)

    }
    /*
    sylph, a species-level metagenome profiler that estimates genome-to-metagenome containment 
    average nucleotide identity (ANI) through zero-inflated Poisson k-mer statistics, enabling 
    ANI-based taxa detection. 
    */
    if (!params.skip_sylph) {
        // Load DB map from params
        // Get database file paths using the function
        def db_paths = selectDatabaseFiles(params.sylph_db_files, params.sylph_databases)

        // Wrap as a channel
        ch_samplesheet
            .map { meta, reads -> tuple(meta, reads, db_paths) }
            .set { ch_sylph_inputs }


        SYLPH_PROFILE(ch_sylph_inputs)
        ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SYLPH_PROFILE.out.profile_out.collect { it[1] })

        // Get database file paths using the function
        def sylphtax_db_paths = selectDatabaseFiles(params.sylphtax_db_files, params.sylphtax_databases)
        // Wrap as a channel
        //ch_sylph_profile_output = SYLPH_PROFILE.out.profile_out

        SYLPH_PROFILE.out.profile_out
            .filter { meta, tsv -> tsv.size() > 0 && tsv.countLines() > 1 }
            .map { meta, tsv ->
                tuple(meta, tsv, sylphtax_db_paths)
            }
            .set { ch_sylphtax_inputs }

        SYLPHTAX_TAXPROF(ch_sylphtax_inputs)
        ch_versions = ch_versions.mix(SYLPHTAX_TAXPROF.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SYLPHTAX_TAXPROF.out.taxprof_output.collect { it[1] })
    }
    /*
    SingleM is a software suite which takes short read metagenomic data as input, 
    and estimates the relative abundance and per-base read coverage of Bacteria and 
    Archaea at each taxonomic level from domain to species. SingleM starts by matching 
    reads to highly conserved regions (’windows’) of 59 single copy marker genes 
    (22 Bacteria-specific, 24 Archaea-specific, 13 targeting both domains). 
    Importantly, reads are matched to these conserved gene windows by searching in 
    amino acid space, using DIAMOND BLASTX(Buchfink et al. 2021), maximising recruitment 
    of reads from divergent lineages. This is in contrast to other marker-based taxonomic 
    profilers, which map the nucleotide sequences of reads to markers directly (e.g. MetaPhlAn, mOTUs).
    */
    if (!params.skip_singlem && params.platform != "nanopore" ) {
        SINGLEM_PIPE(ch_samplesheet, params.singlem_db)
        ch_versions = ch_versions.mix(SINGLEM_PIPE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SINGLEM_PIPE.out.profile_out.collect { it[1] })
    }


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'pathotax_software_' + 'mqc_' + 'versions.yml',
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
