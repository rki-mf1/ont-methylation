#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// terminal prints
println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "

// error codes
if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if ( !workflow.revision ) { 
  println "\033[0;33mWARNING: It is recommended to use a stable release version via -r." 
  println "Use 'nextflow info valegale/ONT_methylation' to check for available release versions.\033[0m\n"
}
// help
if (params.help) { exit 0, helpMSG() }


// input channels conditioned on flow (main, annotation, dmr) 
if (params.flow == 'main') {

    if (!params.bam)   { error "❌ --bam is required for --flow main" }
    if (!params.fasta) { error "❌ --fasta is required for --flow main" }

    if (params.fasta && params.list) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
    } else if (params.fasta) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    }

    if (params.bam && params.list) { bam_input_ch = Channel
        .fromPath( params.bam, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
    } else if (params.bam) { bam_input_ch = Channel
        .fromPath( params.bam, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    }

    if (params.meta) {
        if (!params.bin_folder) { error "--bin_folder must be provided when using --meta" }
        bins = Channel
            .fromPath("${params.bin_folder}/*.{fasta,fa}", checkIfExists: true)
            .ifEmpty { error("No bin FASTA files found in folder: ${params.bin_folder}") }
        fasta_input_ch.count().map { cnt ->
            if (cnt != 1) error "Exactly one FASTA file must be provided when using --meta. Found: ${cnt}"
        }
        bam_input_ch.count().map { cnt ->
            if (cnt != 1) error "Exactly one BAM file must be provided when using --meta. Found: ${cnt}"
        }
    }

} else if (params.flow == 'annotation') {

    if (!params.bed)   { error "❌ --bed is required for --flow annotation" }
    if (!params.fasta) { error "❌ --fasta is required for --flow annotation" }
    if (!params.gff3)  { error "❌ --gff3 is required for --flow annotation" }

    bed_input_ch = Channel
        .fromPath(params.bed, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }

    fasta_input_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }

    gff3_input_ch = Channel
        .fromPath(params.gff3, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }

} else if (params.flow == 'dmr') {

    if (!params.bam)   { error "❌ --bam is required for --flow dmr" }
    if (!params.fasta) { error "❌ --fasta is required for --flow dmr" }
    if (!params.gff3)  { error "❌ --gff3 is required for --flow dmr" }
    // dmr input channels WIP

}

// include workflows
include { MAIN_FLOW }       from './workflows/main_flow'
include { ANNOTATION_FLOW } from './workflows/annotation_flow'

params.flow = 'main'

workflow {
    if (params.flow == 'main') {
        MAIN_FLOW(bam_input_ch, fasta_input_ch)
    } else if (params.flow == 'annotation') {
        ANNOTATION_FLOW(
            bed_input_ch,
            fasta_input_ch,
            gff3_input_ch
        )
    } else if (params.flow == 'dmr') {
        //work in progress
        DMR_FLOW()
    } else {
        error "❌ Unknown --flow '${params.flow}'. Valid options: main, annotation, dmr"
    }
}

// --help
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_red = "\033[0;31m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    Nextflow Pipeline for Methylated Motif Extraction and Statistical Analysis from ONT bacterial data.

    ${c_yellow}Usage example:${c_reset}
    nextflow run valegale/ONT_methylation -r 0.0.1 --fasta '*.fasta' --bam '*.bam' 

    Use the following commands to check for latest pipeline versions:
    
    nextflow pull valegale/ONT_methylation
    nextflow info valegale/ONT_methylation

    ${c_yellow}Input${c_reset}
    ${c_green} --fasta ${c_reset}           '*.fasta'       -> one genome/assembly per file
    ${c_green} --bam ${c_reset}             '*.bam'         -> one sorted BAM matching one reference FASTA

    ${c_dim}  change above input to csv:${c_reset} ${c_green}--list ${c_reset}

    ${c_yellow}IMPORTANT:${c_reset} Unless ${c_green}--list${c_reset} is used, the ${c_yellow}basename${c_reset} of the FASTA and BAM files must match
    (e.g., sample1.fasta <-> sample1.bam). 

    ${c_yellow}General Options:${c_reset}
    --cores             Max cores per process for local use [default: $params.cores]
    --max_cores         Max cores (in total) for local use [default: $params.max_cores]
    --memory            Max memory for local use [default: $params.memory]
    --outdir            Name of the result folder [default: $params.outdir]

    ${c_yellow}Additional Options:${c_reset}${c_reset}
    --filter_threshold_modkit             Filter threshold for modkit [default: $params.filter_threshold_modkit]
    --automatic_threshold_modkit          Enable automatic estimation of the filter threshold by modkit.
                                          When true, modkit will determine an optimal threshold from the data and the value of --filter_threshold_modkit will be ignored.
                                          [default: $params.automatic_threshold_modkit]
    --percent_cutoff_modification_table   Minimum methylation percentage required for genome positions to be reported in the modification tables [default: $params.percent_cutoff_modification_table].
 
    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)
    -resume                  resume a previous calculation w/o recalculating everything (needs the same run command and work dir!)

    ${c_yellow}Caching:${c_reset}
    --singularityCacheDir   Location for storing the Singularity images [default: $params.singularityCacheDir]
    -w                      Working directory for all intermediate results [default: work] 

    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}docker${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
    
    ${c_blue}Engines${c_reset} (choose one):
      docker
      singularity
    
    Per default: -profile local,docker is executed (-profile standard).
    
    ${c_reset}
    """.stripIndent()
}
