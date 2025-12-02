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


// input
// genomes fasta input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  //.view() 
  }
  else if (params.fasta) { fasta_input_ch = Channel
    .fromPath( params.fasta, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
}

// BAM input & --list support
if (params.bam && params.list) { bam_input_ch = Channel
  .fromPath( params.bam, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  //.view() 
  }
  else if (params.bam) { bam_input_ch = Channel
    .fromPath( params.bam, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
}

// bins input & --meta support
if (params.meta) {
      
    if (!params.bin_folder) {
        error "--bin_folder must be provided when using --meta"
    }

    bins = Channel
        .fromPath("${params.bin_folder}/*.{fasta,fa}", checkIfExists: true)
        .ifEmpty { error("No bin FASTA files found in folder: ${params.bin_folder}") }
    
    fasta_input_ch.count()
        .map { cnt ->
            if (cnt != 1) error "Exactly one FASTA file must be provided when using -meta. Found: ${cnt}"
        }

    bam_input_ch.count()
        .map { cnt ->
            if (cnt != 1) error "Exactly one BAM file must be provided when using -meta. Found: ${cnt}"
        }
}

// load modules
include { bam2fastq; zipfastq; minimap2; split_bam_by_bin } from './modules/map_index_bam.nf'
include { modkit_pileup; modkit_pileup_bedgraphs; modkit_find_motifs; custom_bedgraphs; publish_results_meta; publish_results_motifs_meta; publish_results; publish_results_motifs} from './modules/modkit.nf'
include { compute_statistics } from './modules/statistics.nf'


// main workflow
workflow {
    fastq_files = bam2fastq(bam_input_ch)
    zipfastq(fastq_files)

    // combine the fastq files with the reference fasta files 
    fastq_ref_pairs = fastq_files.join(fasta_input_ch)
    fastq_ref_pairs
    .ifEmpty { 
        error "❌ No matching FASTA found for BAM samples (example: sample1.bam and sample1.fasta). " +
              "Check naming or provide --list mapping file with explicit sample-to-genome pairs." 
    }
    mapped_bams = minimap2(fastq_ref_pairs)

    // if meta mode is on, split bams by bins first 
    if (params.meta) {
        bam_bin_pairs = mapped_bams.combine(bins)
        filtered_bams = split_bam_by_bin(bam_bin_pairs)
        
        bed_file = modkit_pileup(filtered_bams)  
        pileup_bedgraphs_ch = modkit_pileup_bedgraphs(filtered_bams)
    } else {
        bed_file = modkit_pileup(mapped_bams)  
        pileup_bedgraphs_ch = modkit_pileup_bedgraphs(mapped_bams)
    }

    motifs_ch = modkit_find_motifs(bed_file)
    custom_bedgraphs_ch = custom_bedgraphs(bed_file)
    statistics_ch = compute_statistics(bed_file)

    publish_input = bed_file.join(pileup_bedgraphs_ch)
                            .join(custom_bedgraphs_ch)
                            .join(statistics_ch) 

    publish_motifs_input = bed_file.join(motifs_ch)

    if (params.meta) {
        publish_results_meta(publish_input)
        publish_results_motifs_meta(publish_motifs_input)
    } else {
        publish_results(publish_input)
        publish_results_motifs(publish_motifs_input)
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

    Nextflow Pipeline for Methylated Motif Extraction and Statistical Analysis from ONT data.

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
