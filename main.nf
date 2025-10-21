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
  println "\033[0;33mWARNING: It is recommended to use a stable relese version via -r." 
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

// Define the list of BAM files and reference files
/*    params.input_files = [
        [ file('../run_dorado_v5/barcode01/calling_barcode1.bam'), file('../run_dorado_v5/barcode01/polished_assembly_01.fasta')],
        [ file('../run_dorado_v5/barcode02/calling_barcode02.bam'), file('../run_dorado_v5/barcode02/polished_assembly_02.fasta')],
        [ file('../run_dorado_v5/barcode03/calling_barcode03.bam'), file('../run_dorado_v5/barcode03/polished_assembly_03.fasta')],
        [ file('../run_dorado_v5/barcode04/calling_barcode4.bam'), file('../run_dorado_v5/barcode04/polished_assembly_04.fasta')],
        [ file('../run_dorado_v5/barcode05/calling_barcode5.bam'), file('../run_dorado_v5/barcode05/polished_assembly_05.fasta')],
        [ file('../run_dorado_v5/barcode06/calling_barcode6.bam'), file('../run_dorado_v5/barcode06/polished_assembly_06.fasta')],
        [ file('../run_dorado_v5/barcode07/calling_barcode7.bam'), file('../run_dorado_v5/barcode07/polished_assembly_07.fasta')],
        [ file('../run_dorado_v5/barcode08/calling_barcode08.bam'), file('../run_dorado_v5/barcode08/polished_assembly_08.fasta')],
        [ file('../run_dorado_v5/barcode09/calling_barcode09.bam'), file('../run_dorado_v5/barcode09/polished_assembly_09.fasta')],
        [ file('../run_dorado_v5/barcode10/calling_barcode10.bam'), file('../run_dorado_v5/barcode10/polished_assembly_10.fasta')],
        [ file('../run_dorado_v5/barcode11/calling_barcode11.bam'), file('../run_dorado_v5/barcode11/polished_assembly_11.fasta')],
        [ file('../run_dorado_v5/barcode12/calling_barcode12.bam'), file('../run_dorado_v5/barcode12/polished_assembly_12.fasta')]
    ]
*/

// load modules
include { bam2fastq; zipfastq; minimap2; index } from './modules/map_index_bam.nf'
include { modkit_pileup; modkit_pileup_bedgraphs; modkit_find_motifs; custom_bedgraphs; split_bed_by_bin} from './modules/modkit.nf'
include { compute_statistics } from './modules/statistics.nf'


// main workflow
workflow {
    
    // combine FASTA and BAM channels to generate a channel: tuple val(sample_id), path(bam_file), path(reference) 
    bam_ref_pairs = bam_input_ch.join(fasta_input_ch)
    
    fastq_files = bam2fastq(bam_ref_pairs)
    zipfastq(fastq_files)
    mapped_bams = minimap2(fastq_files)

    index_bam = index(mapped_bams)

    // Join the channels by sample ID
    bam_and_index = mapped_bams
        .join(index_bam)
    
    bed_file = modkit_pileup(bam_and_index)

    if (params.meta) {
        if (!params.bin_folder) {
            error "--bin_folder should be provided when using --meta"
        }

        bins = Channel.fromPath("${params.bin_folder}/*.fasta")

        motifs = split_bed_by_bin(bed_file, bins)
    } else {
        modkit_pileup_bedgraphs(bam_and_index)
        modkit_find_motifs(bed_file)
        custom_bedgraphs(bed_file)
        compute_statistics(bed_file)
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
