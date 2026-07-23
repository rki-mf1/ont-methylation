include { bam2fastq; zipfastq; minimap2; split_bam_by_bin } from '../modules/map_index_bam.nf'
include { modkit_pileup; modkit_pileup_bedgraphs; modkit_find_motifs; custom_bedgraphs; publish_results_meta; publish_results_motifs_meta; publish_results; publish_results_motifs} from '../modules/modkit.nf'
include { compute_statistics } from '../modules/statistics.nf'

workflow MAIN_FLOW {
    take:
        bam_input_ch
        fasta_input_ch
    main:
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

    emit:
        bed_file
        motifs_ch
        mapped_bams  // for DMR flow late

}