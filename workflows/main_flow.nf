include { bam2fastq; zipfastq; minimap2; split_bam_by_bin } from '../modules/map_index_bam.nf'
include { modkit_pileup; compress_index_modkit_bed; modkit_pileup_bigwigs; modkit_find_motifs; compute_methylation_tracks; methylation_tracks_to_bigwig; publish_results_meta; publish_results_motifs_meta; publish_results; publish_results_motifs} from '../modules/modkit.nf'
include { compute_statistics } from '../modules/statistics.nf'
include { capture_minimap2_samtools_version; capture_modkit_version; write_versions_summary } from '../modules/versions.nf'

workflow MAIN_FLOW {
    take:
        bam_input_ch
        fasta_input_ch
        bins_ch
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
            bam_bin_pairs = mapped_bams.combine(bins_ch)
            filtered_bams = split_bam_by_bin(bam_bin_pairs)
            
            bed_file = modkit_pileup(filtered_bams)
            bigwigs_modkit_ch = modkit_pileup_bigwigs(filtered_bams)
        } else {
            bed_file = modkit_pileup(mapped_bams)
            bigwigs_modkit_ch = modkit_pileup_bigwigs(mapped_bams)
        }

        // bgzip + tabix-index the pileup bed: needed by modkit dmr, and usable directly by the annotation flow too
        bed_gz_ch = compress_index_modkit_bed(bed_file)

        motifs_ch = modkit_find_motifs(bed_file)
        methylation_tracks = compute_methylation_tracks(bed_file)
        bigwigs_custom_ch = methylation_tracks_to_bigwig(methylation_tracks.tracks)
        modifications_tables_ch = methylation_tracks.tables
        statistics_ch = compute_statistics(bed_file)

        publish_input = bed_file.join(bed_gz_ch.map { reference_name, sample_id, bed_gz, bed_gz_tbi -> tuple(reference_name, bed_gz, bed_gz_tbi) })
                                .join(bigwigs_modkit_ch)
                                .join(bigwigs_custom_ch)
                                .join(modifications_tables_ch)
                                .join(statistics_ch)

        publish_motifs_input = bed_file.join(motifs_ch)

        if (params.meta) {
            publish_results_meta(publish_input)
            publish_results_motifs_meta(publish_motifs_input)
        } else {
            publish_results(publish_input)
            publish_results_motifs(publish_motifs_input)
        }

        write_versions_summary(capture_minimap2_samtools_version(), capture_modkit_version())

    emit:
        bed_file
        bed_gz_ch  // bgzip+tabix bed, for the annotation and dmr flows
        motifs_ch
        mapped_bams  // for DMR flow later

}