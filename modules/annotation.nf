process methylation_density {
    label 'annotation'
    publishDir { "${params.outdir}/${sample_id}/annotation/density" }, mode: 'copy', pattern: "density_*.csv"
    publishDir { "${params.outdir}/${sample_id}/annotation/tables" },   mode: 'copy', pattern: "peaks_*.csv"


    input:
    tuple val(sample_id), path(fasta), path(modkit_bed), path(gff3)

    output:
    tuple val(sample_id), path("density_*.csv"),
                          path("peaks_*.csv")

    script:
    """
    methylation_density.py \
        --modkit_bed  ${modkit_bed}  \
        --fasta       ${fasta}       \
        --gff3        ${gff3}        \
        --sample      ${sample_id}   \
        --outdir      . \
        --modifications              ${params.modifications_annotation} \
        --percent_modified_threshold ${params.percent_modified_threshold_annotation} \
        --top_n                      ${params.top_n} \
        --window_size                ${params.window_size} \
        --step_size                  ${params.step_size} \
        --smoothing_window           ${params.smoothing_window} \
        --enrichment_threshold       ${params.enrichment_threshold} \
        --min_sites                  ${params.min_sites_annotation}
    """
    stub:
    """
    touch density_6mA.csv
    touch peaks_6mA.csv
    """
}

process motif_density {
    label 'annotation'
    publishDir { "${params.outdir}/${sample_id}/annotation/density" }, mode: 'copy', pattern: "density_*.csv"
    publishDir { "${params.outdir}/${sample_id}/annotation/tables" },   mode: 'copy', pattern: "peaks_*.csv"

    input:
    tuple val(sample_id), path(fasta), path(gff3), path(motif_file)

    output:
    tuple val(sample_id), path("density_*.csv"),
                          path("peaks_*.csv")

    script:
    """
    motif_density.py \
        --fasta       ${fasta}      \
        --gff3        ${gff3}       \
        --motif_file  ${motif_file} \
        --sample      ${sample_id}  \
        --outdir      . \
        --top_n                  ${params.top_n} \
        --window_size            ${params.window_size} \
        --step_size              ${params.step_size} \
        --smoothing_window       ${params.smoothing_window} \
        --enrichment_threshold   ${params.enrichment_threshold} \
        --min_absolute_count     2
    """
    stub:
    """
    touch density_CGGYCG.csv
    touch peaks_CGGYCG.csv
    """
}

process combine_peaks {
    label 'annotation'
    publishDir { "${params.outdir}/${sample_id}/annotation/tables" }, mode: 'copy'

    input:
    tuple val(sample_id), path(meth_peaks_csvs), path(motif_peaks_csvs)

    output:
    tuple val(sample_id), path("all_peaks_combined.csv")

    script:
    def motif_arg = params.modkit_motifs ? "--motif_peaks ${motif_peaks_csvs}" : ""
    """
    combine_peaks.py \
        --meth_peaks  ${meth_peaks_csvs} \
        ${motif_arg} \
        --outdir      .
    """
    stub:
    """
    touch all_peaks_combined.csv
    """
}

process circular_plot {
    label 'annotation'
    publishDir { "${params.outdir}/${sample_id}/annotation/plots" }, mode: 'copy'

    input:
    tuple val(sample_id), val(modification), path(peaks_genic_csv), path(density_csv)

    output:
    tuple val(sample_id), val(modification),
          path("circular_density_genes_${sample_id}_${modification}.png")

    script:
    """
    circular_plot.py \
        --peaks  ${peaks_genic_csv}  \
        --density      ${density_csv}      \
        --modification ${modification}     \
        --species      ${sample_id}        \
        --top_n        ${params.top_n ?: 30} \
        --output       circular_density_genes_${sample_id}_${modification}.png
    """
    stub:
    """
    touch circular_density_genes_${sample_id}_${modification}.png
    """
}

process promoter_analysis {
    label 'annotation'
    publishDir { "${params.outdir}/${sample_id}/annotation/promoters" }, mode: 'copy'

    input:
    tuple val(sample_id), path(fasta), path(modkit_bed), path(gff3), path(motif_file)

    output:
    tuple val(sample_id), path("promoter_*.csv")

    script:
    """
    promoter_analysis.py \
        --fasta                      ${fasta}      \
        --modkit_bed                 ${modkit_bed} \
        --gff3                       ${gff3}       \
        --motif_file                 ${motif_file} \
        --sample                     ${sample_id}  \
        --outdir                     . \
        --promoter_window            ${params.promoter_window} \
        --percent_modified_promoter ${params.percent_modified_threshold_promoter}
    """
    stub:
    """
    touch promoter_CGGYCG.csv
    """
}