process modkit_pileup {
    label 'modkit'
    // execute the modkit pileup command
    
    input:
    tuple val(sample_id), path(mapped_bam), path(mapped_bam_bai), path(reference)

    output:
    tuple val(reference.baseName), val(sample_id), path(reference), path("modkit_pileup_output.bed")

    script:
    filter_threshold = params.automatic_threshold_modkit ? '' : "--filter-threshold ${params.filter_threshold_modkit}"
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} modkit_pileup_output.bed ${filter_threshold}
    """
    stub:
    """
    touch modkit_pileup_output.bed   
    """
}

process modkit_pileup_bedgraphs {
    label 'modkit'
    // execute the modkit pileup command and convert to per-modification bigWig tracks
    // (`pileup --bedgraph` was removed in modkit v0.6.0; bigWig is the recommended replacement,
    // see https://github.com/nanoporetech/modkit/blob/master/book/src/migrating_060.md)

    input:
    tuple val(sample_id), path(mapped_bam), path(mapped_bam_bai), path(reference)

    output:
    tuple val(reference.baseName), path("bigwigs")

    script:
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} modkit_pileup_bedgraph_output.bed --filter-threshold ${params.filter_threshold_modkit}
    mkdir -p bigwigs
    modkit bedmethyl tobigwig modkit_pileup_bedgraph_output.bed bigwigs/6mA.bw --mod-codes a --header ${mapped_bam} --negative-strand-values -t ${task.cpus}
    modkit bedmethyl tobigwig modkit_pileup_bedgraph_output.bed bigwigs/5mC.bw --mod-codes m --header ${mapped_bam} --negative-strand-values -t ${task.cpus}
    modkit bedmethyl tobigwig modkit_pileup_bedgraph_output.bed bigwigs/4mC.bw --mod-codes 21839 --header ${mapped_bam} --negative-strand-values -t ${task.cpus}
    """
    stub:
    """
    mkdir -p bigwigs
    touch bigwigs/6mA.bw bigwigs/5mC.bw bigwigs/4mC.bw
    """
}

process custom_bedgraphs {
    label 'biopython'
    // run an inhouse script that computes which bases are methylated: modified bases / total bases.
    // additionally, it saves the positions with a high methylation levels (>0.5) in tsv tables, one for each modification (6mA, 5mC and 4mC).

    input:
    tuple val(reference_name), val(sample_id), path(reference), path(bed_file)

    output:
    tuple val(reference_name), path("bedgraphs_customized"), path("modifications_tables")

    script:
    """
    mkdir -p bedgraphs_customized
    mkdir -p modifications_tables

    custom_bedgraphs.py ${bed_file} ${reference} . --percent_cutoff ${params.percent_cutoff_modification_table}
    """
    stub:
    """
    mkdir -p bedgraphs_customized
    mkdir -p modifications_tables
    """
}

process modkit_find_motifs {
    label 'modkit_low'
    // find motifs from the output of modkit pileup
    // `find-motifs` was renamed to `motif search` in modkit v0.6.0 (same flags)
    //errorStrategy 'ignore'

    input:
    tuple val(reference_name), val(sample_id), path(reference), path(bed_file)

    output:
    tuple val(reference_name), path("modkit_motifs.tsv")

    script:
    """
    modkit motif search -t ${task.cpus} --in-bedmethyl ${bed_file} --ref ${reference} -o modkit_motifs.tsv
    """
    stub:
    """
    touch modkit_motifs.tsv
    """ 
}


process publish_results_meta {
    label 'publish'
    publishDir "${params.outdir}/${sample_id}/bins/${reference_name}", mode: 'copy'

    input:
    tuple val(reference_name), val(sample_id), path(reference), 
          path(bed_file), path(pileup_bedgraphs), path(custom_bedgraphs), path(modifications_tables), path(statistics)

    output:
    tuple path(bed_file), path(pileup_bedgraphs),  path(custom_bedgraphs), path(modifications_tables), path(statistics)

    script:
    """
    """
}

process publish_results_motifs_meta {
    label 'publish'
    publishDir "${params.outdir}/${sample_id}/bins/${reference_name}", mode: 'copy'

    input:
    tuple val(reference_name), val(sample_id), path(reference), 
          path(bed_file), path(motifs)

    output:
    path(motifs)

    script:
    """
    """
}

process publish_results {
    label 'publish'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(reference_name), val(sample_id), path(reference), 
          path(bed_file), path(pileup_bedgraphs), path(custom_bedgraphs), path(modifications_tables), path(statistics)

    output:
    tuple path(bed_file), path(pileup_bedgraphs), path(custom_bedgraphs), path(modifications_tables), path(statistics)

    script:
    """
    """
}

process publish_results_motifs {
    label 'publish'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(reference_name), val(sample_id), path(reference), 
          path(bed_file), path(motifs)

    output:
    path(motifs)

    script:
    """
    """
}