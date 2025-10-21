process modkit_pileup {
    label 'modkit'
    // execute the modkit pileup command
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(mapped_bam), path(reference), path(index_bam)

    output:
    tuple val(sample_id), path("${sample_id}/modkit_pileup_output.bed"), path(reference)

    script:
    filter_threshold = params.automatic_threshold_modkit ? '' : "--filter-threshold ${params.filter_threshold_modkit}"
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} ${sample_id}/modkit_pileup_output.bed ${filter_threshold}
    """
}

process modkit_pileup_bedgraphs {
    label 'modkit'
    // execute the modkit pileup command to obtain the bedgraphs
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(mapped_bam), path(reference), path(index_bam)

    output:
    path("${sample_id}/bedgraphs")

    script:
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} --bedgraph ${sample_id}/bedgraphs --filter-threshold ${params.filter_threshold_modkit} 
    """
}

process custom_bedgraphs {
    label 'biopython'
    // run an inhouse script that computes which bases are methylated: modified bases / total bases.
    // additionally, it saves the positions with a high methylation levels in tsv tables, one for each modification (6mA, 5mC and 4mC).
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(bed_file), path(reference)

    output:
    path("${sample_id}/bedgraphs_customized")
    path("${sample_id}/modifications_tables")

    script:
    """
    mkdir -p ${sample_id}/bedgraphs_customized
    mkdir -p ${sample_id}/modifications_tables

    custom_bedgraphs.py ${bed_file} ${reference} ${sample_id} --percent_cutoff ${params.percent_cutoff_modification_table}
    """
}

process modkit_find_motifs {
    label 'modkit_low'
    // find motifs from the output of modkit pileup
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(bed_file), path(reference)

    output:
    path "${sample_id}/modkit_motifs.tsv" 

    script:
    """
    mkdir -p ${sample_id}
    modkit find-motifs -t ${task.cpus} --in-bedmethyl ${bed_file} --ref ${reference} -o ${sample_id}/modkit_motifs.tsv
    """
}

process split_bed_by_bin {
    label 'modkit_low'
    publishDir "${params.outdir}/split_beds", mode: 'copy'

    input:
    path bed_file
    path bin_fa from bins

    output:
    tuple val(bin_base), path("${bin_base}.bed"), path(bin_fa) // use this to call modkit find motifs per bin

    script:
    bin_base = bin_fa.baseName() // not sure you need this

    """
    # here extract the contigs belonging to the current bin 
    awk '/^>/{gsub(/^>/,""); split(\$1,a," "); print a[1]}' ${bin_fa} > ${bin_base}.contigs.txt
    # here you want to filter the bed file according to to the contigs in the current bin
    awk 'NR==FNR {keep[\$1]; next} (\$1 in keep)' ${bin_base}.contigs.txt ${bed_file} > ${bin_base}.bed

    rm -f ${bin_base}.contigs.txt # not sure what is this for
    """
}
