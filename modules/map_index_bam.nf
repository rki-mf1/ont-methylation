process bam2fastq {
    // converts the BAM file from Dorado to FASTQ format 
    label 'minimap2'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.fastq")

    script:
    """
    samtools fastq ${bam_file} -T MM,ML > ${sample_id}.fastq
    """
    stub:
    """
    touch ${sample_id}.fastq    
    """
}

process zipfastq {
    // zip fastq file
    label 'minimap2'
    publishDir  "${params.outdir}/${sample_id}", mode:'copy'
        
    input:
    tuple val(sample_id), path(fastq_file)

    output:
    path "${sample_id}.fastq.gz"

    script:
    """
    gzip -c ${fastq_file} > ${sample_id}.fastq.gz
    """
    stub:
    """
    touch ${sample_id}.fastq.gz
    """
}

process minimap2 {
    label 'minimap2'
    // align with minimap2 
    publishDir  "${params.outdir}/${sample_id}", mode:'copy'

    input:
    tuple val(sample_id), path(fastq_file), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_methylation_mapped.bam"), path("${sample_id}_methylation_mapped.bam.bai"), path(reference)

    script:
    """
    minimap2 -t ${task.cpus} --secondary=no -ax map-ont -y ${reference} ${fastq_file} | \
    samtools view -b | \
    samtools sort -@ ${task.cpus} -o ${sample_id}_methylation_mapped.bam
    samtools index ${sample_id}_methylation_mapped.bam
    """
    stub:
    """
    touch ${sample_id}_methylation_mapped.bam
    touch ${sample_id}_methylation_mapped.bam.bai
    """
}


process split_bam_by_bin {
    label 'minimap2'
    publishDir "${params.outdir}/${sample_id}/bins/${bin_fa.baseName}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(reference), path(bin_fa)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path(bin_fa)

    script:
    """
    # extract contigs from bin fasta
    awk '/^>/{gsub(/^>/,""); split(\$1,a," "); print a[1]}' ${bin_fa} > ${bin_fa.baseName}.contigs.txt

    # filter bam by contigs
    samtools view -b ${bam_file} \$(cat ${bin_fa.baseName}.contigs.txt) -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
    stub:
    """
    touch ${sample_id}.bam
    touch ${sample_id}.bam.bai
    """
}
