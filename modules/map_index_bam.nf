process bam2fastq {
    // converts the BAM file from Dorado to FASTQ format 
    label 'minimap2'

    input:
    tuple val(sample_id), path(bam_file), path(reference) 

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}.fastq"), path(reference) 
    script:
    """
    mkdir -p ${sample_id}
    samtools fastq ${bam_file} -T MM,ML > ${sample_id}/${sample_id}.fastq
    """
}

process zipfastq {
    // zip fastq file
    label 'minimap2'
    publishDir  params.outdir, mode:'copy'
        
    input:
    tuple val(sample_id), path(fastq_file), path(reference)

    output:
    path "${sample_id}/${sample_id}.fastq.gz"

    script:
    """
    mkdir -p ${sample_id}
    gzip -c ${fastq_file} > ${sample_id}/${sample_id}.fastq.gz
    """
}

process minimap2 {
    label 'minimap2'
    // align with minimap2 
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(fastq_file), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}/methylation_mapped.bam"), path(reference)

    script:
    """
    mkdir -p ${sample_id}
    minimap2 -t ${task.cpus} --secondary=no -ax map-ont -y ${reference} ${fastq_file} | \
    samtools view -b | \
    samtools sort -@ ${task.cpus} -o ${sample_id}/methylation_mapped.bam
    """
}

process index {
    label 'minimap2'
    // index an aligned and sorted BAM
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(mapped_bam), path(reference)

    output:
    tuple val(sample_id), path ("${sample_id}/${mapped_bam}.bai")

    script:
    """
    mkdir -p ${sample_id}
    mv ${mapped_bam} ${sample_id}
    samtools index ${sample_id}/${mapped_bam} 
    
    """
}
