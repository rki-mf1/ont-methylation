process bam2fastq {
    // converts the BAM file from Dorado to FASTQ format 
    input:
    path bam_file

    output:
    path "calling.fastq" 

    script:
    """
    samtools fastq ${bam_file} -T MM,ML > calling.fastq

    """
}

process minimap2 {
    // align with minimap2 

    cpus 8
     
    input:
    path fastq_file
    path reference

    output:
    path "methylation_mapped.bam" 

    script:
    """
    minimap2 -t $task.cpus --secondary=no -ax map-ont -y ${reference} ${fastq_file} | \
    samtools view -b | \
    samtools sort -@ $task.cpus -o methylation_mapped.bam

    """
}


process index {
    // index an aligned and sorted BAM
    input:
    path bam_file

    output:
    path "${bam_file}.bai"

    script:
    """
    samtools index ${bam_file}
    """
}
