process compute_statistics {
    label 'biopython'
    // run an inhouse script that computes which bases are methylated: modified bases / total bases 
    
    input:
    tuple val(reference_name), val(sample_id), path(reference), path(bed_file)

    output:
    tuple val(reference_name), path("methylation_statistics")

    script:
    """
    compute_statistics.py ${bed_file} ${reference} methylation_statistics
    """
    stub:
    """
    mkdir -p methylation_statistics    
    """
}