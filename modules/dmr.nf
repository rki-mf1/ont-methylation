process verify_same_reference {
    label 'minimap2'
    // when the user hands us already-compressed pileup beds (instead of BAMs), we skip
    // mapping entirely -- but modkit dmr pair silently produces garbage if the samples
    // weren't actually pileup'd against the same reference, so check contig names up front.

    input:
    tuple val(sample_id), path(bed_gz), path(bed_gz_tbi), path(reference)

    output:
    tuple val(sample_id), path(bed_gz), path(bed_gz_tbi)

    script:
    """
    samtools faidx ${reference}
    cut -f1 ${reference}.fai | sort -u > ref_contigs.txt
    tabix -l ${bed_gz} | sort -u > bed_contigs.txt
    if [ -s bed_contigs.txt ] && ! comm -23 bed_contigs.txt ref_contigs.txt | grep -q .; then
        echo "OK: all contigs in ${bed_gz} are present in ${reference}"
    else
        echo "❌ ${bed_gz} has contigs not found in ${reference} -- it was likely pileup'd against a different reference." >&2
        echo "Contigs only in ${bed_gz}:" >&2
        comm -23 bed_contigs.txt ref_contigs.txt >&2
        exit 1
    fi
    """
    stub:
    """
    echo stub
    """
}

process reclassify_coverage {
    label 'biopython'
    // "strict" coverage mode: fold below-threshold/ambiguous reads into the unmodified
    // count before modkit ever sees the pileup, so its own statistics get computed against
    // total coverage instead of modkit's default confident-reads-only (valid) coverage.
    // Empirically this is a strict subset of the valid-coverage hit list (~95% fewer sites,
    // but 100% of what survives was already in the valid-coverage set) -- see
    // TODO_dmr_reimplementation.md for the full writeup.

    input:
    tuple val(sample_id), path(bed_gz), path(bed_gz_tbi)

    output:
    tuple val(sample_id), path("${sample_id}_strict.bed")

    script:
    """
    reclassify_coverage.py ${bed_gz} ${sample_id}_strict.bed
    """
    stub:
    """
    touch ${sample_id}_strict.bed
    """
}

process compress_index_reclassified {
    label 'minimap2'

    input:
    tuple val(sample_id), path(bed)

    output:
    tuple val(sample_id), path("${bed}.gz"), path("${bed}.gz.tbi")

    script:
    """
    bgzip -c ${bed} > ${bed}.gz
    tabix -p bed ${bed}.gz
    """
    stub:
    """
    touch ${bed}.gz
    touch ${bed}.gz.tbi
    """
}

process dmr_pair {
    label 'modkit'
    // pairwise, single-site DMR comparison for one modification base between two samples

    input:
    tuple val(sample_a), path(bed_a, stageAs: "a.bed.gz"), path(tbi_a, stageAs: "a.bed.gz.tbi"),
          val(sample_b), path(bed_b, stageAs: "b.bed.gz"), path(tbi_b, stageAs: "b.bed.gz.tbi"),
          val(venn_label), val(base), val(coverage_mode), path(reference)

    output:
    tuple val(sample_a), val(sample_b), val(base), val(venn_label), val(coverage_mode),
          path("dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}.bed")

    script:
    """
    modkit dmr pair \
        -a ${bed_a} \
        -b ${bed_b} \
        -o dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}.bed \
        --ref ${reference} \
        --base ${base} \
        -t ${task.cpus} \
        --log-filepath dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}.log
    """
    stub:
    """
    touch dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}.bed
    """
}

process filter_dmr {
    label 'biopython'
    // the raw dmr_pair bed is 400MB+ per comparison on a real genome -- only publish the
    // small filtered table, the raw file just passes through for volcano_plot/annotate_dmr
    publishDir { "${params.outdir}/dmr/${sample_a}_vs_${sample_b}/${coverage_mode}" }, mode: 'copy', pattern: "*_filtered.tsv"

    input:
    tuple val(sample_a), val(sample_b), val(base), val(venn_label), val(coverage_mode), path(raw_bed)

    output:
    tuple val(sample_a), val(sample_b), val(base), val(venn_label), val(coverage_mode), path(raw_bed),
          path("dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_filtered.tsv")

    script:
    def score_arg = params.dmr_min_score ? "--min-score ${params.dmr_min_score}" : ""
    """
    filter_dmr.py ${raw_bed} dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_filtered.tsv \
        --min-coverage ${params.dmr_min_coverage} \
        --max-pvalue ${params.dmr_max_pvalue} \
        --min-effect ${params.dmr_min_effect} \
        ${score_arg}
    """
    stub:
    """
    touch dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_filtered.tsv
    """
}

process volcano_plot {
    label 'annotation'
    publishDir { "${params.outdir}/dmr/${sample_a}_vs_${sample_b}/${coverage_mode}" }, mode: 'copy'

    input:
    tuple val(sample_a), val(sample_b), val(base), val(venn_label), val(coverage_mode), path(raw_bed), path(filtered_tsv)

    output:
    path("dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_volcano.png")

    script:
    """
    dmr_volcano_plot.py \
        --raw ${raw_bed} \
        --filtered ${filtered_tsv} \
        --output dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_volcano.png \
        --title "${sample_a} vs ${sample_b} (${base}, ${coverage_mode})" \
        --min-coverage ${params.dmr_min_coverage} \
        --min-effect ${params.dmr_min_effect}
    """
    stub:
    """
    touch dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_volcano.png
    """
}

process annotate_dmr {
    label 'annotation'
    publishDir { "${params.outdir}/dmr/${sample_a}_vs_${sample_b}/${coverage_mode}" }, mode: 'copy'

    input:
    tuple val(sample_a), val(sample_b), val(base), val(venn_label), val(coverage_mode), path(raw_bed), path(filtered_tsv), path(reference), path(gff3)

    output:
    tuple val(sample_a), val(sample_b), val(base), val(venn_label), val(coverage_mode),
          path("dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_annotated.tsv"),
          path("dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_annotated_intergenic.tsv"),
          path("dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_meme.fasta")

    script:
    """
    annotate_dmr.py \
        --dmr ${filtered_tsv} \
        --gff3 ${gff3} \
        --fasta ${reference} \
        --output dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_annotated.tsv \
        --meme dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_meme.fasta \
        --window 50
    """
    stub:
    """
    touch dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_annotated.tsv
    touch dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_annotated_intergenic.tsv
    touch dmr_${sample_a}_${sample_b}_${base}_${coverage_mode}_meme.fasta
    """
}

process combine_dmr_results {
    label 'annotation'
    publishDir { "${params.outdir}/dmr/combined/${coverage_mode}" }, mode: 'copy'

    input:
    tuple val(base), val(coverage_mode), val(labels), path(filtered_tsvs)

    output:
    path("dmr_combined_${base}_${coverage_mode}.tsv")
    path("dmr_overlap_${base}_${coverage_mode}_summary.tsv")
    path("dmr_overlap_${base}_${coverage_mode}.png")

    script:
    """
    combine_dmr_results.py \
        --labels ${labels.join(' ')} \
        --tsvs ${filtered_tsvs} \
        --base ${base}_${coverage_mode} \
        --outdir .
    """
    stub:
    """
    touch dmr_combined_${base}_${coverage_mode}.tsv
    touch dmr_overlap_${base}_${coverage_mode}_summary.tsv
    touch dmr_overlap_${base}_${coverage_mode}.png
    """
}

process rank_dmr_genes {
    label 'biopython'
    publishDir { "${params.outdir}/dmr/combined/${coverage_mode}" }, mode: 'copy'
    // which genes carry the most differential methylation, across every comparison and
    // modification base -- ranked by number of distinct DMR sites they contain.

    input:
    tuple val(coverage_mode), path(annotated_tsvs)

    output:
    path("dmr_gene_ranking_${coverage_mode}.tsv")
    path("dmr_gene_ranking_${coverage_mode}_by_name.tsv")

    script:
    """
    rank_dmr_genes.py --annotated ${annotated_tsvs} --output dmr_gene_ranking_${coverage_mode}.tsv
    """
    stub:
    """
    touch dmr_gene_ranking_${coverage_mode}.tsv
    touch dmr_gene_ranking_${coverage_mode}_by_name.tsv
    """
}

process discover_motifs {
    label 'meme'
    publishDir { "${params.outdir}/dmr/${sample_a}_vs_${sample_b}" }, mode: 'copy'
    // de novo motif discovery (streme, from the MEME Suite) on the sequence context around
    // this comparison's DMR sites, one run per comparison -- kept separate rather than pooled
    // across comparisons, since different comparisons can be driven by different underlying
    // motifs/mechanisms and mixing them together risks washing out or confusing both signals.
    // (Mapping a *known* motif's genome-wide positions is a separate, main-flow concern --
    // see bin/map_motif.py -- this step is specifically about finding new candidate motifs
    // from the DMR sites themselves, which only exist in this flow.)
    //
    // --objfun cd (Central Distance) instead of streme's default (Differential Enrichment):
    // annotate_dmr.py's get_sequence_context() always puts the actual DMR base at the exact
    // center of every sequence, so we reward motifs that consistently occur there rather than
    // just "over-represented anywhere in the window vs. a shuffled background".

    input:
    tuple val(sample_a), val(sample_b), val(base), path(meme_fasta)

    output:
    tuple val(sample_a), val(sample_b), val(base), path("streme_${sample_a}_${sample_b}_${base}")

    script:
    """
    n_sites=\$(grep -c "^>" ${meme_fasta} || true)
    if [ "\$n_sites" -lt ${params.dmr_motif_min_sites} ]; then
        mkdir -p streme_${sample_a}_${sample_b}_${base}
        echo "Skipped motif discovery for ${sample_a} vs ${sample_b} (${base}): only \$n_sites DMR site sequence(s), need >= ${params.dmr_motif_min_sites} (--dmr_motif_min_sites)." > streme_${sample_a}_${sample_b}_${base}/SKIPPED.txt
    else
        streme --p ${meme_fasta} \
            --dna \
            --objfun cd \
            --minw ${params.dmr_motif_minw} \
            --maxw ${params.dmr_motif_maxw} \
            --oc streme_${sample_a}_${sample_b}_${base}
        summarize_streme_motifs.py \
            --xml streme_${sample_a}_${sample_b}_${base}/streme.xml \
            --output streme_${sample_a}_${sample_b}_${base}/motifs_summary.tsv
    fi
    """
    stub:
    """
    mkdir -p streme_${sample_a}_${sample_b}_${base}
    touch streme_${sample_a}_${sample_b}_${base}/streme.txt
    touch streme_${sample_a}_${sample_b}_${base}/motifs_summary.tsv
    """
}

