// modkit dmr pair's --base C pools 5mC and 4mC together (its own a_counts/b_counts columns
// show calls split by code, e.g. "m:316,21839:0" at a single site -- a "C" DMR result can be
// driven by either or both). --single-code narrows a --base run down to one modification code
// in isolation, which is how we get real 5mC-vs-4mC separation -- it's a filter on top of
// --base, not a replacement for it ("Error! need to specify at least 1 modified base" if
// --base is left out, confirmed against modkit 0.6.3 directly).
MOD_CODE = ["6mA": "a", "5mC": "m", "4mC": "21839"]
PRIMARY_BASE = ["6mA": "A", "5mC": "C", "4mC": "C"]

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
    // pairwise, single-site DMR comparison for one modification code between two samples

    input:
    tuple val(sample_a), path(bed_a, stageAs: "a.bed.gz"), path(tbi_a, stageAs: "a.bed.gz.tbi"),
          val(sample_b), path(bed_b, stageAs: "b.bed.gz"), path(tbi_b, stageAs: "b.bed.gz.tbi"),
          val(venn_label), val(modification), val(coverage_mode), path(reference)

    output:
    tuple val(sample_a), val(sample_b), val(modification), val(venn_label), val(coverage_mode),
          path("dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}.bed")

    script:
    def mod_code = MOD_CODE[modification]
    def primary_base = PRIMARY_BASE[modification]
    """
    modkit dmr pair \
        -a ${bed_a} \
        -b ${bed_b} \
        -o dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}.bed \
        --ref ${reference} \
        --base ${primary_base} \
        --single-code ${mod_code} \
        -t ${task.cpus} \
        --log-filepath dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}.log
    """
    stub:
    """
    touch dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}.bed
    """
}

process filter_dmr {
    label 'biopython'
    // this IS the DMR sites table (modkit dmr pair's raw output is every modifiable
    // position genome-wide, mostly uninteresting -- this keeps only the ones that pass
    // coverage/p-value/effect-size thresholds). When --gff3 is given, annotate_dmr adds
    // gene names on top of this same table and that becomes the one published "sites"
    // table instead, so this publish only fires as the fallback when there's no gff3.
    // a table that's just a header and zero rows (e.g. a rare modification like 4mC under
    // the strict coverage mode routinely has nothing survive filtering) isn't published --
    // checked here in bash, where paths resolve against the task's own work dir, rather
    // than in a publishDir saveAs closure, which resolves file() against the pipeline's
    // launch directory instead and crashes trying to read a file that isn't there.
    publishDir { "${params.outdir}/dmr_analysis/dmr_sites/tables/${coverage_mode}" }, mode: 'copy', pattern: "${sample_a}_vs_${sample_b}_${modification}_${coverage_mode}.tsv", enabled: !params.gff3

    input:
    tuple val(sample_a), val(sample_b), val(modification), val(venn_label), val(coverage_mode), path(raw_bed)

    output:
    tuple val(sample_a), val(sample_b), val(modification), val(venn_label), val(coverage_mode), path(raw_bed),
          path("dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_filtered.tsv"), emit: dmr
    path("${sample_a}_vs_${sample_b}_${modification}_${coverage_mode}.tsv"), optional: true, emit: published

    script:
    def score_arg = params.dmr_min_score ? "--min-score ${params.dmr_min_score}" : ""
    """
    filter_dmr.py ${raw_bed} dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_filtered.tsv \
        --min-coverage ${params.dmr_min_coverage} \
        --max-pvalue ${params.dmr_max_pvalue} \
        --min-effect ${params.dmr_min_effect} \
        ${score_arg}
    n_lines=\$(wc -l < dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_filtered.tsv)
    if [ "\$n_lines" -gt 1 ]; then
        cp dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_filtered.tsv ${sample_a}_vs_${sample_b}_${modification}_${coverage_mode}.tsv
    fi
    """
    stub:
    """
    touch dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_filtered.tsv
    """
}

process volcano_plot {
    label 'annotation'
    publishDir { "${params.outdir}/dmr_analysis/dmr_sites/plots/${coverage_mode}" }, mode: 'copy'

    input:
    tuple val(sample_a), val(sample_b), val(modification), val(venn_label), val(coverage_mode), path(raw_bed), path(filtered_tsv)

    output:
    path("dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_volcano.png")

    script:
    """
    dmr_volcano_plot.py \
        --raw ${raw_bed} \
        --filtered ${filtered_tsv} \
        --output dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_volcano.png \
        --title "${sample_a} vs ${sample_b} (${modification}, ${coverage_mode})" \
        --min-coverage ${params.dmr_min_coverage} \
        --min-effect ${params.dmr_min_effect}
    """
    stub:
    """
    touch dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_volcano.png
    """
}

process annotate_dmr {
    label 'annotation'
    // the sites table: filtered DMR sites (see filter_dmr) with gene annotation added --
    // this IS the DMR results, not an intermediate -- one row per site (one per overlapping
    // gene for genic sites, which already includes intergenic sites tagged as such, so
    // there's no separate intergenic file to also look at).
    // see filter_dmr for why emptiness is checked in bash rather than a publishDir saveAs
    publishDir { "${params.outdir}/dmr_analysis/dmr_sites/tables/${coverage_mode}" }, mode: 'copy', pattern: "${sample_a}_vs_${sample_b}_${modification}_${coverage_mode}.tsv"

    input:
    tuple val(sample_a), val(sample_b), val(modification), val(venn_label), val(coverage_mode), path(raw_bed), path(filtered_tsv), path(reference), path(gff3)

    output:
    tuple val(sample_a), val(sample_b), val(modification), val(venn_label), val(coverage_mode),
          path("annotated_${sample_a}_${sample_b}_${modification}_${coverage_mode}.tsv"),
          path("dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_meme.fasta"), emit: annotated
    path("${sample_a}_vs_${sample_b}_${modification}_${coverage_mode}.tsv"), optional: true, emit: published

    script:
    """
    annotate_dmr.py \
        --dmr ${filtered_tsv} \
        --gff3 ${gff3} \
        --fasta ${reference} \
        --output annotated_${sample_a}_${sample_b}_${modification}_${coverage_mode}.tsv \
        --meme dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_meme.fasta \
        --window 25
    n_lines=\$(wc -l < annotated_${sample_a}_${sample_b}_${modification}_${coverage_mode}.tsv)
    if [ "\$n_lines" -gt 1 ]; then
        cp annotated_${sample_a}_${sample_b}_${modification}_${coverage_mode}.tsv ${sample_a}_vs_${sample_b}_${modification}_${coverage_mode}.tsv
    fi
    """
    stub:
    """
    touch annotated_${sample_a}_${sample_b}_${modification}_${coverage_mode}.tsv
    touch dmr_${sample_a}_${sample_b}_${modification}_${coverage_mode}_meme.fasta
    """
}

process combine_dmr_results {
    label 'annotation'
    // see filter_dmr for why emptiness is checked in bash rather than a publishDir saveAs.
    // nothing downstream in DMR_FLOW consumes these three files (they're pure publish
    // outputs), so it's safe to just not create them at all when there's nothing to show.
    publishDir { "${params.outdir}/dmr_analysis/dmr_sites/tables/${coverage_mode}/combined" }, mode: 'copy', pattern: "*.tsv"
    publishDir { "${params.outdir}/dmr_analysis/dmr_sites/plots/${coverage_mode}/combined" }, mode: 'copy', pattern: "*.png"

    input:
    tuple val(modification), val(coverage_mode), val(labels), path(filtered_tsvs)

    output:
    path("dmr_combined_${modification}_${coverage_mode}.tsv"), optional: true
    path("dmr_overlap_${modification}_${coverage_mode}_summary.tsv"), optional: true
    path("dmr_overlap_${modification}_${coverage_mode}.png"), optional: true

    script:
    """
    combine_dmr_results.py \
        --labels ${labels.join(' ')} \
        --tsvs ${filtered_tsvs} \
        --base ${modification}_${coverage_mode} \
        --outdir .
    n_lines=\$(wc -l < dmr_combined_${modification}_${coverage_mode}.tsv)
    if [ "\$n_lines" -le 1 ]; then
        rm -f dmr_combined_${modification}_${coverage_mode}.tsv \
              dmr_overlap_${modification}_${coverage_mode}_summary.tsv \
              dmr_overlap_${modification}_${coverage_mode}.png
    fi
    """
    stub:
    """
    touch dmr_combined_${modification}_${coverage_mode}.tsv
    touch dmr_overlap_${modification}_${coverage_mode}_summary.tsv
    touch dmr_overlap_${modification}_${coverage_mode}.png
    """
}

process rank_dmr_genes {
    label 'biopython'
    publishDir { "${params.outdir}/dmr_analysis/dmr_sites/tables/${coverage_mode}/combined" }, mode: 'copy'
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
    label 'biopython'
    // center-aware motif detection (bin/detect_motifs.py) on the sequence context around
    // this comparison's DMR sites, one run per comparison -- kept separate rather than
    // pooled across comparisons, since different comparisons can be driven by different
    // underlying motifs/mechanisms and mixing them together risks washing out or confusing
    // both signals. (Mapping a *known* motif's genome-wide positions is a separate,
    // main-flow concern -- see bin/map_motif.py -- this step is specifically about finding
    // new candidate motifs from the DMR sites themselves, which only exist in this flow.)
    //
    // Replaces streme (MEME Suite, dropped entirely -- see TODO_dmr_reimplementation.md).
    // Every training sequence here is already centered on the real modified base, but
    // streme has no notion of that and has to blindly infer both a motif's width and its
    // alignment from scratch -- which repeatedly produced distorted, truncated motifs (e.g.
    // CTGGCTGC instead of the real CCWGG, missing the base that makes it recognizable).
    // detect_motifs.py only ever searches windows anchored at the known center -- see its
    // own module docstring for the full algorithm. Pure stdlib, no MEME Suite dependency.
    publishDir { "${params.outdir}/dmr_analysis/motifs" }, mode: 'copy', pattern: "${sample_a}_vs_${sample_b}_${modification}_motifs.tsv"

    input:
    tuple val(sample_a), val(sample_b), val(modification), path(meme_fasta)

    output:
    path("${sample_a}_vs_${sample_b}_${modification}_motifs.tsv"), optional: true

    script:
    """
    n_sites=\$(grep -c "^>" ${meme_fasta} || true)
    if [ "\$n_sites" -lt ${params.dmr_motif_min_sites} ]; then
        echo "Skipped motif discovery for ${sample_a} vs ${sample_b} (${modification}): only \$n_sites DMR site sequence(s), need >= ${params.dmr_motif_min_sites} (--dmr_motif_min_sites)."
    else
        detect_motifs.py \
            --fasta ${meme_fasta} \
            --output ${sample_a}_vs_${sample_b}_${modification}_motifs.tsv \
            --minw ${params.dmr_motif_minw} \
            --maxw ${params.dmr_motif_maxw}
        n_lines=\$(wc -l < ${sample_a}_vs_${sample_b}_${modification}_motifs.tsv)
        if [ "\$n_lines" -le 1 ]; then
            rm -f ${sample_a}_vs_${sample_b}_${modification}_motifs.tsv
        fi
    fi
    """
    stub:
    """
    touch ${sample_a}_vs_${sample_b}_${modification}_motifs.tsv
    """
}

