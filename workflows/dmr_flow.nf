include { bam2fastq; zipfastq; minimap2 } from '../modules/map_index_bam.nf'
include { modkit_pileup; compress_index_modkit_bed } from '../modules/modkit.nf'
include { verify_same_reference; reclassify_coverage; compress_index_reclassified;
          dmr_pair; filter_dmr; volcano_plot; annotate_dmr; combine_dmr_results;
          discover_motifs; rank_dmr_genes } from '../modules/dmr.nf'

// builds the all-vs-one pairing (every non-reference sample vs the reference) for one
// coverage mode's set of samples. "venn_label" is just the non-reference sample's own id --
// under all-vs-one every sample maps to exactly one comparison, so the overlap/Venn step can
// unambiguously name each circle after the sample itself ("does X differ from the reference").
def buildPairs(samples) {
    samples = samples.sort { it[0] }
    def n = samples.size()
    if (n < 2) {
        error "❌ --flow dmr needs at least 2 samples, found ${n}"
    }
    def ids = samples.collect { it[0] }
    def dupes = ids.findAll { ids.count(it) > 1 }.unique()
    if (dupes) {
        error "❌ Duplicate sample id(s) among --flow dmr inputs: ${dupes.join(', ')}. " +
              "This happens if multiple input files share the same basename (e.g. globbing " +
              "'results/*/modkit_pileup_output.bed.gz', which is identical per sample) -- use --list with an explicit sample_id,path CSV instead."
    }
    def ref_id = params.dmr_reference_sample ?: ids[0]
    def ref_sample = samples.find { it[0] == ref_id }
    if (!ref_sample) {
        error "❌ --dmr_reference_sample '${ref_id}' not found among input samples (${ids.join(', ')})"
    }
    def pairs = []
    samples.findAll { it[0] != ref_id }.each { s -> pairs << [ref_sample, s, s[0].toString()] }
    log.info "\033[0;33mDMR: ${n} samples, comparing all against reference '${ref_id}' -> ${pairs.size()} comparison(s)\033[0m"
    pairs
}

workflow DMR_FLOW {
    take:
        sample_input_ch   // tuple(sample_id, path) -- a BAM or an already-compressed pileup bed.gz, depending on input_is_bam
        fasta_ch          // value channel: the single reference every sample must share
        input_is_bam      // boolean

    main:
        if (input_is_bam) {
            fastq_files = bam2fastq(sample_input_ch)
            zipfastq(fastq_files)
            fastq_ref_pairs = fastq_files.combine(fasta_ch)
            mapped_bams = minimap2(fastq_ref_pairs)
            pileup_ch = modkit_pileup(mapped_bams)
            bed_gz_ch = compress_index_modkit_bed(pileup_ch)
                .map { reference_name, sample_id, bed_gz, bed_gz_tbi -> tuple(sample_id, bed_gz, bed_gz_tbi) }
        } else {
            bed_with_tbi_ch = sample_input_ch.map { sample_id, bed_gz ->
                def tbi = file("${bed_gz}.tbi")
                if (!tbi.exists()) {
                    error "❌ Missing tabix index for ${bed_gz} (expected ${bed_gz}.tbi next to it). bgzip + tabix your pileup beds before using --flow dmr with --modkit_bed."
                }
                tuple(sample_id, bed_gz, tbi)
            }
            // already-compressed input skips mapping entirely, so we can't be sure every
            // sample was actually pileup'd against the same reference -- check contig names.
            bed_gz_ch = verify_same_reference(bed_with_tbi_ch.combine(fasta_ch))
        }

        // "valid" coverage mode = modkit's own native behavior (confidently-classified
        // reads only). "strict" = reclassify ambiguous reads as unmodified first, so modkit's
        // statistics get computed against total coverage instead -- a stricter subset of the
        // valid-mode hits, see TODO_dmr_reimplementation.md. Both run side by side; only
        // "valid" feeds de novo motif discovery (see below).
        bed_gz_strict_ch = compress_index_reclassified(reclassify_coverage(bed_gz_ch))

        bases_ch = Channel.fromList(params.dmr_bases.split(",").collect { it.trim() })

        // DSL2 forbids calling the same process twice in one workflow scope, so rather than
        // running dmr_pair/filter_dmr/... once per coverage mode, tag every sample with its
        // mode up front and run each process exactly once over the merged stream -- the mode
        // just rides along as a value in the tuple and fans back out into separate
        // publishDir/filenames naturally (see modules/dmr.nf).
        bed_gz_tagged_ch = bed_gz_ch.map { sample_id, bed_gz, tbi -> tuple("valid", sample_id, bed_gz, tbi) }
            .mix(bed_gz_strict_ch.map { sample_id, bed_gz, tbi -> tuple("strict", sample_id, bed_gz, tbi) })

        pairs_ch = bed_gz_tagged_ch
            .map { mode, sample_id, bed_gz, tbi -> tuple(mode, [sample_id, bed_gz, tbi]) }
            .groupTuple(by: 0)
            .flatMap { mode, samples -> buildPairs(samples).collect { pair -> tuple(mode, pair[0], pair[1], pair[2]) } }
            .map { mode, a, b, venn_label -> tuple(a[0], a[1], a[2], b[0], b[1], b[2], venn_label, mode) }

        dmr_raw_ch = dmr_pair(
            pairs_ch.combine(bases_ch)
                .map { sample_a, bed_a, tbi_a, sample_b, bed_b, tbi_b, venn_label, mode, base ->
                    tuple(sample_a, bed_a, tbi_a, sample_b, bed_b, tbi_b, venn_label, base, mode) }
                .combine(fasta_ch)
        )
        dmr_filtered_ch = filter_dmr(dmr_raw_ch)

        volcano_plot(dmr_filtered_ch)

        if (!params.gff3) {
            println "\033[0;33mNote: --gff3 not provided -- DMR sites will not be annotated with gene names, de novo motif discovery and gene ranking (which both need that step's output) will be skipped.\033[0m"
        } else {
            gff3_ch = Channel.value(file(params.gff3, checkIfExists: true))
            annotated_ch = annotate_dmr(dmr_filtered_ch.combine(fasta_ch).combine(gff3_ch))

            // de novo motif discovery only makes sense on the "valid" (more sensitive) site
            // set -- it's an exploratory step looking for candidate motifs, not a confirmatory
            // one, and running it twice on largely-overlapping/subset data wastes a genuinely
            // slow step for no extra information.
            motif_input_ch = annotated_ch
                .filter { sample_a, sample_b, base, venn_label, coverage_mode, annotated_tsv, intergenic_tsv, meme_fasta -> coverage_mode == "valid" }
                .map { sample_a, sample_b, base, venn_label, coverage_mode, annotated_tsv, intergenic_tsv, meme_fasta -> tuple(sample_a, sample_b, base, meme_fasta) }

            discover_motifs(motif_input_ch)

            // gene ranking runs per coverage mode, so you can compare the "strict" ranking
            // against the full "valid" one the same way we did by hand earlier.
            gene_rank_input_ch = annotated_ch
                .map { sample_a, sample_b, base, venn_label, coverage_mode, annotated_tsv, intergenic_tsv, meme_fasta -> tuple(coverage_mode, annotated_tsv) }
                .groupTuple(by: 0)

            rank_dmr_genes(gene_rank_input_ch)
        }

        // combine every pairwise comparison's filtered sites (per modification base and
        // coverage mode) into one long table plus a reproducibility/overlap check
        combine_input_ch = dmr_filtered_ch
            .map { sample_a, sample_b, base, venn_label, coverage_mode, raw_bed, filtered_tsv -> tuple(base, coverage_mode, venn_label, filtered_tsv) }
            .groupTuple(by: [0, 1])

        combine_dmr_results(combine_input_ch)
}
