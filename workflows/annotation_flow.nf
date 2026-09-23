include { methylation_density; motif_density; combine_peaks; circular_plot; promoter_analysis } from '../modules/annotation'

workflow ANNOTATION_FLOW {
    take:
        fasta_ch
        modkit_bed_ch
        gff3_ch

    main:
        if (!params.modkit_motifs) {
            println "\033[0;33mNote: --modkit_motifs not provided — motif-based density and promoter analysis will be skipped. Include --modkit_motifs (from modkit motif search) to enable them.\033[0m"
        }

        meth_input_ch = fasta_ch
            .combine(modkit_bed_ch)
            .combine(gff3_ch)

        methylation_density(meth_input_ch)

        if (params.modkit_motifs) {

            motif_input_ch = fasta_ch
                .combine(gff3_ch)
                .combine(Channel.fromPath(params.modkit_motifs, checkIfExists: true))

            motif_density(motif_input_ch)

            motif_peaks_ch = motif_density.out
                .map { sample_id, density_csvs, peaks_csvs ->
                    tuple(sample_id, peaks_csvs)
                }

        } else {
            motif_peaks_ch = methylation_density.out
                .map { sample_id, density_csvs, peaks_csvs ->
                    tuple(sample_id, [])
                }
        }

    
        combine_input_ch = methylation_density.out
            .map { sample_id, density_csvs, peaks_csvs ->
                tuple(sample_id, peaks_csvs)
            }
            .join(motif_peaks_ch, by: 0)


        combine_input_ch.view()
        combine_peaks(combine_input_ch)

        circular_meth_ch = methylation_density.out
            .flatMap { sample_id, density_csvs, peaks_csvs ->
                def d_list = density_csvs instanceof List ? density_csvs : [density_csvs]
                def p_list = peaks_csvs   instanceof List ? peaks_csvs   : [peaks_csvs]

                def d_map = d_list.collectEntries { f ->
                    def mod = (f.name =~ /density_(.+)\.csv/)[0][1]
                    [mod, f]
                }
                def p_map = p_list.collectEntries { f ->
                    def mod = (f.name =~ /peaks_(.+)\.csv/)[0][1]
                    [mod, f]
                }
                d_map.keySet().intersect(p_map.keySet()).collect { mod ->
                    tuple(sample_id, mod, p_map[mod], d_map[mod])
                }
            }

        if (params.modkit_motifs) {
            circular_motif_ch = motif_density.out
                .flatMap { sample_id, density_csvs, peaks_csvs ->
                    def d_list = density_csvs instanceof List ? density_csvs : [density_csvs]
                    def p_list = peaks_csvs   instanceof List ? peaks_csvs   : [peaks_csvs]

                    def d_map = d_list.collectEntries { f ->
                        def mod = (f.name =~ /density_(.+)\.csv/)[0][1]
                        [mod, f]
                    }
                    def p_map = p_list.collectEntries { f ->
                        def mod = (f.name =~ /peaks_(.+)\.csv/)[0][1]
                        [mod, f]
                    }
                    d_map.keySet().intersect(p_map.keySet()).collect { mod ->
                        tuple(sample_id, mod, p_map[mod], d_map[mod])
                    }
                }
            circular_ch = circular_meth_ch.mix(circular_motif_ch)
        } else {
            circular_ch = circular_meth_ch
        }
        circular_plot(circular_ch)

        // promoter_analysis: requires both modkit_motifs and promoter_analysis flag
        if (params.promoter_analysis && params.modkit_motifs) {
            promoter_input_ch = fasta_ch
                .combine(modkit_bed_ch)
                .combine(gff3_ch)
                .combine(Channel.fromPath(params.modkit_motifs, checkIfExists: true))
            promoter_analysis(promoter_input_ch)
        } else if (params.promoter_analysis && !params.modkit_motifs) {
            println "\033[0;33mNote: --promoter_analysis was requested but --modkit_motifs is missing — skipping promoter analysis.\033[0m"
        }


    emit:
        density  = methylation_density.out
        combined = combine_peaks.out
        circular = circular_plot.out
}