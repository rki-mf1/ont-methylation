process capture_minimap2_samtools_version {
    label 'minimap2'
    // capture the actual minimap2/samtools versions from the container or conda env used at runtime

    output:
    path("minimap2_samtools_versions.txt")

    script:
    """
    echo "minimap2 \$(minimap2 --version)" > minimap2_samtools_versions.txt
    samtools --version | head -1 >> minimap2_samtools_versions.txt
    """
    stub:
    """
    touch minimap2_samtools_versions.txt
    """
}

process capture_modkit_version {
    label 'modkit'
    // capture the actual modkit version from the container or conda env used at runtime

    output:
    path("modkit_version.txt")

    script:
    """
    modkit --version > modkit_version.txt
    """
    stub:
    """
    touch modkit_version.txt
    """
}

process capture_basecalling_info {
    label 'minimap2'

    input:
    tuple val(sample_id), path(bam)

    output:
    path("${sample_id}_basecalling.tsv")

    script:
    """
    samtools view -H ${bam} > header.txt
    dorado=\$(awk -F'\t' '\$1=="@PG" && \$2=="ID:basecaller" {for(i=1;i<=NF;i++) if(\$i ~ /^VN:/) {sub(/^VN:/,"",\$i); sub(/[+].*/,"",\$i); print \$i; exit}}' header.txt)
    model=\$(grep '^@RG' header.txt | tr '\t ' '\n\n' | awk -F= '\$1=="basecall_model" {print \$2; exit}')
    mods=\$(grep '^@RG' header.txt | tr '\t ' '\n\n' | awk -F= '\$1=="modbase_models" {print \$2; exit}')
    printf "%s\t%s\t%s\t%s\n" "${sample_id}" "\${dorado:-unknown}" "\${model:-unknown}" "\${mods:-unknown}" > ${sample_id}_basecalling.tsv
    """
    stub:
    """
    printf "${sample_id}\tunknown\tunknown\tunknown\n" > ${sample_id}_basecalling.tsv
    """
}

process write_versions_summary {
    label 'publish'
    publishDir "${params.outdir}", mode: 'copy'
    // merge the captured tool versions into one human-readable Markdown summary for the run

    input:
    path(minimap2_samtools_versions)
    path(modkit_version)
    path(basecalling_tsvs)

    output:
    path("pipeline_versions.md")

    script:
    engine = workflow.containerEngine ?: (workflow.profile.tokenize(',').find { it in ['conda', 'mamba'] } ?: 'none (host tools)')
    """
    {
        echo "# Pipeline tool versions"
        echo ""
        echo "- Generated: \$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
        echo "- Nextflow version: ${nextflow.version}"
        echo "- Pipeline revision: ${workflow.revision ?: 'local (no -r used)'}"
        echo "- Execution engine: ${engine}"
        echo "- Execution profile: ${workflow.profile}"
        echo ""
        echo "| Tool | Version |"
        echo "|------|---------|"
        while read -r tool ver; do
            echo "| \$tool | \$ver |"
        done < ${minimap2_samtools_versions}
        read -r tool ver < ${modkit_version}
        echo "| \$tool | \$ver |"
        echo ""
        echo "## Basecalling (from BAM headers)"
        echo ""
        cat ${basecalling_tsvs} > basecalling_all.tsv
        n=\$(wc -l < basecalling_all.tsv)
        if [ \$(cut -f2- basecalling_all.tsv | sort -u | wc -l) -eq 1 ]; then
            IFS=\$'\t' read -r s dorado model mods < basecalling_all.tsv
            echo "| Setting | Value |"
            echo "|---------|-------|"
            echo "| dorado | \$dorado |"
            echo "| Basecall model | \$model |"
            echo "| Modification models | \$mods |"
            echo "| Samples | \$n (all identical) |"
        else
            echo "**Multiple models** - samples were not all basecalled with the same dorado version/models:"
            echo ""
            echo "| Sample | dorado | Basecall model | Modification models |"
            echo "|--------|--------|----------------|---------------------|"
            sort basecalling_all.tsv | awk -F'\t' '{print "| "\$1" | "\$2" | "\$3" | "\$4" |"}'
        fi
    } > pipeline_versions.md
    """
    stub:
    """
    touch pipeline_versions.md
    """
}
