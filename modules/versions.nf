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

process write_versions_summary {
    label 'publish'
    publishDir "${params.outdir}", mode: 'copy'
    // merge the captured tool versions into one human-readable Markdown summary for the run

    input:
    path(minimap2_samtools_versions)
    path(modkit_version)

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
    } > pipeline_versions.md
    """
    stub:
    """
    touch pipeline_versions.md
    """
}
