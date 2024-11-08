# Nextflow pipeline for extracting bacterial methylation sites and motifs from Nanopore data

![](https://img.shields.io/badge/nextflow->=22.01.0-brightgreen)
![](https://img.shields.io/badge/can_use-docker-blue.svg)
![](https://img.shields.io/badge/can_use-singularity-orange.svg)
![](https://img.shields.io/badge/can_use-conda/mamba-yellow.svg)
![](https://img.shields.io/badge/licence-GLP3-lightgrey.svg)

This pipeline processes basecalled reads generated by [**Dorado**](https://github.com/nanoporetech/dorado), with DNA modification basecalling required. It extracts methylated positions along a reference genome and identifies specific motifs using [**Modkit**](https://github.com/nanoporetech/modkit).

Schematic overview of the pipeline:

![Alt text](images/workflow.png)

## Input

- **Basecalled reads**:  Outputs from Dorado in BAM file format, which include basecalled DNA modifications. See [Basecalling with Dorado](#basecalling-with-dorado) for details.
- **Reference file**: A reference genome or sequence against which the reads will be aligned. The assembly obtained from the BAM files is recommended (e.g., by _de novo_ assembly with [flye](https://github.com/mikolmogorov/Flye)). 

## Output

- **Motif analysis with Modkit**: Identified motifs based on DNA modifications using Modkit.
- **Statistical analysis**: Summary statistics of the methylation status for each contig in the reference genome.
- **List of methylated positions**: A list of methylated positions relative to the reference genome, including only those with a confidence level above a custom threshold.
- **BedGraph files for IGV visualization** BedGraph files produced by Modkit and custom refined BedGraph files are provided, displaying methylation confidence for each base. These files can be loaded into [IGV](https://igv.org/) for visualization as annotation tracks.


## Requirements

### **Workflow management and dependencies** 

You need to install [Nextflow](https://www.nextflow.io/docs/latest/install.html) to run the pipeline. We recommend installing Nextflow using `curl`:

```bash
curl -s https://get.nextflow.io | bash
```

However, if this does not work for you, you can also install Nextflow via [conda](https://docs.anaconda.com/miniconda/) or [mamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

To avoid installing all necessary tool dependencies manually, we recommend using [Docker](https://docs.docker.com/get-started/get-docker/) or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html). Nextflow will then handle all your dependencies using the Docker or Singularity profile (see below). 

**Not recommended!** Alternatively, the necessary dependencies can be also installed via [conda](https://docs.anaconda.com/miniconda/) or [mamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) outside of the Nextflow workflow:

```bash
conda create -n ONT_methylation -c bioconda samtools minimap2 biopython python=3.8 pandas
conda activate ONT_methylation
```

Modkit can be installed from the [github repository](https://github.com/nanoporetech/modkit) or via conda as well

```bash
conda install -c bioconda ont-modkit
```


### Preparing input: Basecalling with Dorado 

Your Nanopore data needs to be basecalled in a methyltion-aware way to be usable as input for the pipeline. [Dorado](https://github.com/nanoporetech/dorado) is the official open-source basecaller for Oxford Nanopore reads. To basecall pod5 files, it is strongly recommended to use the v5 models in super high accuracy mode (SUP), as these include detection for 6mA, 5mC, and 4mC modifications.

The following command will download the most recent models with super accuracy mode, using the methylation models for 6mA, 4mC and 5mC as well, and run Dorado. 

```bash
dorado basecaller sup,6mA,4mC_5mC <pod5_folder_path> >  results.bam
```
   
Alternatively, you can download specific models individually. In this case, you would need to download the basecalling model and the methylation models separately, then modify the Dorado command accordingly:
   
```bash
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1
```

Once the models are downloaded, you can run Dorado using the command below, specifying the basecalling model and adding the methylation models for modified bases:

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 <pod5_folder_path> --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1,dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1  > results.bam
```

#### Using Dorado duplex

Another option is to basecall with Dorado in duplex mode, which has recently been updated to support modification detection.

```bash
dorado duplex sup,6mA,4mC_5mC <pod5_folder_path> > results.bam
```

The `dx` tag in the BAM record allows you to distinguish between simplex and duplex reads. You can use this tag to filter out redundant reads if necessary.

However, it is not yet clear if duplex basecalling offers significant benefits for downstream analysis, as modification calls are still made on only one DNA strand.

## How to install and update the pipeline

To install the pipeline, simply use Nextflow:

```bash
nextflow pull valegale/ONT_methylation
# check available release versions and branches:
nextflow info valegale/ONT_methylation
# show the help message for a certain pipeline release version:
nextflow run valegale/ONT_methylation -r 0.0.1 --help
# update the pipeline simply via pulling the code again:
nextflow pull valegale/ONT_methylation
```

**Check to use the latest pipeline release version**. To have reproducible results, use the same version. 

You can also `git clone` this repository and run the pipeline via `nextflow run main.nf` - but we do not recommend this. 

## How to Run

To run the pipeline with a single reference and BAM file, use the following command (**adjust the `-r` version as necessary**):

```bash
nextflow run valegale/ONT_methylation -r 0.0.1 --fasta sample_test.fasta --bam sample_test.bam
```
   
- **fasta** specifies the reference genome file in FASTA format.
- **bam** specifies the basecalled BAM file output from Dorado (see above).

If you have multiple references and BAM files, you can provide them using wildcard patterns (`*`). For example:

```bash
nextflow run valegale/ONT_methylation -r 0.0.1 --fasta '*.fasta' --bam '*.bam'
```

**Don't forget the single ticks `'...'`!**

In this case, ensure that the base names of the FASTA and BAM files match exactly. This matching allows the pipeline to correctly associate each reference with its corresponding BAM file.

For complex cases with numerous references and BAM files, you can also use the `--list` parameter to specify CSV files as input instead of the direct paths to the files. This approach offers flexibility by allowing each reference and BAM file pair to be listed explicitly.

For example:

```bash
nextflow run valegale/ONT_methylation -r 0.0.1 --list --fasta references.csv --bam mappings.csv
```

The `references.csv` file should contain the sample name and path to each reference FASTA file:

```csv
sample1,/path/to/reference1.fasta
sample2,/path/to/reference2.fasta
sample3,/path/to/another/reference3.fasta
```

The `mappings.csv` file should list the sample name (matching the `references.csv`) and path to each BAM file:

```csv
sample1,/path/to/mapping1.bam
sample2,/path/to/mapping2.bam
sample3,/path/to/another/mapping3.bam
```


### Running with Containers

To run the pipeline using Docker, use the following command:

```bash
nextflow run valegale/ONT_methylation -r 0.0.1 --fasta sample_test.fasta --bam sample_test.bam -profile docker
```

To run the pipeline with SLURM and Singularity, use this command:
```bash
nextflow run nextflow pull valegale/ONT_methylation --fasta sample_test.fasta --bam sample_test.bam -profile slurm,singularity
```
   
## Results interpretation

### Motifs from Modkit

The results include a table generated by Modkit, named `modkit_motifs_tsv`, which summarizes the detected motifs. It is recommended to compare the detected motifs and their frequencies to determine if the same motif is modified for multiple modifications (also check the reverse complement of the motifs). Note that the code "21839" corresponds to the 4mC modification.

**Important:** Modkit applies a [confidence threshold](https://github.com/nanoporetech/modkit/blob/master/book/src/filtering.md) to filter modified positions. Only bases with a methylation confidence value above this threshold will be considered methylated. You can adjust this value using the `--filter_threshold_modkit` parameter (see also the `--help`). The default is set to 0.75, as values lower than this tend to produce a high number of false positives, while higher values may be overly stringent.

The intermediate output file, `modkit_pileup_output.bed`, is a tab-separated file that summarizes all reference positions along with their aggregated results from individual reads. This file serves as input for the next steps in the analysis.

### Methylation Level Statistics

Methylation statistics are stored in the `methylation_statistics` folder. Separate tables are created for each modification, with the percentage of methylation calculated by dividing the number of methylated bases (those exceeding the Modkit threshold) by the total number of relevant bases (A for 6mA, C for 4mC and 5mC).

### List of Methylated Positions and BedGraph files for Visualization

To extract the list of likely methylated positions, the output from Modkit (`modkit_pileup_output.bed`) must undergo filtering and refining (see [How Percent modified is computed](#How-percent-modified-is-computed) for more information). The tables containing these lists are located in the `modification_tables` folder.

Only positions with coverage greater than 10 and a methylation confidence above a specified threshold are included (default = 0.5). This threshold can be adjusted using the parameter `--percent_cutoff_modification_table`.

The BedGraph files generated directly from Modkit can be found in the `bedgraphs` folder, while the preprocessed BedGraph files, which reflect the same values as in the `modification_tables`, are located in the `bedgraphs_customized` folder.

## How "Percent Modified" is Computed

Modkit pileup aggregates the information from all reads for all positions of the reference genome and computes a value called "fraction modified," which represents the percentage of reads with a methylated base for each position. This calculation considers only the reads where the base has passed the confidence threshold; thus, if a position has a high number of reads with bases that didn't meet this threshold, those reads are excluded from the count.

In Modkit (check [Description of bedMethyl output](https://github.com/nanoporetech/modkit#Description-of-bedMethyl-output)), "fraction modified" is defined as N<sub>mod</sub> / N<sub>valid_cov</sub>. In our analysis, we introduce a new metric called "percent modified," computed as: N<sub>mod</sub> / (N<sub>valid_cov</sub> + N<sub>fail</sub> + N<sub>diff</sub>). This approach helps to prevent positions with only a few valid reads from being incorrectly classified as modified. See [Description of bedMethyl output](https://github.com/nanoporetech/modkit#Description-of-bedMethyl-output) for details regarding the before mentioned variables and their definition.


## Running MicrobeMod

[MicrobeMod](https://github.com/cultivarium/MicrobeMod) is a popular tool for motif annotation and identification. We recommend comparing the results from Modkit (generated using our pipeline) with those obtained from MicrobeMod for a comprehensive analysis.

MicrobeMod requires the user to align the BAM files generated by Dorado to the reference genome. Since the read extraction from the basecalled BAM format and the alignment to a reference FASTA are the first steps in our pipeline, you can directly run MicrobeMod using the intermediate file `methylation_mapped.bam` as follows:

```bash
MicrobeMod call_methylation -b <path_to_the_results_folder>/methylation_mapped.bam -r genome_reference.fasta
```
