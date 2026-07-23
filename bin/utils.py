import re 
import pandas as pd

def parse_gff3_gene_names(gff3_file):
    """
    Parse bakta GFF3 and extract the gene names.

    Returns:
        list of tuples: (contig, start, end, strand, gene_name)
    """

    all_regions = []

    with open(gff3_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]

            if feature_type in ['CDS', 'rRNA', 'tRNA', 'oriC']:

                contig = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                gene_name = None

                # Priority order
                for key in ["Name=", "gene=", "gene_name=", "locus_tag="]:
                    match = re.search(f"{key}([^;]+)", attributes)
                    if match:
                        gene_name = match.group(1)
                        break

                if gene_name is None:
                    gene_name = "unknown"

                if feature_type == 'rRNA':
                    gene_name = f"rRNA_{gene_name}"
                elif feature_type == 'tRNA':
                    gene_name = f"tRNA_{gene_name}"
                elif feature_type == 'oriC':
                    gene_name = "oriC"

                all_regions.append((contig, start, end, strand, gene_name))

    return all_regions

def read_modkit(modkit_output, percent_modified_threshold, modification):
    """Parse the output of Modkit and save it in a pd dataframe. 
    """

    d = pd.read_csv(modkit_output, sep="\t", header=None)

    ## Rename columns
    d.columns = [
        "Contig",
        "Position",
        "End",
        "Modification",
        "drop0",
        "Strand",
        "drop1",
        "drop2",
        "drop3",
        "Valid_coverage",  # Nmod + Nother_mod + Ncanonical
        "drop5",
        "Modified_bases",
        "Unmodified_bases",
        "Other_mod_base",
        "drop6",
        "Modification_below_threshold",
        "Other_bases",
        "drop7",
    ]
    d = d[d.columns[~d.columns.str.contains("drop")]]
    
    # total coverage is valid coverage + bases that didn't pass the modification threshold. Used instead of valid!
    d["Total_coverage"] = (
        d["Modified_bases"]
        + d["Unmodified_bases"]
        + d["Other_mod_base"]
        + d["Modification_below_threshold"]
        + d["Other_bases"]
    )

    d["Percent_modified"] = d["Modified_bases"] / d["Total_coverage"]  # here computing the percent modified based on total coverage
    d["SNP_Position"] = d["Contig"].astype(str) + ":" + d["Position"].astype(str)

    d["SNP_Position_Strand"] = d["SNP_Position"] + "-" + d["Strand"]

    d = d[(d.Percent_modified >= percent_modified_threshold) & (d.Modification == modification) & (d.Total_coverage > 0)].drop("End", axis=1) 

    return d
