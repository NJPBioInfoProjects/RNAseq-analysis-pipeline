#!/usr/bin/env python3

import argparse
import re

# Initialize argparse
parser = argparse.ArgumentParser()

# Define input and output arguments
parser.add_argument("-i", "--input", help="Input GTF file", dest="input", required=True)
parser.add_argument("-o", "--output", help="Output file name and path with Ensembl Gene IDs and Gene Symbols", dest="output", required=True)

# Parse arguments
args = parser.parse_args()

# Open input GTF and output file
with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
    outfile.write("Ensembl_Gene_ID\tGene_Symbol\n")

    for line in infile:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if fields[2] != "gene":
            continue

        # Extract gene_id and gene_name
        gene_id_match = re.search(r'gene_id "([^"]+)"', fields[8])
        gene_name_match = re.search(r'gene_name "([^"]+)"', fields[8])

        if gene_id_match and gene_name_match:
            outfile.write(f"{gene_id_match.group(1)}\t{gene_name_match.group(1)}\n")

print(f"Parsing complete. Output saved to {args.output}")