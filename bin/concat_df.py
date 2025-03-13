#!/usr/bin/env python3

import argparse
import pandas as pd
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Concatenate multiple exon count files into a single matrix.")
parser.add_argument("-i", "--input", nargs="+", required=True, help="List of exon count files to concatenate.")
parser.add_argument("-o", "--output", required=True, help="Output CSV file name.")
args = parser.parse_args()

# Sort files to maintain consistent order
sorted_files = sorted(args.input)

# Read and merge exon count files
df_list = []
for file in sorted_files:
    df = pd.read_csv(file, sep="\t", index_col=0, header=0)  # Adjust separator if necessary
    df.columns = [os.path.basename(file).replace(".exon.txt", "")]  # Use sample name as column header
    df_list.append(df)

# Concatenate along columns
merged_df = pd.concat(df_list, axis=1)

# Save output
merged_df.to_csv(args.output)

print(f"Concatenation complete. Output saved to {args.output}")