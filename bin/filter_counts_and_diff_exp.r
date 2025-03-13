#this is to filter the counts matrix

library(DESeq2)
library(tidyverse)

# Load the counts matrix, ensuring gene names are set as row names
counts <- read_csv("results/verse_concat.csv") %>%
  column_to_rownames(var = "gene")  # Ensure gene names are used as row names

# Convert counts to numeric (in case they are read as characters)
counts <- counts %>%
  mutate(across(everything(), as.numeric))

# Count genes before filtering
genes_before <- nrow(counts)

# Filtering step 1: Remove genes with all-zero counts
filtered_counts <- counts[rowSums(counts) > 0, ]

# Filtering step 2: Keep genes expressed in at least 3 samples with at least 1 count
filtered_counts <- filtered_counts[rowSums(filtered_counts > 1) >= 3, ]

# Count genes after filtering
# genes_after <- nrow(filtered_counts)

# Report results
# cat("Genes before filtering:", genes_before, "\n")
# cat("Genes after filtering:", genes_after, "\n")

# par(mfrow=c(1,2))
# hist(log10(rowSums(counts) + 1), main="Before Filtering", xlab="Log10 Sum", col="steelblue")
# hist(log10(rowSums(filtered_counts) + 1), main="After Filtering", xlab="Log10 Sum", col="darkred")

sample_info <- data.frame(
  row.names = colnames(filtered_counts),
  condition = c("control", "control", "control", "exp", "exp", "exp")  # Adjust based on samples
)

# Ensure 'condition' is a factor (not character)
sample_info$condition <- as.factor(sample_info$condition)

# Create DESeq2 dataset from filtered counts
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = sample_info, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Convert results to dataframe and order by padj (adjusted p-value)
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene") %>%
  arrange(padj)  # Sort by adjusted p-value (most significant first)

id_map <- read_tsv("results/id2name.txt", col_names = c("gene", "gene_name"))

res_df_named <- res_df %>%
  left_join(id_map, by = "gene") %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene, gene_name))


# Select only the relevant columns (gene name, log2FoldChange, padj)
top10_genes_named <- res_df_named %>%
  slice(1:10) %>%
  select(gene_name, padj)  # Removed 'gene' column

# Save updated top 10 genes table
write_csv(top10_genes_named, "top10_significant_genes_named.csv")

# Print the table
print(top10_genes_named)

# next, select a p value threshold and report the number of genes that meet the criterion
padj_threshold <- 0.05

num_significant_genes <- sum(res_df_named$padj < 0.05, na.rm = TRUE)

#report the results
cat("Number of significant genes at padj <", padj_threshold, ":", num_significant_genes, "\n")

# Extract significant genes (padj < 0.05)
significant_genes <- res_df_named %>%
  filter(padj < 0.05) %>%
  pull(gene_name)  # Extract only gene names

# Save list as a text file (one gene per line)
write_lines(significant_genes, "significant_genes_list.txt")

cat("Significant genes saved to significant_genes_list.txt\n")