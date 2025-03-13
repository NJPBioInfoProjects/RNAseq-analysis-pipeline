library(tidyverse)
library(DESeq2)

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

sample_info <- data.frame(
  row.names = colnames(filtered_counts),
  condition = c("control", "control", "control", "exp", "exp", "exp")  # Adjust based on samples
)

# Ensure 'condition' is a factor (not character)
sample_info$condition <- as.factor(sample_info$condition)

# Create DESeq2 dataset from filtered counts
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = sample_info, design = ~ condition)

# Apply variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)  # or use rlog(dds, blind = TRUE)

# Extract normalized counts
normalized_counts <- assay(vsd)

# Save normalized counts to CSV
# write_csv(as.data.frame(normalized_counts) %>% rownames_to_column("gene"),
#           "normalized_counts.csv")

# Now, let's do PCA
# Load necessary libraries
library(ggplot2)

# Transpose the normalized counts matrix (PCA requires samples as rows)
pca_data <- prcomp(t(normalized_counts), scale. = TRUE)

# Extract PCA results
pca_df <- as.data.frame(pca_data$x)
pca_df$Sample <- rownames(pca_df)  # Add sample names

# Add condition labels (modify accordingly)
pca_df$Condition <- ifelse(grepl("control", pca_df$Sample), "Control", "Experiment")

# View PCA variance explained
print(summary(pca_data))

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  labs(title = "PCA of RNA-seq Data",
       x = paste0("PC1: ", round(100 * summary(pca_data)$importance[2,1], 2), "% Variance"),
       y = paste0("PC2: ", round(100 * summary(pca_data)$importance[2,2], 2), "% Variance")) +
  theme_minimal() +
  theme(legend.position = "top")

print(p)