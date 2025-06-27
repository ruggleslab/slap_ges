# SLAP-GES: Scoring Lupus Activity Transcriptomic Signature

## Overview

SLAP-GES (Systemic Lupus Activity Prediction - Gene Expression Signature) is a transcriptomic signature designed to assess lupus activity from gene expression data. This repository provides a tutorial for scoring the SLAP-GES signature using the `singscore` R package.

## What is SLAP-GES?

SLAP-GES is a bidirectional gene expression signature that captures the transcriptomic patterns associated with systemic lupus erythematosus (SLE) activity. The signature consists of:

- **Up-regulated genes**: 313 genes that are typically increased in active lupus
- **Down-regulated genes**: 235 genes that are typically decreased in active lupus

## Repository Structure

```text
slap_ges/
├── data/
│   ├── slap-ges_up.txt      # Up-regulated gene list (313 genes)
│   └── slap-ges_down.txt    # Down-regulated gene list (235 genes)
└── README.md               # This file
```

## Prerequisites

### Required R Packages

Install the following R packages before running the analysis:

```r
```r
# Check if packages are installed and install if needed
required_packages <- c("tidyverse", "magrittr")
bioc_packages <- c("singscore", "SummarizedExperiment")

# Check and install CRAN packages
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
    install.packages(missing_packages)
}

# Check and install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc)
}
```

### Data Requirements

Your gene expression data should be:

- A count matrix or normalized expression matrix
- Genes as rows, samples as columns
- Gene names/symbols as row names
- Preferably in a `SummarizedExperiment` object

## Tutorial: Scoring SLAP-GES Signature

### Step 1: Load Required Libraries

```r
library(tidyverse)
library(singscore)
library(SummarizedExperiment)
```

### Step 2: Load the SLAP-GES Gene Signatures

```r
# Load the up-regulated genes
up_genes <- read.table('data/slap-ges_up.txt', 
                       header = FALSE, 
                       stringsAsFactors = FALSE)$V1

# Load the down-regulated genes  
down_genes <- read.table('data/slap-ges_down.txt', 
                         header = FALSE, 
                         stringsAsFactors = FALSE)$V1

# Check signature sizes
cat("Up-regulated genes:", length(up_genes), "\n")
cat("Down-regulated genes:", length(down_genes), "\n")
```

### Step 3: Prepare Your Expression Data

```r
# Assuming you have a SummarizedExperiment object called 'se'
# If you have a count matrix, convert it to SummarizedExperiment:
# se <- SummarizedExperiment(assays = list(counts = your_count_matrix),
#                           colData = your_sample_metadata)

# Check gene overlap with signature
genes_in_data <- rownames(se)
up_overlap <- sum(up_genes %in% genes_in_data)
down_overlap <- sum(down_genes %in% genes_in_data)

cat("Up genes found in data:", up_overlap, "/", length(up_genes), 
    "(", round(up_overlap/length(up_genes)*100, 1), "%)\n")
cat("Down genes found in data:", down_overlap, "/", length(down_genes), 
    "(", round(down_overlap/length(down_genes)*100, 1), "%)\n")
```

### Step 4: Score the SLAP-GES Signature

```r
# Rank genes based on expression
ranked_genes <- rankGenes(se)

# Calculate signature scores using singscore
slap_ges_scores <- simpleScore(ranked_genes, 
                               upSet = up_genes, 
                               downSet = down_genes)

# Extract the total scores
signature_scores <- slap_ges_scores$TotalScore
names(signature_scores) <- colnames(se)
```

### Step 5: Add Scores to Sample Metadata

```r
# Add scores to colData
colData(se)$SLAP_GES_score <- signature_scores
colData(se)$SLAP_GES_up_score <- slap_ges_scores$UpScore
colData(se)$SLAP_GES_down_score <- slap_ges_scores$DownScore

# Convert to data frame for easier manipulation
sample_data <- as.data.frame(colData(se))
```

### Step 6: Normalize Scores (Optional)

If you want to normalize scores to a 0-1 range:

```r
# Min-max normalization
normalize_scores <- function(scores) {
  (scores - min(scores)) / (max(scores) - min(scores))
}

sample_data$SLAP_GES_score_normalized <- normalize_scores(sample_data$SLAP_GES_score)
```

### Step 7: Visualize Results

```r
library(ggplot2)

# Histogram of signature scores
p1 <- ggplot(sample_data, aes(x = SLAP_GES_score)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Distribution of SLAP-GES Scores",
       x = "SLAP-GES Score",
       y = "Count") +
  theme_bw()

# Boxplot by condition (assuming you have a 'condition' column)
# p2 <- ggplot(sample_data, aes(x = condition, y = SLAP_GES_score, fill = condition)) +
#   geom_boxplot() +
#   labs(title = "SLAP-GES Scores by Condition",
#        x = "Condition",
#        y = "SLAP-GES Score") +
#   theme_bw() +
#   theme(legend.position = "none")

print(p1)
```

## Understanding the Results

### Score Interpretation

- **Higher SLAP-GES scores**: Indicate higher lupus activity (up-regulated genes are highly expressed, down-regulated genes are lowly expressed)
- **Lower SLAP-GES scores**: Indicate lower lupus activity
- **Score range**: Typically ranges from negative to positive values, with the exact range depending on your dataset

### Score Components

- **TotalScore**: Combined score incorporating both up and down gene sets
- **UpScore**: Score based only on up-regulated genes
- **DownScore**: Score based only on down-regulated genes (note: this is typically negative)

## Troubleshooting

### Common Issues

1. **Low gene overlap**: If many signature genes are missing from your data:
   - Check gene naming conventions (HGNC symbols vs Ensembl IDs)
   - Consider gene mapping/conversion
   - Evaluate if your platform captures the signature genes

2. **Unexpected score distributions**:
   - Verify data normalization
   - Check for batch effects
   - Ensure proper sample QC

3. **Missing genes**: The scoring will work with partial gene sets, but consider the impact on signature performance

## Citation

Please cite the original SLAP-GES publication and the singscore package:

```text
# SLAP-GES citation (replace with actual publication details)
# [Citation for SLAP-GES signature]

# singscore citation
Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., & Davis, M. J. (2018). 
Single sample scoring of molecular phenotypes. BMC bioinformatics, 19(1), 1-10.
```

## Contact

For questions or issues with this tutorial, please contact [your contact information].

## Session Information

```r
sessionInfo()
```
