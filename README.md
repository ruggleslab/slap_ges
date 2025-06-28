# SLE Activity Platelet - Gene Expression Signature

## Overview

SLAP-GES (SLE Activity Platelet - Gene Expression Signature) is a transcriptomic signature designed to assess lupus activity from gene expression data. This repository provides a tutorial for scoring the SLAP-GES signature using the `singscore` R package.

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
required_packages <- c("tidyverse")
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

- **Higher SLAP-GES scores**: Associated with higher lupus activity, increased platelet aggregation and impaired macro- and micro-vascular function. Altogether, higher SLAP-GES suggest higher cardiovascular risk (up-regulated genes are highly expressed, down-regulated genes are lowly expressed).
- **Score range**: Typically ranges from negative to positive values, with the exact range depending on your dataset

### Score Components

- **TotalScore**: Combined score representing **SLAP-GES**
- **UpScore**: Score based only on up-regulated genes
- **DownScore**: Score based only on down-regulated genes (note: this is typically negative)

## Troubleshooting

### Common Issues

1. **Low gene overlap**: If all signature genes are missing from your data:
   - Check gene naming conventions (HGNC symbols vs Ensembl IDs)
   - Make sure mapping is to hg38
   - Evaluate if your platform captures the signature genes

2. **Unexpected score distributions**:
   - Verify data normalization (variance stabilization from DESeq2 is recommended)
   - Check for batch effects
   - Ensure proper sample QC

3. **Missing genes**: The scoring will work with partial gene sets, but consider the impact on signature performance

## Citation

Please cite the original SLAP-GES publication and the singscore package:

```text
# SLAP-GES citation
# [TBD]

# singscore citation
Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., & Davis, M. J. (2018). 
Single sample scoring of molecular phenotypes. BMC bioinformatics, 19(1), 1-10.
```

## Contact

For questions or issues with this tutorial, please contact [mm12865@nyu.edu].

## Session Information

```r
R version 4.4.3 (2025-02-28)
Platform: aarch64-apple-darwin24.2.0
Running under: macOS Sequoia 15.5

Matrix products: default
BLAS:   /opt/homebrew/Cellar/openblas/0.3.29/lib/libopenblasp-r0.3.29.dylib 
LAPACK: /opt/homebrew/Cellar/r/4.4.3_1/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] SummarizedExperiment_1.36.0 Biobase_2.66.0             
 [3] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
 [5] IRanges_2.40.1              S4Vectors_0.44.0           
 [7] BiocGenerics_0.52.0         MatrixGenerics_1.18.1      
 [9] matrixStats_1.5.0           singscore_1.26.0           
[11] lubridate_1.9.4             forcats_1.0.0              
[13] stringr_1.5.1               dplyr_1.1.4                
[15] purrr_1.0.4                 readr_2.1.5                
[17] tidyr_1.3.1                 tibble_3.2.1               
[19] ggplot2_3.5.1               tidyverse_2.0.0            

loaded via a namespace (and not attached):
 [1] KEGGREST_1.46.0         gtable_0.3.6            lattice_0.22-6         
 [4] tzdb_0.4.0              vctrs_0.6.5             tools_4.4.3            
 [7] generics_0.1.3          AnnotationDbi_1.68.0    RSQLite_2.3.9          
[10] blob_1.2.4              pkgconfig_2.0.3         Matrix_1.7-2           
[13] graph_1.84.1            lifecycle_1.0.4         GenomeInfoDbData_1.2.13
[16] compiler_4.4.3          Biostrings_2.74.1       statmod_1.5.0          
[19] munsell_0.5.1           pillar_1.10.1           crayon_1.5.3           
[22] limma_3.62.2            DelayedArray_0.32.0     cachem_1.1.0           
[25] abind_1.4-8             locfit_1.5-9.11         tidyselect_1.2.1       
[28] stringi_1.8.4           reshape2_1.4.4          fastmap_1.2.0          
[31] grid_4.4.3              SparseArray_1.6.1       colorspace_2.1-1       
[34] cli_3.6.5               magrittr_2.0.3          S4Arrays_1.6.0         
[37] XML_3.99-0.17           GSEABase_1.68.0         edgeR_4.4.2            
[40] withr_3.0.2             scales_1.3.0            UCSC.utils_1.2.0       
[43] bit64_4.6.0-1           timechange_0.3.0        XVector_0.46.0         
[46] httr_1.4.7              bit_4.6.0               png_0.1-8              
[49] hms_1.1.3               memoise_2.0.1           rlang_1.1.6            
[52] Rcpp_1.0.14             xtable_1.8-4            glue_1.8.0             
[55] DBI_1.2.3               annotate_1.84.0         jsonlite_2.0.0         
[58] plyr_1.8.9              R6_2.6.1                zlibbioc_1.52.0   ```
