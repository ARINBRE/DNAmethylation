# DNAmethylation

## Overview
This project is designed to preprocess, analyze, and visualize DNA methylation data using the Illumina Infinium methylation array 450k. The workflow includes downloading raw data from the Gene Expression Omnibus (GEO), preprocessing the data, performing differential methylation analysis, and generating plots to visualize the results.

## Workflow Steps
1. **Data Acquisition**:
   - Download raw DNA methylation data using the `GEOquery` package.
   - Decompress and prepare the dataset for analysis.

2. **Data Preprocessing**:
   - Read and preprocess the raw `.idat` files using the `minfi` package.
   - Perform quality control, including probe quality assessment and removal of poor-quality probes.
   - Normalize data and compute methylation metrics.

3. **Differential Methylation Analysis**:
   - Annotate probes and perform differential methylation analysis using the `limma` package.
   - Extract statistically significant probes based on adjusted p-values and log fold changes.

4. **Visualization**:
   - Generate a volcano plot to highlight significant methylation changes.
   - Create a PCA scatter plot to visualize the data distribution among samples.

5. **Result Compilation**:
   - Save the analysis results in a CSV file for further examination and validation.

## Setup and Installation
Ensure that R is installed on your machine along with the required libraries. You can install the necessary R packages using the following commands:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GEOquery", "minfi", "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "limma"))
