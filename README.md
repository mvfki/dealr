# dealr: Differentially Expression Aware Ligand-Receptor Inference

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Cell-cell interaction inference from single-cell transcriptomics data has been popular, but the comparison of the signaling between conditions has been an awkward field. 

Basically, we aim to derive the statistics based on DE tests so it provides evident inference. Please see function manuals for detail description.

This package can only be installed from github for now. Run the command below in R:

```r
if (!requireNamespace("remotes", quietly = TRUE)) 
  install.packages("remotes")
remotes::install_github("mvfki/dealr")
```

## Usage

### 1. Prepare Pseudo-bulk DESeq2 result

This tool currently only targets at scRNAseq dataset, with the study design 
consisting of multiple conditions each having multiple replicates. The 
prerequisite is the pseudo-bulk DESeq2 tests comparing two conditions within 
each group of cell. For example,

1. Optionally, subset the dataset to only the two conditions of interests if 
there are more;
2. Split the subset by cell type and do the following step for each cell type;
3. Aggregate (sum up) raw counts by replicate. This builds the pseudo-bulk. This 
naturally yields replicate-level metadata of for the conditions of interests;
4. Run DESeq2 Wald test on the pseudo-bulk data with the design on the condition.

**Make sure to have all tests contrasting the same direction.**

The procedure above should results in a number of DESeq2 result data frames. 
This function expect a named list object that gathers all these data frames. 
List names are for the cell types. **DO NOT FILTER** the DESeq2 result data 
frames before passing to the analysis.

### 2. Run DEALR inference

1. Load the package

```r
library(dealr)
```

2. Load the ligand-receptor database

```r
data(db)
```

3. Run LR-pair level inference

```r
dealr_result <- dealr(degList, db$mouse)
```

4. Proceed to pathway level inference

```r
dealr_pe_result <- pathwayEnrich(dealr_result)
```

### 3. Visualize the result

Show differential signaling LR-pairs sent from a certain cell type, received by all possible cell types

```r
plotLRPairDot(lr, sender_use = 'Mac')
```

Show differential signaling LR-pairs received by a certain cell type, sent from all possible cell types

```r
plotLRPairDot(lr, receiver_use = 'Mac')
```

Show all differential signaling communication involving a cell type of interests

```r
plotLRPairDot(lr, focus = 'Mac')
```

Show the enriched differential signaling pathways sent from a certain cell type, received by all possible cell types

```r
plotPathwayEnrichDot(dealr_pe_result, sender_use = 'Mac')
```

Show the enriched differential signaling pathways sent from all possible cell types, received by a certain cell type

```r
plotPathwayEnrichDot(dealr_pe_result, receiver_use = 'Mac')
```

Show all enriched differential signaling pathways involving a cell type of interests

```r
plotPathwayEnrichDot(dealr_pe_result, focus = 'Mac')
```
