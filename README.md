# dealr: Differential Expression Aware Ligand-Receptor Inference

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Cell-cell interaction inference from single-cell transcriptomics data has been 
popular, but the comparison of the signaling between conditions has been an 
awkward field. Current tools tend to score LR-pairs within single samples, which
allows in-sample comparison, while it is less convincing to directly compare 
the scores across conditions consisting of multiple samples.

DEALR starts from DEG analysis and take the statistics for differential LR-pair
signaling analysis. The DEG statistics provides confidence for the inference
of signaling changes. We provide wrapped-up pipeline to perform the DEG analysis
that fit into DEALR inference, as well as visualization method to examine if 
the differential signaling is true to intuition.

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

We have implemented a DE test function solely for this purpose. it works for
either a plain matrix or list of matrices bundled with metadata variables, or
commonly seen container objects like `Seurat` and `liger`.

For A Seurat or liger object:

```r
library(dealr)

degList <- pseudobulkDE(
    data = seurat_or_liger_Obj,
    condVar = 'treatment', condTest = 'treated', condCtrl = 'control', 
    splitVar = 'cell_type', 
    replicateVar = 'sample'
)
```

`condVar`, `splitVar` and `replicateVar` must be available in the metadata slot.
`splitVar` and `replicateVar` will be set to conventional variables if left as
default. That is `Idents(seuratObj)` and `seuratObj$orig.ident` for Seurat, and,
`defaultCluster(ligerObj)` and `ligerObj$dataset` for liger, respectively.

Please see function manual for more details when using matrix or list of 
matrices pulled out from other forms of data structure.

### 2. Run DEALR inference

1. Load the ligand-receptor database. 

```r
data(db)
```

The object `db` is a list object with two elements `db$mouse` and `db$human`.
They are identical to `CellChat::CellChatDB.mouse$interaction` and 
`CellChat::CellChatDB.human$interaction`, respectively (copied from v2.2.0).

2. Run LR-pair level inference

```r
dealr_result <- dealr(degList, db$mouse)
```

`degList` is currently expected to be a named list of data.frame objects, each
returned by a DESeq2 run. See previous section for detail. 

Parameter `baseMeanThresh` can be adjusted to filter out lowly expressed genes
to deny that they are expressed enough for forming LR interaction. This is 
currently based on the `baseMean` field of DESeq2 result, which is the average 
of the normalized counts across all samples. In our case, the average of the
DESeq2-normalized pseudo-bulk counts from all samples involved in a test. 

Parameter `abslogfcThresh` can also be adjusted to filter out LR-pairs without
largely differentially expressed components. 

>This might be improved in the future introducing expression percentage obtained
from single-cell data instead of just pseudo-bulk data.

3. Proceed to pathway level inference

```r
dealr_pe_result <- pathwayEnrich(dealr_result)
```

### 3. Visualize the result

To explore the result in a wider range, we provide `plotLRPairDot()` and 
`plotPathwayEnrichDot()` to roughly show the top significant LR-pairs in a
range of sender-receiver pairs. For a precise evidence at single-cell gene 
expression level, we provide `plotLRGeneHeatmap()` to exactly show the 
expression change of the LR-pair components.

#### Broad exploration with dot plot of statistics

- `plotLRPairDot(dealr_result)`: Dot plot of LR-pair significance between 
possible sender-receiver pair. Larger dot size indicates more significant 
LR-pair. Red color (default for `upreg_col`) indicates up-regulated LR-pair 
signaling and blue color (default for `downreg_col`) indicates down-regulated
LR-pair signaling. The color intensity indicates the significance level as well.
Note that there is no difference between the significance level represented by
size and color. Suggestions for more metrics to be shown are very welcomed.
- `plotPathwayEnrichDot(dealr_pe_result)`: Dot plot of pathway significance 
between possible sender-receiver pair. Larger dot size indicates more LR-pairs
are significant within the pathway, while metric for this can vary. Use 
`size_by = 'overlap'` (default) for exact number of LR-pairs or `'enrichment'` 
for using a relative enrichment ratio. Red color (default for `upreg_col`) 
indicates up-regulated pathway signaling and blue color (default for
`downreg_col`) indicates down-regulated pathway signaling. Darker color 
indicates more significant pathways. 

The interpretaion of the result is quite complex, given the network between 
three levels of information: senders, receivers and LR-pairs/pathways. The 
functions above visualize all significant inference though, we have the 
filtering logic to filter the senders and receivers to allow focused
panels for specific identities. 

- `sender_use`: Include sender-receiver pairs with the sender(s) of interest. 
One or more.
- `receiver_use`: Include sender-receiver pairs with the receiver(s) of 
interest. One or more.
- `focus`: Include sender-receiver pairs with having this identity as either 
the sender or the receiver. Only one allowed.

For example, when having cell type A, B and C in the analysis, a parameter set
of `sender_use = c('A', 'B'), receiver_use = c('A', 'B'), focus = 'A'` will
show you the significan LR-pairs/pathways in A->A, A->B and B->A. B->B, A->/<-C
and B->/<-C will be filtered out.

More examples are shown below.

Show differential signaling LR-pairs sent from a certain cell type, received by 
all possible cell types

```r
plotLRPairDot(dealr_result, sender_use = 'Mac')
```

Show all differential signaling communication involving a cell type of interests

```r
plotLRPairDot(dealr_result, focus = 'Mac')
```

Show the enriched differential signaling pathways sent from all possible cell 
types, received by a certain cell type

```r
plotPathwayEnrichDot(dealr_pe_result, receiver_use = 'Mac')
```

Show all enriched differential signaling pathways involving a cell type of 
interests, but only between this cell type and other immune cells.

```r
immunes <- c('Mono', 'Mac', 'DC', 'B', 'T4', 'T8', 'NK') # ... as many as you have
plotPathwayEnrichDot(dealr_pe_result, focus = 'Mac', sender_use = immunes, 
                     receiver_use = immunes)
```

#### Precise examination of gene expression between a certain pair of sender and receiver

- `plotLRGeneHeatmap(data_obj, dealr_result, ...)`: This makes a heatmap of
gene expression in the original single cells to show expression change of 
components of an LR-pair between conditions and from sender to receiver. We put
three panels of heatmaps together: one on the left for the sender cell type,
combining cells from all samples and grouped by conditions, one on the right for
the receiver cell type, also grouped by conditions, and one in the middle 
showing the statistics of the LR-pairs. For the left or right pannel, each row 
is a gene, and each column is a cell. Genes on the same row belong to the same 
LR-pair. LR-pairs with complex components, e.g. A-(B+C), are expanded to 
multiple rows, e.g. A-B and A-C. LR-pairs on the rows are further grouped by 
pathways.

For this function, users can only set one `sender_use` and one `receiver_use`. 
Since this visualization method incorporates the original single-cell data, the 
cell-level variables used for the psuedo-bulk DEG must be provided identically (
i.e. use the same `splitVar`, `condVar`, `condTest` and `condCtrl` if you used
the `pseudobulkDE()` pipeline, or otherwise equivalent if you constructed 
`degList` yourself). 

> Engineering of the software might be improved in the future to eliminate the 
need that users have to make sure they use the identical input. 

To show 1. the statistics of top significant LR-pairs between specific 
sender-receiver pair of cell types; 2. the corresponding ligand gene expression 
from the sender cell type, with contrast between conditions; 2. and the receptor 
gene expression from the receiver cell type, also with contrast between 
conditions, all at the same time:

```r
plotLRGeneHeatmap(
    data = seurat_or_liger_Obj,
    dealr = dealr_result,
    splitVar = 'cell_type', sender_use = 'DC', receiver_use = 'T',
    condVar = 'treatment', condTest = 'treated', condCtrl = 'control'
)
```
