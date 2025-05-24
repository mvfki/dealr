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

We have implemented a DE test function solely for this purpose. it works for
either a plain matrix or list of matrices bundled with metadata variables, or
commonly seen container objects like `Seurat` and `liger`.

For A Seurat or liger object:

```r
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

1. Load the package

```r
library(dealr)
```

2. Load the ligand-receptor database. 

```r
data(db)
```

The object `db` is a list object with two elements `db$mouse` and `db$human`.
They are identical to `CellChat::CellChatDB.mouse$interaction` and 
`CellChat::CellChatDB.human$interaction`, respectively (copied from v2.2.0).

3. Run LR-pair level inference

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

>This might be improved in the future introducing expression percentage obtained
from single-cell data instead of just pseudo-bulk data.

4. Proceed to pathway level inference

```r
dealr_pe_result <- pathwayEnrich(dealr_result)
```

### 3. Visualize the result

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
