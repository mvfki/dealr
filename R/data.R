#' CellChatDB LR-pair databases
#' @description
#' The LR-pair interaction database directly taken from CellChat v2.2.0.
#' Information is minimized to only necessary columns for the analysis.
#' @format
#' A named list of two data.frames - \code{db$mouse} and \code{db$human}. Both
#' has the \code{rownames()} as unique LR-pair identifiers. Both have the same
#' columns:
#' \itemize{
#' \item{\code{pathway_name} - Name of the pathway that LR-pair belongs to.}
#' \item{\code{ligand.symbol} - A string concatenating all gene symbols of the
#' ligand, separated by \code{', '}}
#' \item{\code{receptor.symbol} - A string concatenating all gene symbols of the
#' receptor, separated by \code{', '}}
#' }
'db'

#' Example DEG list input
#' @description
#' This is a minimum subset of real analysis result from our in-house dataset.
#' Pseudo-bulk DESeq2 tests were performed between treated and control samples
#' within each cell type. This object contains the result for only three cell
#' types (CD4_T, CD8_T and NK), for minimum demonstration purpose. The DESeq2
#' result objects were coerced to plain `data.frame` and only contain test
#' statistics for 5 genes. These genes were chosen based on small pathways
#' (IL16 and TNF) that generally show significant intensity changes across the
#' three cell types.
#' @format A list of three data.frames, named as `c('CD4_T', 'CD8_T', 'NK')`.
#' Each data.frame has the same 5 rows, corresponding to the 5 genes. Each
#' data.frame has the same columns and these are the exact output of DESeq2.
"deg_mini"
