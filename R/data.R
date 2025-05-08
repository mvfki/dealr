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
