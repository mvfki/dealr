% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{db}
\alias{db}
\title{CellChatDB LR-pair databases}
\format{
A named list of two data.frames - \code{db$mouse} and \code{db$human}. Both
has the \code{rownames()} as unique LR-pair identifiers. Both have the same
columns:
\itemize{
\item{\code{pathway_name} - Name of the pathway that LR-pair belongs to.}
\item{\code{ligand.symbol} - A string concatenating all gene symbols of the
ligand, separated by \code{', '}}
\item{\code{receptor.symbol} - A string concatenating all gene symbols of the
receptor, separated by \code{', '}}
}
}
\usage{
db
}
\description{
The LR-pair interaction database directly taken from CellChat v2.2.0.
Information is minimized to only necessary columns for the analysis.
}
\keyword{datasets}
