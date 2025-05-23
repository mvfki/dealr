% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dealr.R
\name{pathwayEnrich}
\alias{pathwayEnrich}
\title{Calculate pathway enrichment from DEALR result}
\usage{
pathwayEnrich(dealr, p_thresh = 0.01)
}
\arguments{
\item{dealr}{A \code{dealr} object. DO NOT FILTER.}

\item{p_thresh}{Numeric threshold on individual inference to determine if
the LR-pair is significant between the sender and receiver. Default
\code{0.01}.}
}
\value{
A \code{dealr_pe} object, inherited from \code{tbl_df}, \code{dbl} and
\code{data.frame}. It contains the following columns:
\itemize{
\item{\code{sender} - Group name where the ligand statistics are from}
\item{\code{receiver} - Group name where the receptor statistics are from}
\item{\code{pathway} - Pathway name}
\item{\code{stat} - Combined z-score of the pathway between the sender and
receiver}
\item{\code{overlap} - Number of significant LR-pairs in the
pathway, between the sender and receiver}
\item{\code{signif_n_intr} - Number of significant LR-pairs between the
sender and receiver}
\item{\code{pathway_size} - Number of LR-pairs in the pathway, summarized
from \code{db} but not yielded from the analysis}
\item{\code{p} - P-value of the combined z-score}
\item{\code{enrichment} - Enrichment score of the pathway}
}
}
\description{
Downstream of the DEALR analysis, this function further summarizes if any
pathway is significantly enriched between the sender and receiver. It
combines the z-scores of the LR-pairs in the same pathway to get the z-score
of the pathway. It also calculates the enrichment score of the pathway
between each pair of sender and receiver with
\deqn{\frac{(n\>significant\>LRpairs\>in\>pathway)/(size\>of\>pathway)}
{(n\>significant\>LRpairs)/(size\>of\>database)}\times{\sqrt{size\>of\>pathway}}}
}
\examples{
db_mini <- db$mouse[db$mouse$pathway_name \%in\% c('IL16', 'TNF'),]
lr <- dealr(deg_mini, db_mini)
pe <- pathwayEnrich(lr)
}
