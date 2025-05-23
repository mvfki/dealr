% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dealr.R
\name{print.dealr}
\alias{print.dealr}
\title{Show collapsed information of DEALR result}
\usage{
\method{print}{dealr}(x, trunc = 10, ...)
}
\arguments{
\item{x}{A \code{dealr} object.}

\item{trunc}{Number of maximum list element to show.}

\item{...}{Additional arguments (not used).}
}
\description{
Show collapsed information of DEALR result
}
\examples{
db_mini <- db$mouse[db$mouse$pathway_name \%in\% c('IL16', 'TNF'),]
lr <- dealr(deg_mini, db_mini)
print(lr)
}
