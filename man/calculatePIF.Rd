% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rRIF.R
\name{calculatePIF}
\alias{calculatePIF}
\title{calculates phenotype impact factor (PIF) from transcriptomic data}
\usage{
calculatePIF(eSet = eSet, target.factor = target.factor, A, B)
}
\arguments{
\item{eSet}{an ExpressionSet}

\item{target.factor}{a factor of interest}

\item{A}{a level of a factor (e.g., condition 1)}

\item{B}{a level of a factor (e.g., condition 2)}
}
\value{
numeric
}
\description{
This function calculates phenotype impact factor (PIF).
}
\examples{
data(ToniData)

formula <- geno1~geno2
samp <- list(A=as.character(formula[[2]]), B=as.character(formula[[3]]))
target.factor <- "Genotype"
PIF <- calculatePIF(eSet=ToniData, target.factor=target.factor,
                    samp$A, samp$B)

}
\references{
Reverter A et al. Bioinformatics 26:896 (2010)
(\href{https://www.ncbi.nlm.nih.gov/pubmed/20144946}{PubMed})
}
\author{
Kevin Rue-Albrecht \url{https://github.com/kevinrue/HudsonRIF}
cosmetical changes by Atsushi Fukushima
}
