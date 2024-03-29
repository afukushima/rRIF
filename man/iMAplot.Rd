% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{iMAplot}
\alias{iMAplot}
\title{an interactive MA-plot of rRIF object}
\usage{
iMAplot(rRIF)
}
\arguments{
\item{rRIF}{"rRIF" class object}
}
\value{
none
}
\description{
This function provides an interactive MA-plot of rRIF results
Differential expression of all genes plotted against the average abundance
}
\examples{
data(ToniData)
data(ToniData.DEGs)
data(ToniData.TFs)

res <- rRIF(eSet = ToniData,
            formula = geno1~geno2,
            target.factor = "Genotype",
            DEGs = ToniData.DEGs,
            cor.method="pearson",
            regulator.list = ToniData.TFs
            )
iMAplot(res)
}
