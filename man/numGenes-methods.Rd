\name{numGenes-methods}
\docType{methods}
\alias{numGenes}
\alias{numGenes-methods}
\alias{numGenes,GenesRanking-method}
\title{Number of genes in the genesRanking.}
\description{
Provides the number of genes in the genesRanking.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "GenesRanking")}}{
}
}}
\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's class (\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}}) help page.
}
\examples{
data(leukemiasClassifier)
numGenes(leukemiasClassifier@genesRanking)
}
\keyword{methods}