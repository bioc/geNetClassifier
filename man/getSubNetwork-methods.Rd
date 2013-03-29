\name{getSubNetwork-methods}
\docType{methods}
\alias{getSubNetwork}
\alias{getSubNetwork-methods}
\alias{getSubNetwork,GenesNetwork-method}
\alias{getSubNetwork,GeNetClassifierReturn-method}
\alias{getSubNetwork,list-method}
\alias{getSubNetwork,NULL-method}
\title{Get a sub-network.}
\description{
Returns the sub-network formed by the given genes.
}
\usage{
getSubNetwork(network, genes, showWarnings=TRUE)
}
\arguments{
 \item{network}{
GenesNetwork or GenesNetwork list containing the whole network.
}
 \item{genes}{
GenesRanking or character vector. Genes in the new network.
}
 \item{showWarnings}{
Logical. If true, shows warnings if the given genes are not in the network.
}
}
\value{
A \code{\link[geNetClassifier:GenesNetwork-class]{GenesNetwork}} or list of networks between the given genes.
}
				
\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's class (\code{\link[geNetClassifier:GenesNetwork-class]{GenesNetwork}}) help page.
}
\examples{
data(leukemiasClassifier)
clGenesSubNet <- getSubNetwork(leukemiasClassifier@genesNetwork, 
leukemiasClassifier@classificationGenes)
getSubNetwork(leukemiasClassifier@genesNetwork, getTopRanking(leukemiasClassifier@genesRanking, numGenesClass=30))
}
\keyword{methods}





