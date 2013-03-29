\name{numSignificantGenes-methods}
\docType{methods}
\alias{numSignificantGenes}
\alias{numSignificantGenes-methods}
\alias{numSignificantGenes,GenesRanking-method}
\alias{numSignificantGenes.GenesRanking}
\title{Number of ranked genes over the posterior probability threshold.}
\description{
Provides the number of ranked genes over the posterior probability threshold
}
% \usage{
% \method{numSignificantGenes}{GenesRanking}(object, lpThreshold=0.75, numSignificantGenesType="ranked")
% }
\section{Methods}{
\describe{
\code{numSignificantGenes(object, lpThreshold=0.75, numSignificantGenesType="ranked") 
}
}}
\arguments{
 \item{object}{a GenesRanking}
 \item{lpThreshold}{Posterior probability threshold}
 \item{numSignificantGenesType}{"ranked" or "global". \cr
 Ranked will show the count of genes on the ranking of each class. Each gene will be counted only once, since it is only kept in the class for which it had better ranking. \cr
 Global counts the genes over the threshold before assigning them to a class. i.e. a gene might have 0.3 for one class, and 0.25 for another, if we are taking a thershold of 0.20, it will be counted on both classes.}
}
\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's class (\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}}) help page.
}
\examples{
data(leukemiasClassifier)
# Total number of genes in the ranking:
numGenes(leukemiasClassifier@genesRanking)
# Number of genes ovher the posterior probability threshold
# Default:  lpThreshold=0.75, numSignificantGenesType="ranked"
numSignificantGenes(leukemiasClassifier@genesRanking)
numSignificantGenes(leukemiasClassifier@genesRanking, numSignificantGenesType="global")
numSignificantGenes(leukemiasClassifier@genesRanking, lpThreshold=0.90)
}
\keyword{methods}










