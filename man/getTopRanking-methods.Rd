\name{getTopRanking-methods}
\docType{methods}
\alias{getTopRanking}
\alias{getTopRanking-methods}
\alias{getTopRanking,GenesRanking-method}
\alias{getTopRanking.GenesRanking}
\title{Gets a new ranking with the given top genes.}
\description{
Returns a new ranking containing only the top genes of each class.
}
% \usage{
% \method{getTopRanking}{GenesRanking}(object, numGenesClass)
% }
\section{Methods}{
\describe{
\code{getTopRanking(object, numGenesClass)
}
}}
\arguments{
 \item{object}{a GenesRanking}
 \item{numGenesClass}{integer. Number of genes per class.}
 }
\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's class (\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}}) help page.
}
\examples{
data(leukemiasClassifier) # Sample classifier

# Sub-ranking with the top 10 genes:
getTopRanking(leukemiasClassifier@genesRanking, 10)
}
\keyword{methods}








