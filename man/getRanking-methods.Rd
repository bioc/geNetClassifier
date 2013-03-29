\name{getRanking-methods}
\docType{methods}
\alias{getRanking}
\alias{getRanking-methods}
\alias{getRanking,GenesRanking-method}
\alias{getRanking.GenesRanking}
\title{Shows the genes ranking.}
\description{
Shows the ranking as matrix: Ranked genes by classes.
}
% \usage{
%\method{getRanking}{GenesRanking}(object, showGeneLabels=TRUE, showGeneID=FALSE)
% }
\arguments{
\item{object}{a GenesRanking}
\item{showGeneID}{boolean. If TRUE, the genes will be shown with the gene IDs used in the expressionSet. This matrix will be \code{...$geneID} in the returned list.}
\item{showGeneLabels}{boolean.  If TRUE, and if the ranking contains gene labels, the ranking matrix will use them. This matrix will be \code{...$geneLabels} in the returned list.}
}
\value{
The method returns a list with one or two matrices: \code{...$geneLabels} and \code{...$geneID}.
}

\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's class (\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}}) help page.
}
\examples{
data(leukemiasClassifier)
getRanking(leukemiasClassifier@classificationGenes)

# Top 7 genes (two ways):
getRanking(leukemiasClassifier@genesRanking)$geneLabels[1:7,]
getRanking(getTopRanking(leukemiasClassifier@genesRanking, 7))

# Show gene ID and select a class:
getRanking(leukemiasClassifier@classificationGenes, showGeneID=TRUE
)$geneID[,"CML", drop=FALSE]
}
\keyword{methods}





