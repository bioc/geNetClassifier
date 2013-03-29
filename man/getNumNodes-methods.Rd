\name{getNumNodes-methods}
\docType{methods}
\alias{getNumNodes}
\alias{getNumNodes-methods}
\alias{getNumNodes,GenesNetwork-method}
\alias{numNodes,GenesNetwork-method}
\title{Number of nodes in the network.}
\description{
Returns the number of nodes (genes) in the network.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "GenesNetwork")}}{
%%  ~~describe this method here~~
}
}}
\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's class (\code{\link[geNetClassifier:GenesNetwork-class]{GenesNetwork}}) help page.
}
\examples{
data(leukemiasClassifier)
getNumNodes(leukemiasClassifier@genesNetwork$AML)
}
\keyword{methods}

