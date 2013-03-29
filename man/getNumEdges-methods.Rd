\name{getNumEdges-methods}
\docType{methods}
\alias{getNumEdges}
\alias{getNumEdges-methods}
\alias{getNumEdges,GenesNetwork-method}
\alias{numEdges,GenesNetwork-method}
\title{Number of edges in the network.}
\description{
Returns the number of edges (relationships) in the network.
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
getNumEdges(leukemiasClassifier@genesNetwork$AML)
}
\keyword{methods}

