\name{getEdges-methods}
\docType{methods}
\alias{getEdges}
\alias{getEdges-methods}
\alias{getEdges,GenesNetwork-method}
\alias{edges,GenesNetwork-method}
\title{Edges in the network.}
\description{
	Returns the network's edges (relations between genes).
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
getEdges(leukemiasClassifier@genesNetwork$AML)[1:5,]
}
\keyword{methods}

