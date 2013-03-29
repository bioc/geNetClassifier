\name{getNodes-methods}
\docType{methods}
\alias{getNodes}
\alias{getNodes-methods}
\alias{getNodes,GenesNetwork-method}
\alias{nodes,GenesNetwork-method}
\title{Nodes in the network.}
\description{
Returns the network's nodes (genes).
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
getNodes(leukemiasClassifier@genesNetwork$AML)[1:5]
}
\keyword{methods}

