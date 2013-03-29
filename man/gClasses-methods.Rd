\name{gClasses-methods}
\docType{methods}
\alias{gClasses}
\alias{gClasses-methods}
\alias{gClasses,GenesRanking-method}
\title{Classes in the ranking.}
\description{
Returns the names of the classes in a GenesRanking
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
gClasses(leukemiasClassifier@genesRanking)
}
\keyword{methods}

