\name{overview-methods}
\docType{methods}
\alias{overview}
\alias{overview-methods}
\alias{overview,GenesNetwork-method}
\alias{overview,GenesRanking-method}
\alias{overview,GeNetClassifierReturn-method}
\alias{overview,GeneralizationError-method}
\title{Overview}
\description{
Provides an overview of all the slots in the object.
}
\section{Methods}{It can be applied to the following classes:
\describe{
\item{\code{signature(object = "GenesNetwork")}}{}
\item{\code{signature(object = "GenesRanking")}}{}
\item{\code{signature(object = "GeNetClassifierReturn")}}{}
\item{\code{signature(object = "GeneralizationError")}}{}
}}
\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's classes help pages:\cr
\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}} \cr
\code{\link[geNetClassifier:GenesNetwork-class]{GenesNetwork}} \cr
\code{\link[geNetClassifier:GeNetClassifierReturn-class]{GeNetClassifierReturn}} \cr
\code{\link[geNetClassifier:GeneralizationError-class]{GeneralizationError}} \cr
}
\examples{
library(geNetClassifier)
data(leukemiasClassifier)
# geNetClassifier return:
overview(leukemiasClassifier)
# Generalization Error and stats estimated by cross-validation:
overview(leukemiasClassifier@generalizationError)
# A GenesNetwork:
# (a class has to be selected, otherwise it is a list)
overview(leukemiasClassifier@genesNetwork$ALL)

# For a GenesRanking,  we recommend to use genesDetails() instead:
genesDetails(leukemiasClassifier@classificationGenes)$AML
}
\keyword{methods}



