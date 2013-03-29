\name{plot.GenesRanking}
\alias{plotGenesRanking}
\alias{plot.GenesRanking}
\alias{plot,GenesRanking-method}

\title{
Plot GenesRanking
}
\description{
Plots the posterior probability of the genes ordered by class ranking.
}
\usage{
\method{plot}{GenesRanking}(x, y="missing", numGenesPlot=1000, 
    plotTitle="Significant genes", lpThreshold = 0.75, ...)
}
\arguments{
  \item{x}{GenesRanking.
}
  \item{numGenesPlot}{Numeric. Number of genes to plot.
}
  \item{plotTitle}{Character. Plot main title.
}
  \item{lpThreshold}{ Numeric between 0 and 1. Required posterior probability value to consider a gene 'significant'.
}
  \item{y}{Not required.}
  \item{...}{Not required}
}
\details{
Significant genes: Genes with posterior probability over 'lpThreshold'.\cr
More significant genes may mean:
 \itemize{
\item{}{Very different class}
\item{}{More systemic disease}
}
Plot lines represet the posterior probability of genes, sorted by rank from left to right.\cr\cr

In order to find genes that diferentiate the classes from each other, the function ranks the genes bassed on their posterior probability for each class. \cr
The posterior probability represents how well a gene differentiates samples from a class, from samples from other classes. Therefore, Genes with high posterior probability are good to differentiate a class from all the others.\cr
This posterior probability is calculated by \code{\link[EBarrays:emfit]{emfit (pkg:EBarrays)}}, an expectation-maximization (EM) algorithm for gene expression mixture model.\cr
}
\value{
Posterior probability plot of the top genes.
}
\examples{
# Load or calculate a ranking (or a classifier with geNetClassifier)
data(leukemiasClassifier) # Sample trained classifier, @genesRanking

# Default plot:
plot(leukemiasClassifier@genesRanking)

# Changing options:
plot(leukemiasClassifier@genesRanking, 
    numGenesPlot=5000, plotTitle="Leukemias", lpThreshold=0.9)
}

\keyword{classif }


