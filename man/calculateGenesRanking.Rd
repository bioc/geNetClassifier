\name{calculateGenesRanking}
\alias{calculateGenesRanking}

\title{
Calculate GenesRanking
}
\description{
Calculates the genes ranking and/or plots the posterior probability of the genes ordered by class ranking.
}
\usage{
calculateGenesRanking(eset=NULL, sampleLabels=NULL, 
numGenesPlot=1000, plotTitle="Significant genes", plotLp=TRUE, 
lpThreshold = 0.95, numSignificantGenesType="ranked", 
returnRanking="full", nullHiphothesisFilter=0.95,  nGenesExprDiff=1000, 
geneLabels=NULL, precalcGenesRanking=NULL, IQRfilterPercentage= 0, 
verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
ExpressionSet or Matrix. Gene expression of the train samples.
}
  \item{sampleLabels}{
Character. PhenoData variable (column name) containing the train samples class labels. \cr
Matrix or Factor. Class labels of the train samples.
}
  \item{numGenesPlot}{
Integer. Number of genes to plot.
}
  \item{plotTitle}{
Character. Plot title.
}
  \item{plotLp}{
Logical. If FALSE no plot is drawn.
}
  \item{lpThreshold}{
Numeric between 0 and 1. Required posterior probability value to consider a gene 'significant'.
}
 \item{numSignificantGenesType}{
Character. Type of count for number of genes over lpThreshold.
\itemize{
\item{"global"}. Counts all genes of a class with posterior probability over lpThreshold, even if in the final ranking they were assigned to another class.
\item{"ranked"}. Counts only genes assigned to each class.
}
}
  \item{returnRanking}{
Character. Type of ranking to return:
\itemize{
\item{"full"}. Ranking of all available genes.
\item{"lp"/"significant"/"lpThreshold"/TRUE}. Ranking of the significant genes (genes with posterior probability over lpThreshold).
\item{FALSE/NULL}. No ranking is returned.
}
}
 \item{nullHiphothesisFilter}{
Numeric between 0 and 1. Genes with a Null Hipothesis with a posterior probability over this threshold will be removed from the ranking. \cr
Null Hipothesis: They don't represent any class.
}
  \item{nGenesExprDiff}{
Numeric. Number of top genes to calculate the differencial expression for.
}
  \item{geneLabels}{
Vector or Matrix. Gene name, ID or label which should be shown in the returned results and plots.
}
  \item{IQRfilterPercentage}{
Integer. InterQuartile Range (IQR) filter applied to the initial data. Not recommended for more than two classes.
}
  \item{precalcGenesRanking}{
Allows providing a \code{genesRanking} provided by \code{geNetClassifier} or by a previous execution for the same data and parameters.
}
  \item{verbose}{
Logical. If TRUE, messages indicating the execution progress will be printed on screen.
}
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
 \itemize{
\item{\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}} }{Optional. Requested genes ranking.}
\item{Plot }{Optional. Plot of the posterior probability of the top genes.}
}
}
\seealso{\code{\link[geNetClassifier:plot.GenesRanking]{plot.GenesRanking}} is a shortcut to plotting a previusly calculated genes ranking. \cr
i.e. \code{plot(genesRanking)}}
\examples{

# Load an expressionSet:
library(leukemiasEset)
data(leukemiasEset)

# Select the train samples: 
trainSamples<- c(1:10, 13:22, 25:34, 37:46, 49:58) 
# summary(leukemiasEset$LeukemiaType[trainSamples])

\dontrun{
######
# Calculate/plot the significant genes (+ info) of a dataset 
# without training classifier/calculating network
######
# Return only significant genes ranking (default)
signGenesRanking <- calculateGenesRanking(leukemiasEset[,trainSamples], 
    sampleLabels="LeukemiaType")
numGenes(signGenesRanking)

# Return the full genes ranking:
fullRanking <- calculateGenesRanking(leukemiasEset[,trainSamples], 
    sampleLabels="LeukemiaType", returnRanking="full")
numGenes(fullRanking)
numSignificantGenes(fullRanking)
# The significant genes can then be extracted from it:
signGenesRanking2  <- getTopRanking(fullRanking, 
    numGenesClass=numSignificantGenes(fullRanking))
numGenes(signGenesRanking2)

# Changing the posterior probability required to consider genes significant:
signGenesRanking90 <- calculateGenesRanking(leukemiasEset[,trainSamples], 
    sampleLabels="LeukemiaType", lpThreshold=0.9)
numGenes(signGenesRanking90)
}
######
# Ploting previously calculated rankings:
######
# Load or calculate a ranking (or a classifier with geNetClassifier)
data(leukemiasClassifier) # Sample trained classifier, @genesRanking

# Default plot:
# - equivalent to plot(leukemiasClassifier@genesRanking)
# - in this case, the previously calculated 'fullRanking' 
#   is equivalent to 'leukemiasClassifier@genesRanking'
calculateGenesRanking(precalcGenesRanking=leukemiasClassifier@genesRanking)

# Changing arguments:
calculateGenesRanking(precalcGenesRanking=leukemiasClassifier@genesRanking, 
    numGenesPlot=5000, plotTitle="Leukemias", lpThreshold=0.9)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif }

