\name{GenesRanking-class}
\Rdversion{1.1}
\docType{class}
\alias{GenesRanking}
\alias{GenesRanking-class}
\alias{extractGenes,GenesRanking-method}
%\alias{getTopRanking,GenesRanking-method}
%\alias{gClasses,GenesRanking-method}
%\alias{genesDetails,GenesRanking-method}
%\alias{getRanking,GenesRanking-method}
\alias{initialize,GenesRanking-method}
%\alias{nGenes,GenesRanking-method}
%\alias{numSignificantGenes,GenesRanking-method}
%\alias{overview,GenesRanking-method}
%\alias{plot,GenesRanking,missing-method}
%\alias{setProperties,GenesRanking-method}
\alias{show,GenesRanking-method}

\title{Class "GenesRanking"}
\description{
Contains a genes ranking and the genes info calculated by \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}.\cr
(Slots \code{@classificationGenes} and \code{@genesRanking} from \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} output)
}
\section{Methods}{
  \describe{
    \item{genesDetails}{\code{signature(object = "GenesRanking")}: Returns data.frames with information about the genes. }
    \item{getRanking}{\code{signature(object = "GenesRanking")}: Returns a matrix containing the ranked genes. }
    \item{getTopRanking}{\code{signature(object = "GenesRanking", numGenesClass)}: Returns a new GenesRanking object containing only the top genes of each class. }
    \item{gClasses}{\code{signature(object = "GenesRanking")}: Returns the classes for which the genes are ranked.}
    \item{numGenes}{\code{signature(object = "GenesRanking")}: Returns the number of available ranked genes per class. }
    \item{numSignificantGenes}{\code{signature(object = "GenesRanking")}: Returns the number of significant genes per class (genes over the given posterior probability threshold).}
    \item{plot}{\code{signature(x = "GenesRanking", y = "missing")}: Plots the genes' posterior probability. Wrapper of \code{\link[geNetClassifier:calculateGenesRanking]{calculateGenesRanking}}.}
 }
}

\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
\seealso{

For more information on how the ranking is calculated and how to interpret the given information, see the package vignette.\cr

Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}\cr
Plot the ranking genes's posterior probability: \code{\link[geNetClassifier:plot.GenesRanking]{plot.GenesRanking}}
}
\examples{
######
# Calculate a genesRanking
######

# Load an expressionSet:
library(leukemiasEset)
data(leukemiasEset)

# Select the train samples: 
trainSamples<- c(1:10, 13:22, 25:34, 37:46, 49:58) 
# summary(leukemiasEset$LeukemiaType[trainSamples])

# Calculate the genesRanking with calculateGenesRanking()
\dontrun{
genesRanking <- calculateGenesRanking(leukemiasEset[,trainSamples], 
    sampleLabels="LeukemiaType", returnRanking="full")}

# geNetClassifier() also calculates a genes ranking
# Sample output: 
data(leukemiasClassifier) 
genesRanking <- leukemiasClassifier@genesRanking

######
# Exploring the rankings
######
# Number of available genes in the ranking:
numGenes(genesRanking)

# Number of significant genes (genes with posterior probability over the threshold. 
# Default: lpThreshold=0.95):
numSignificantGenes(genesRanking)

# Top 10 genes of CML:
genesDetails(genesRanking)$CML[1:10,]

# To get a sub ranking with the top 10 genes:
getTopRanking(genesRanking, 10)

# Genes details of the top 10 genes:
genesDetails(getTopRanking(genesRanking, 10))	

######
# Exploring the genes used for training the classifier
######
numGenes(leukemiasClassifier@classificationGenes)
leukemiasClassifier@classificationGenes
#genesDetails(leukemiasClassifier@classificationGenes)  # List by classes
genesDetails(leukemiasClassifier@classificationGenes)$AML # Show a class genes
# If your R console wraps the table rows, try widening your display width: 
# options(width=200)

######
# Creating a GenesRanking object
# i.e. To use geNetClassifier() with a ranking based on another algorithm
######

### 1. Calculate gene scores 
# Two classes:
geneScore <- matrix(sample(seq(0,1,by=0.01), size=100, replace=TRUE))
colnames(geneScore) <- "BothClasses"
rownames(geneScore) <- paste("Gene", 1:100, sep="")

# More than two classes:
geneScore <- matrix(sample(seq(0,1,by=0.01), size=300, replace=TRUE), ncol=3)
colnames(geneScore) <- paste("Class", 1:3, sep="")
rownames(geneScore) <- paste("Gene", 1:100, sep="")

### 2. Create object
postProb <- geneScore
ord <- apply(postProb, 2, function(x) order(x, decreasing=TRUE))
numGenesClass <- apply(postProb, 2, function(x) sum(!is.na(x)))
customRanking <- new("GenesRanking", postProb=postProb, ord=ord, numGenesClass=numGenesClass)

# GenesRanking object ready:
customRanking
genesDetails(customRanking)
customRanking@numGenesClass
numSignificantGenes(customRanking)

# geNetClassifier(..., precalcGenesRanking = customRanking)

}
\keyword{classes}
