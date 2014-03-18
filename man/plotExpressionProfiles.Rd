\name{plotExpressionProfiles}
\alias{plotExpressionProfiles}

\title{
Expression profiles plot.
}
\description{
Plots the expression profiles of the given genes.
}
\usage{
plotExpressionProfiles(eset, genes=NULL, fileName=NULL, geneLabels = NULL, 
    sampleLabels = NULL,  sampleColors=NULL, labelsOrder = NULL, 
    showSampleNames = FALSE, showMean = FALSE, sameScale = TRUE, 
    verbose = TRUE)

}
\arguments{
  \item{eset}{
ExpressionSet or Matrix. Gene expression of the samples.
}
  \item{genes}{
Vector or Matrix. IDs of the genes to plot. \cr
If matrix: genes should be ordered by classes. Columns should be named after the classes. If not provided, all available genes will be plot.
}
  \item{fileName}{
Character. File name to save the plots. If not provided, up to 20 genes will be shown on screen.
}
  \item{geneLabels}{
Vector or Matrix. Gene name, ID or label which should be shown in the returned results and plots.
}
  \item{sampleLabels}{
Character. PhenoData variable (column name) containing the train samples class labels. \cr
Matrix or Factor. Class labels of the train samples.
}
  \item{sampleColors}{
Character. Colors for the lines of the samples.
}
  \item{labelsOrder}{
Vector or Factor. Order in which the labels should be shown in the returned results and plots.
}
  \item{showSampleNames}{
Logical. If TRUE, the sample names are shown in the plot. Not recommended for big datasets.
}
  \item{showMean}{
Logical. If TRUE, plots the class expression mean.
}
  \item{sameScale}{
Logical. If TRUE, plots all the genes in the same expression scale.
}
sameScale
  \item{verbose}{
Logical. If TRUE, messages indicating the execution progress will be printed on screen.
}
}
\value{
The expression profiles plot, saved in the working directory as 'fileName.pdf'.
}
\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
\examples{
######
# Load libraries and expression data
######
library(geNetClassifier)
# Load an expressionSet:
library(leukemiasEset)
data(leukemiasEset)


######
# Generic expression profile plot
######
# Default settings: Receiving only an expressionSet.
# Attempts to plot all available genes 
# (If fileName is not provided, it shows up to 20 genes on screen)
plotExpressionProfiles(leukemiasEset)

# Plot expression of specific genes:
myGenes <- c("ENSG00000169575","ENSG00000078399","ENSG00000005381","ENSG00000154511")
plotExpressionProfiles(leukemiasEset, genes=myGenes, sampleLabels="LeukemiaType")

# Color samples:
plotExpressionProfiles(leukemiasEset, genes="ENSG00000078399", sampleLabels="LeukemiaType", 
sampleColors=c("red", "blue")[(leukemiasEset$LeukemiaType =="AML")+1])

######
# Expression profiles related to a classifier
######
# Train a classifier or load a trained one:
trainSamples<- c(1:10, 13:22, 25:34, 37:46, 49:58) 
# summary(leukemiasEset$LeukemiaType[trainSamples])
# leukemiasClassifier <- geNetClassifier(leukemiasEset[,trainSamples], 
#    sampleLabels="LeukemiaType", plotsName="leukemiasClassifier") 
data(leukemiasClassifier) # Sample trained classifier

# Plot expression of the train samples:
plotExpressionProfiles(leukemiasEset[,trainSamples], leukemiasClassifier, 
    sampleLabels="LeukemiaType", fileName="leukExprs.pdf")

# Plot expression of all the samples of a specific class:
classGenes <- getRanking(leukemiasClassifier@classificationGenes, 
    showGeneID=TRUE)$geneID[,"CLL"] # Feel free to modify
plotExpressionProfiles(leukemiasEset, genes=classGenes, sampleLabels="LeukemiaType", 
    fileName="cllExprs.pdf")

# Plot (on screen) the expression of the top ranked genes of each class
plotExpressionProfiles(leukemiasEset, leukemiasClassifier, sampleLabels="LeukemiaType")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif }

