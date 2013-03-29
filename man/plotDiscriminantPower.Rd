\name{plotDiscriminantPower}
\alias{plotDiscriminantPower}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the genes' Discriminant Power.
}
\description{
Calculates and plots the Discriminant Power of the genes in the given classifier. 
}
\usage{
plotDiscriminantPower(classifier, classificationGenes = NULL, 
    geneLabels = NULL, classNames = NULL, plotDP = TRUE, 
	fileName = NULL, returnTable = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classifier}{
Classifier returned by \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}. (@classifier)
}
  \item{classificationGenes}{
Vector or Matrix. IDs of the genes to plot. \cr
If matrix: genes should be ordered by classes. Columns should be named after the classes.
}
  \item{geneLabels}{
Vector or Matrix. Gene name, ID or label which should be shown in the returned results and plots.
}
  \item{classNames}{
Named vector. Short version of the class names if different from the ones used to train the classifier.
}
  \item{plotDP}{
Logical. If TRUE, plots the discriminant power of the given genes.
}
  \item{fileName}{
Character. File name to save the plot with. If not provided, the plots will be shown through the standard output device.
}
  \item{returnTable}{
Logical. If TRUE, returns a table with the genes discriminant power.
}
  \item{verbose}{
Logical. If TRUE, messages indicating the execution progress will be printed on screen.
}
}
\details{
The Discriminant Power represents the weight the (SVM) classifier gives each gene to separate the classes. It is calculated based on the coordinates of the support vectors.
Genes with a high Discriminant Power are better for identifying samples from the class.
}
\value{
 \itemize{
\item{Data frame} {Optional. Data.frame containing the genes and their Discriminant Power.}
\item{Discriminant Power plot} {Optional. Shown throught the standard output devide or saved in the working directory as 'fileName.pdf' if \code{fileName} was provided.}
}
}
\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}

\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}
}
\examples{
######
# Load data and train a classifier
######
library(geNetClassifier)
# Load an expressionSet:
library(leukemiasEset)
data(leukemiasEset)

# Select the train samples: 
trainSamples<- c(1:10, 13:22, 25:34, 37:46, 49:58) 
# summary(leukemiasEset$LeukemiaType[trainSamples])

# Train a classifier or load a trained one:
# leukemiasClassifier <- geNetClassifier(leukemiasEset[,trainSamples], 
#    sampleLabels="LeukemiaType", plotsName="leukemiasClassifier") 
data(leukemiasClassifier) # Sample trained classifier

######
# Discriminant Power
######
# Default (plots up to 20 genes)
plotDiscriminantPower(leukemiasClassifier)
# Plot a specific gene:
plotDiscriminantPower(leukemiasClassifier, classificationGenes="ENSG00000169575")
# Plot top5 genes of a class, and return their discriminant power:
# Note: The discriminant Power can only be calculated for 'classificationGenes' 
#            (genes chosen for training the classifier)
genes <- getRanking(leukemiasClassifier@classificationGenes, 
    showGeneID=TRUE)$geneID[1:5,"AML",drop=FALSE] # Top 5 genes of AML
discPowerTable2 <- plotDiscriminantPower(leukemiasClassifier, 
    classificationGenes=genes, returnTable=TRUE)

# For plotting more than 20 genes or saving the plots as .pdf, provide a fileName
plotDiscriminantPower(leukemiasClassifier, 
     fileName="leukemiasClassifier_DiscriminantPower.pdf")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif }

