\name{GeneralizationError-class}
\Rdversion{1.1}
\docType{class}
\alias{GeneralizationError}
\alias{GeneralizationError-class}
\alias{initialize,GeneralizationError-method}
%\alias{overview,GeneralizationError-method}
\alias{show,GeneralizationError-method}
\title{Class "GeneralizationError" (slot of GeNetClassifierReturn)}
\description{
Contains the estimation of the Generalization Error and the gene stats for \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} executed with the given data and parameters. \\
Calculated by 5-fold cross-validation.
}
\section{Slots}{
  \describe{
    \item{\code{accuracy}:}{\code{"Matrix"}. Accuracy and call rate.  }
    \item{\code{sensitivitySpecificity}:}{\code{"Matrix"}. Sensitivity, Specificity, Matthews Correlation Coefficient and Call Rate for each
of the classes.}
    \item{\code{confMatrix}:}{\code{"Matrix"}. Confussion matrix. }
    \item{\code{probMatrix}:}{\code{"Matrix"}. Probabilities of belonging to each class for the assigned samples. Helps identifying
where errors are likely to occur even though there were not actual errors in the cross-validation.}
    \item{\code{querySummary}:}{\code{"List"}.  Stats regarding the probability and number of assigned test samples to each class.}
    \item{\code{classificationGenes.stats}:}{\code{"List"}. Some basic statistics regarding the chosen genes. }
    \item{\code{classificationGenes.num}:}{\code{"Matrix"}. Number of genes used for each of the 5 cross-validaton classifiers.}
  }
}
\section{Methods}{
  \describe{
    \item{overview}{\code{signature(object = "GeneralizationError")}: Shows an overview of all the slots in the object. }
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
# Note: Required 'estimateGError=TRUE' 
# leukemiasClassifier <- geNetClassifier(leukemiasEset[,trainSamples], 
#    sampleLabels="LeukemiaType", plotsName="leukemiasClassifier", 
#    estimateGError=TRUE) 
data(leukemiasClassifier) # Sample trained classifier

# Global view of the returned object and its structure:
leukemiasClassifier
names(leukemiasClassifier)

#########
# Exploring the cross validation stats
# Note: Required 'estimateGError=TRUE' in geNetClassifier()
#########
# Generalization Error estimated by cross-validation:
leukemiasClassifier@generalizationError
overview(leukemiasClassifier@generalizationError)
	# i.e. probabilityMatrix:
	leukemiasClassifier@generalizationError@probMatrix
	# i.e. statistics of the genes chosen in any of the CV loops for for AML:
	leukemiasClassifier@generalizationError@classificationGenes.stats$AML
}
\keyword{classes}
