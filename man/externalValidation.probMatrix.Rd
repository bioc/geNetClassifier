\name{externalValidation.probMatrix}
\alias{externalValidation.probMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability matrix.
}
\description{
Generates the probability matrix.
}
\usage{
externalValidation.probMatrix(queryResult, realLabels, numDecimals=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{queryResult}{
Object returned by \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}} 
}
  \item{realLabels}{
Factor. Actual/real class of the samples.
}
  \item{numDecimals}{
Integer. Number of decimals to return.
}
}
\details{
A probability matrix contains the probabilities of assigning each assigned sample to each class. 
They help identifying where errors are likelly to occur even though there were not actual errors in the external/cross validation.
}
\value{
The probability matrix.
}
\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
\seealso{
Main package function and classifier training: \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
Query the classifier: \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}} \cr
Query summary: \code{\link[geNetClassifier:querySummary]{querySummary}} \cr
External validation stats: \code{\link[geNetClassifier:externalValidation.stats]{externalValidation.stats}}
}
\examples{
##########################
## Classifier training
##########################
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

##########################
## External Validation
##########################
# Select the samples to query the classifier 
#   - External validation: samples not used for training
testSamples <- c(1:60)[-trainSamples]         

# Make a query to the classifier:
queryResult <- queryGeNetClassifier(leukemiasClassifier, leukemiasEset[,testSamples])

# Obtain the probability matrix for the assigned samples:
externalValidation.probMatrix(queryResult, leukemiasEset[,testSamples]$LeukemiaType)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif }

