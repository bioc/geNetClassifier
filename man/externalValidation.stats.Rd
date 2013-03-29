\name{externalValidation.stats}
\alias{externalValidation.stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistics of the external validation.
}
\description{
Taking as input the confussion matrix resulting from external validation calculates the global Accuracy, Call Rate, Sensitivity, Specificity and Matthews Correlation Coefficient.
}
\usage{
externalValidation.stats(confussionMatrix, numDecimals = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{confussionMatrix}{
Confussion matrix containing the real class as rows and the assigned class as columns.
}
  \item{numDecimals}{
Integer. Number of decimals to show on the statistics.
}
}
\value{
List:
 \itemize{
 \item{global }{General classifier stats.\cr
Accuracy: Percentage of correctly assigned samples within all assigned samples. \cr
CallRate: Percentage of samples wich were assigned to a class.}
\item{byClass }{Stats by class.\cr
Sensitivity: Percentage of samples of each class which were correctly identified (Rate of true positives)\cr
Specificity: Percentage of samples assigned to a given class that really belonged to the class (Rate of true negatives)\cr
MCC (Matthews Correlation Coefficient): Measure wich takes into account both, true and false positives and negatives. (100\%: Perfect assignments)
\item{confMatrix }{Confussion matrix.}
}
}
}
\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
\seealso{
Main package function and classifier training: \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
Querying the classifier: \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}} \cr
Generating the probability matrix: \code{\link[geNetClassifier:externalValidation.probMatrix]{externalValidation.probMatrix}} \cr
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
## External Validation:
##########################
# Select the samples to query the classifier 
#   - External validation: samples not used for training
testSamples <- c(1:60)[-trainSamples]         

# Make a query to the classifier:
queryResult <- queryGeNetClassifier(leukemiasClassifier, leukemiasEset[,testSamples])

# Create the confusion matrix
confMatrix <- table(leukemiasEset[,testSamples]$LeukemiaType,queryResult$class)

# Calculate its accuracy, call rate, sensitivity and specificity:
externalValidation.stats(confMatrix)	

}
\keyword{classif }

