\name{queryGeNetClassifier}
\alias{queryGeNetClassifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Queries the classifier trained with geNetClassifier.
}
\description{
Queries the classifier trained by geNetClassifier in order to find out the class of new samples.
}
\usage{
queryGeNetClassifier(classifier, eset, minProbAssignCoeff = 1,
    minDiffAssignCoeff = 0.8, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classifier}{
Classifier returned by geNetClassifier. (@classifier)
}
  \item{eset}{
ExpressionSet or Matrix. Gene expression matrix of the new samples.
}
  \item{minProbAssignCoeff}{
Numeric. Coefficient to modify the minimum probability requird to assign a sample to a class. Reduce to improve call rate. Increase to reduce error.
0: Removes this restriction. The sample will always be assigned to the class with the highest probability.
between 0 and 1: Reduces the required probability to assign a sample to a class.
>1: Increases the required probability. Warning: if \code{minProbCoef} is equal to \code{2*number of classes}, all the samples will be left as 'NotAssigned'.
}
  \item{minDiffAssignCoeff}{
Numeric. Coefficient to modify the required difference between the two most likelly classes. Reduce to improve call rate. Increase to reduce error.
0: Removes this restriction. The probability of the second most-likely class will not be taken into account.
between 1 and 1: Reduces the required difference to assign the sample.
>1: Increases the required difference. Warning: if \code{minDiffAssignCoeff} is equal to the number of classes, all the samples will be left as 'NotAssigned'.
}
  \item{verbose}{
Logical. If TRUE, messages indicating the execution progress will be printed on screen.
}
}
\details{
By default, in order to assign a sample two contitions must be met:
\itemize{
  \item{if \code{minProbAssignCoeff = 1}}{The probability of belonging to the class should be at least double of the random probability.}
  \item{if \code{minDiffAssignCoeff = 0.8}}{The difference of probabilities between the most likely class and the second most likely class should be more than 80% the random probability. Does not apply for 2-class classifiers.}
}
This means, that in a 4-class classifier, in order to assing a sample, the highest probabiity should be at least 0.5 (2x0.25), and the next most-likely-class should have a probability at least 0.2 (80% of 0.25) lower. \cr
If these conditions are not met, the sample will be left as \code{notAssigned}.\cr
\cr
Modify the arguments values in order to modify these assignment conditions. Setting \code{minProbAssignCoeff = 0} and \code{minDiffAssignCoeff = 0} all samples will be assigned to the most likely class without any further restrictions.
}
}
\value{
List: 
\itemize{
\item{call }{Command used to execute the function.}
\item{classes }{Classes to wich each of the samples were asigned to.}
\item{probabilities }{Probabilities to the 2 classes each sample is most likelly to belong to.}
}
}
\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Main package function and classifier training: \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
Query summary: \code{\link[geNetClassifier:querySummary]{querySummary}} \cr
External validation stats: \code{\link[geNetClassifier:externalValidation.stats]{externalValidation.stats}} and \code{\link[geNetClassifier:externalValidation.probMatrix]{externalValidation.probMatrix}}
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
# There should be the same number of samples from each class.
trainSamples<- c(1:10, 13:22, 25:34, 37:46, 49:58) 
# summary(leukemiasEset$LeukemiaType[trainSamples])

# Train a classifier or load a trained one:
# leukemiasClassifier <- geNetClassifier(leukemiasEset[,trainSamples], 
#    sampleLabels="LeukemiaType", plotsName="leukemiasClassifier") 
data(leukemiasClassifier) # Sample trained classifier

##########################
## Classifier Query
##########################
# Select the samples to query the classifier 
#   - Real use: samples whose class we want to known
querySamples <- "GSM330154.CEL"
#   - External validation: samples not used for training
querySamples <- c(1:60)[-trainSamples]         

#### Make a query to the classifier ("ask" about what class the new samples are):
queryResult <- queryGeNetClassifier(leukemiasClassifier, leukemiasEset[,querySamples])

# See the class it assigned to each sample:
queryResult$class[1:5]
# Or the samples which it wasn't sure about:
t(queryResult$probabilities[,queryResult$class=="NotAssigned"])

# Obtain an overview of the results
querySummary(queryResult)

#### Optional: Modify assignment conditions
# (minDiffCoef=0, minProbCoef=0: All samples will be assigned to the most likely class)
queryResult_AssignAll <- queryGeNetClassifier(leukemiasClassifier, 
    leukemiasEset[,querySamples], minDiffAssignCoeff=0, minProbAssignCoeff=0)
# No samples are left as "NotAssigned":
queryResult$probabilities[,queryResult_AssignAll$class=="NotAssigned"]

#### External validation:
# Confusion matrix:
confMatrix <- table(leukemiasEset[,querySamples]$LeukemiaType, 
    queryResult_AssignAll$class)
# New accuracy, call rate, sensitivity and specificity:
externalValidation.stats(confMatrix)	
# Probability matrix for the assigned samples
externalValidation.probMatrix(queryResult, leukemiasEset[,querySamples]$LeukemiaType)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif }

