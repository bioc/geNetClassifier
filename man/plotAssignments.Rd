\name{plotAssignments}
\alias{plotAssignments}
\title{
Plot assignment probabilities
}
\description{
Plots the assignment probabilities of a previous query.
}
\usage{
plotAssignments(queryResult, realLabels, 
    minProbAssignCoeff = 1, minDiffAssignCoeff = 0.8, 
	totalNumberOfClasses = NULL, pointSize=0.8)
}
\arguments{
  \item{queryResult}{
Object returned by \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}}
}
  \item{realLabels}{
Factor. Actual/real class of the samples.
}
  \item{minProbAssignCoeff}{
Numeric. See \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}} for details.
}
  \item{minDiffAssignCoeff}{
Numeric. See \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}} for details.
}
  \item{totalNumberOfClasses}{
Numeric. Total number of classes the classifier was trained with. The assignment probability is determined bassed on it. It is not needed if there are samples of all the training classes.
}
\item{pointSize}{
Numeric. Point size modifier.
}
}

\value{
Plot. If the plot is shown on the usual X11 or quartz devices, clicking on a point will identify the sample the point represents.  Press ESC or  right-click on the plot screen to finish.
}

\seealso{
Main package function and classifier training: \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
Querying the classifier: \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}} \cr
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

##########################
## Plot:
##########################
plotAssignments(queryResult, realLabels=leukemiasEset[,testSamples]$LeukemiaType)

}

\keyword{classif }
