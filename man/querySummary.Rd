\name{querySummary}
\alias{querySummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of the query.
}
\description{
Counts the number of samples assigned to each class and calculates basic statistics regarding the assignment probabilities.
}
\usage{
querySummary(queryResult, showNotAssignedSamples = TRUE, numDecimals = 2, 
    verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{queryResult}{
Object returned by \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}}
}
  \item{showNotAssignedSamples}{
Logical. Shows the two most likely classes for the NotAssigned samples and the probabilities of belonging to each of them.
}
  \item{numDecimals}{
Integer. Number of decimals to show on the statistics.
}
  \item{verbose}{
Logical. If TRUE, messages indicating the execution progress will be printed on screen.
}
}
\value{
Returns a list with the following fields:
 \itemize{
\item{callRate }{ Count and percentage of assigned samples.}
\item{assigned }{ Number of samples assigned to each class and mean and SD of the assignment probabilities.}
\item{notAssignedSamples }{ Optional. Most likely classes for the Not Assigned samples.}
}
}

\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
\seealso{
Main package function and classifier training: \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
Query the classifier: \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}}
}
\examples{
##########################
## Classifier training
##########################

# Load an expressionSet:
library(leukemiasEset)
data(leukemiasEset)

# Select the train samples: 
trainSamples <- c(1:10, 13:22, 25:34, 37:46, 49:58)

# Train a classifier or load a trained one:
# leukemiasClassifier <- geNetClassifier(leukemiasEset[,trainSamples], 
#    sampleLabels="LeukemiaType", plotsName="leukemiasClassifier") 
data(leukemiasClassifier) # Sample trained classifier

##########################
## Classifier query
##########################
# Select the samples to query the classifier 
#   - Real use: samples whose class we want to known
querySamples <- "GSM330154.CEL"
#   - External validation: samples not used for training
querySamples <- c(1:60)[-trainSamples]         

# Make a query to the classifier:
queryResult <- queryGeNetClassifier(leukemiasClassifier, leukemiasEset[,querySamples])


##########################
## Query Summary
##########################
# Obtain an overview of the results
querySummary(queryResult)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif }

