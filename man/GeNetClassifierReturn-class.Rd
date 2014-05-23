\name{GeNetClassifierReturn-class}
\Rdversion{1.1}
\docType{class}
\alias{GeNetClassifierReturn}
\alias{GeNetClassifierReturn-class}
\alias{names,GeNetClassifierReturn-method}
\alias{getSlots,GeNetClassifierReturn-method}
\alias{initialize,GeNetClassifierReturn-method}
%\alias{overview,GeNetClassifierReturn-method}
\alias{show,GeNetClassifierReturn-method}

\title{Class "GeNetClassifierReturn"}
\description{
Object wich wraps all the items returned by \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}. It usually contains the classifier, the genes ranking and information, the network and any other requested statistics.
}
\section{Methods}{
\describe{
\item{names}{\code{signature(x = "GeNetClassifierReturn")}: Shows the available slots in the object.  }
\item{overview}{\code{signature(object = "GeNetClassifierReturn")}: Shows an overview of all the slots in the object.  }
}
}

\section{Slots}{
Available slots deppends on the arguments used to call \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier()}}: 
  \describe{
 \item{\code{call}:}{
call. Always available.
\cr
}
 \item{\code{classifier}:}{
list.  SVM classifier. Only available if geNetClassifier() was called with option \code{buildClassifier=TRUE} (default settings).
}
 \item{\code{classificationGenes}:}{
\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}}. Genes used to train the classifier. Only available if geNetClassifier() was called with option \code{buildClassifier=TRUE} (default settings).\cr
}

 \item{\code{generalizationError}:}{
\code{\link[geNetClassifier:GeneralizationError-class]{GeneralizationError}}. Statistics calculated for the current training set and options. \cr
Only available if geNetClassifier() was called with option \code{estimateGError=TRUE} (False by default).
}

 \item{\code{genesRanking}:}{
\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}}. Whole genes ranking (if  \code{returnAllGenesRanking=TRUE}) or significant genes ranking (if  \code{returnAllGenesRanking=FALSE}, includes only the genes with posterior probability over \code{lpThreshold})
}
 \item{\code{genesRankingType}:}{
character. "all", "significant" or "significantNonRedundant" \cr
}

 \item{\code{genesNetwork}:}{
List of \code{\link[geNetClassifier:GenesNetwork-class]{GenesNetwork}}. Only available if geNetClassifier() was called with option \code{calculateNetwork=TRUE} (default settings). 
}
 \item{\code{genesNetworkType}:}{
character. At the moment, only "topGenes" available.\cr
}
}
}
\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
\seealso{
Main package function and classifier training: \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}\cr
\code{\link[geNetClassifier:plot.GeNetClassifierReturn]{plot.GeNetClassifierReturn}}
}
\examples{
######
# Load data and train a classifier
######

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
# Explore the returned object
######
# Global view of the object and its structure:
leukemiasClassifier
names(leukemiasClassifier)

### Depending on the available slots:
# Call and acess to the classifier:
leukemiasClassifier@call
leukemiasClassifier@classifier

# Genes used for training the classifier:
numGenes(leukemiasClassifier@classificationGenes)
leukemiasClassifier@classificationGenes
# Show de tetails of the genes of a class
genesDetails(leukemiasClassifier@classificationGenes)$AML 
# If your R console wraps the table rows, try widening your display width: 
# options(width=200)

# Generalization Error estimated by cross-validation:
leukemiasClassifier@generalizationError
overview(leukemiasClassifier@generalizationError)
# i.e. probabilityMatrix:
leukemiasClassifier@generalizationError@probMatrix
# i.e. statistics of the genes chosen in any of the CV loops for for AML:
leukemiasClassifier@generalizationError@classificationGenes.stats$AML
	
# List of Networks by classes:
leukemiasClassifier@genesNetwork
# Access to the nodes or edges of each network:
	getEdges(leukemiasClassifier@genesNetwork$AML)
	getNodes(leukemiasClassifier@genesNetwork$AML)	
		
# Genes ranking:
leukemiasClassifier@genesRanking
	# Number of available genes in the ranking:
	numGenes(leukemiasClassifier@genesRanking)
	# Number of significant genes 
	# (genes with posterior probability over lpThreshold, default=0.75)
	numSignificantGenes(leukemiasClassifier@genesRanking)		
	# Top 10 genes of CML:
	genesDetails(leukemiasClassifier@genesRanking)$CML[1:10,]
	# To get a sub ranking with the top 10 genes:
	getTopRanking(leukemiasClassifier@genesRanking, 10)
	# Genes details of the top 10 genes:
	genesDetails(getTopRanking(leukemiasClassifier@genesRanking, 10))
}
\keyword{classes}

