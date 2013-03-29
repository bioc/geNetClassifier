\name{leukemiasClassifier}
\docType{data}
\alias{leukemiasClassifier}
\title{Sample leukemias classifier}
\description{
A sample of the object returned by geNetClassifier. Containins the classifier, the network, and the gene statistics.
}
\usage{
data(leukemiasClassifier)
}
\format{
 \code{\link[geNetClassifier:GeNetClassifierReturn-class]{GeNetClassifierReturn} object} \cr
}

\examples{
data(leukemiasClassifier)
# Global view of the object and its structure:
leukemiasClassifier
names(leukemiasClassifier)

# Call and acess to the classifier:
leukemiasClassifier@call
leukemiasClassifier@classifier

# Genes used for training the classifier:
numGenes(leukemiasClassifier@classificationGenes)
leukemiasClassifier@classificationGenes
genesDetails(leukemiasClassifier@classificationGenes)
		
# Generalization Error estimated by cross-validation:
# 	leukemiasClassifier@generalizationError
#		overview(leukemiasClassifier@generalizationError)
	
# List of Networks by classes:
	leukemiasClassifier@genesNetwork

# Access to the nodes or edges of each network:
	getEdges(leukemiasClassifier@genesNetwork$AML)[1:5,]
	getNodes(leukemiasClassifier@genesNetwork$AML)[1:50]	
		
# Global genes ranking:
	leukemiasClassifier@genesRanking
		numGenes(leukemiasClassifier@genesRanking)
		numSignificantGenes(leukemiasClassifier@genesRanking)
#		getTopRanking(leukemiasClassifier@genesRanking, 10)
}
\keyword{classif}
\keyword{leukemia}



