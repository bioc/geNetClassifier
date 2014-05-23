\name{GenesNetwork-class}
\Rdversion{1.1}
\docType{class}
\alias{GenesNetwork}
\alias{GenesNetwork-class}
%\alias{edges,GenesNetwork-method}
%\alias{getSubNetwork,GenesNetwork-method}
%\alias{getEdges,GenesNetwork-method}
% \alias{getNodes,GenesNetwork-method}
% \alias{getNumEdges,GenesNetwork-method}
% \alias{getNumNodes,GenesNetwork-method}
\alias{initialize,GenesNetwork-method}
%\alias{nodes,GenesNetwork-method}
%\alias{numEdges,GenesNetwork-method}
%\alias{numNodes,GenesNetwork-method}
% \alias{overview,GenesNetwork-method}
\alias{show,GenesNetwork-method}

\title{Class "GenesNetwork"}
\description{
Contains the network returned by \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}. (Slot: \code{@genesNetwork})
}
\section{Methods}{
\describe{
\item{getNodes}{\code{signature(object = "GenesNetwork")}: Returns the network nodes (genes). }
\item{getEdges}{\code{signature(object = "GenesNetwork")}: Returns the network edges (relationships). }
\item{getNumNodes}{\code{signature(object = "GenesNetwork")}: Returns the number of nodes (genes) in the network.}
\item{getNumEdges}{\code{signature(object = "GenesNetwork")}: Returns the number of edges (relationships) in the network, }
\item{getSubNetwork}{\code{signature(network = "GenesNetwork")}: Returns a new network containing only the given genes. }
\item{network2txt}{\code{signature(network = "GenesNetwork")}: Exports the network as text file. }
}
}
\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
\seealso{
Main package function and classifier training: \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}
Plot network or export as iGraph: \code{\link[geNetClassifier:plotNetwork]{plotNetwork}}

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
names(leukemiasClassifier)

# List of Networks by classes:
leukemiasClassifier@genesNetwork
# Access to the nodes or edges of each network:
getEdges(leukemiasClassifier@genesNetwork$AML)[1:5,]
getNodes(leukemiasClassifier@genesNetwork$AML)[1:50]
		
	
######
# Plotting
######
# Example: Plotting the sub-network of a class classificationGenes
# Get the sub-network containing only the classification genes:
subNet <- getSubNetwork(leukemiasClassifier@genesNetwork, 
    leukemiasClassifier@classificationGenes)
# Get the classification genes' info/details:
clGenesInfo <- genesDetails(leukemiasClassifier@classificationGenes)

# Plot the network of the class "ALL"
plotNetwork(subNet$ALL, genesInfo=clGenesInfo, 
    plotOnlyConnectedNodesNetwork=FALSE)	
}
\keyword{classes}
