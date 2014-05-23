\name{plotNetwork}
\alias{plotNetwork}
\alias{plotNetwork,GenesNetwork-method}
\alias{plot,GenesNetwork-method}
\alias{plot.GenesNetwork}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot GenesNetwork
}
\description{
Plots the coexpression and/or mutual information network for the given genes.
}
\usage{
plotNetwork(genesNetwork, classificationGenes=NULL, genesRanking=NULL, 
genesInfo=NULL, geneLabels=NULL, returniGraphs=FALSE, 
plotType="dynamic", fileName=NULL, 
plotAllNodesNetwork=TRUE, plotOnlyConnectedNodesNetwork=FALSE,  
plotClassifcationGenesNetwork=FALSE, 
labelSize=0.5, vertexSize=NULL, width=NULL, height=NULL, verbose=TRUE)
}
\arguments{
 \item{genesNetwork}{
List of GenesNetwork returned by \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}. (\code{@genesNetwork})
}
  \item{classificationGenes}{
Matrix or classificationGenes returned by \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}. (\code{@classificationGenes})
}
  \item{genesRanking}{
Matrix or genesRanking returned by \code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}. (\code{@genesRanking})
}
  \item{genesInfo}{
List or data.frame with the properties of the genes to plot: \code{genesDetails(_@genesRanking)}
}
  \item{geneLabels}{
Vector or Matrix. Gene name, ID or label which should be shown in the returned results and plots.
}
  \item{returniGraphs}{
deprecated. A list with the plotted networks as igraph objects is always returned (see \code{\link[base:invisible]{invisible}}), assign it to a variable if needed.
}
  \item{plotType}{
Character. "dynamic": Interactive plot. "static": One canvas split for the different networks. "pdf": All the networks are saved into a pdf file.
}
  \item{fileName}{
Character. File name to save the plot with. If not provided, the plots will be shown through the standard output device.
}
  \item{plotAllNodesNetwork}{
Logical. If TRUE, plots a network only with all the available genes
}
  \item{plotOnlyConnectedNodesNetwork}{
Logical. If TRUE, plots a network only with the connected nodes/genes
}
  \item{plotClassifcationGenesNetwork}{
Logical. If TRUE, plots a network only with the classification genes
}
  \item{labelSize}{
Integer. Gene/node label size for static and pdf plots.
}
  \item{vertexSize}{
Integer. Vertex minimum size.
}
  \item{width}{
Numeric. Dinamic or pdf plot width.
}
  \item{height}{
Numeric. Dinamic or pdf plot height.
}
 \item{verbose}{
Logical. If TRUE, messages indicating the execution progress will be shown. 
} 
}
\value{
  \item{Graph list}{
	List with the plotted \code{\link[igraph:igraph]{igraph}} objects.
	}
  \item{Network plots}{
	Shown throught the standard output devide or saved in the working directory as 'fileName.pdf' if \code{fileName} was provided.
	}
}
\references{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}}\cr

Package \code{\link[igraph:igraph]{igraph}}
}
\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
\seealso{
\code{plot.GenesNetwork()} is an alias to this function. It can allso be called as i.e. \code{plotNetwork(clGenesSubNet$ALL)}
Note: The slot @genesNetwork returned by geNetClassifier is a List of GenesNetworks!
}
\examples{

data(leukemiasClassifier)

# Step 1: Select a network or sub network
# Sub-network containing only the classification genes:
clGenesSubNet <- getSubNetwork(leukemiasClassifier@genesNetwork, 
    leukemiasClassifier@classificationGenes)
# Step 2: Select the details/info about the genes to plot
# Classification genes' info:
clGenesInfo <- genesDetails(leukemiasClassifier@classificationGenes)

# Step 3: Plot the network
# Network plots can be interactive or plotted as PDF file.
#  - - Use plotType="pdf" to save the network as a static pdf file. 
#       This option is recommended for getting an overview of several networks.
#  - - To get an interactive network, just skip this argument. 

# Plot ALL network:
plotNetwork(clGenesSubNet$ALL, genesInfo=clGenesInfo)

# Plot AML network containing only the conected nodes:
plotNetwork(clGenesSubNet$ALL, genesInfo=clGenesInfo, 
 plotAllNodesNetwork=FALSE, plotOnlyConnectedNodesNetwork=TRUE)

# The equivalent code to the plot geNetClassifier creates by default is:
topRanking <- getTopRanking(leukemiasClassifier@genesRanking, numGenesClass=100)
netTopGenes <- getSubNetwork(leukemiasClassifier@genesNetwork, 
 getRanking(topRanking, showGeneID=TRUE)$geneID)
plotNetwork(netTopGenes,  classificationGenes=leukemiasClassifier@classificationGenes, 
 genesRanking=topRanking, plotAllNodesNetwork=TRUE, 
 plotOnlyConnectedNodesNetwork=TRUE, plotType="pdf", 
 labelSize=0.3, fileName="leukemiasClassifier")

# In order to save the network as text file, you can use:
network2txt(leukemiasClassifier@genesNetwork, filePrefix="leukemiasNetwork")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif }

