# Usage example for geNetClassifier
# R code from the vignette

# Installation
source("http://bioconductor.org/biocLite.R")
biocLite("geNetClassifier")
biocLite("leukemiasEset")

# Loading data:
library(leukemiasEset)
data(leukemiasEset)
summary(leukemiasEset$LeukemiaType) 
?leukemiasEset
# Any other data set could be loaded using the function load: 
load("...")

# Alternative gene labels and annotation (optional):
# (Download from http://bioinfow.dep.usal.es/xgate/mapping/mapping.php?content=annotationfiles)
geneSymbols <- eval(as.name(load("genes-human-annotation.R")[1]))[,"gene_symbol",drop=F]	
leukEset_protCoding <- leukemiasEset[featureNames(leukemiasEset) %in% rownames(genes.human.Annotation[genes.human.Annotation$biotype %in% "protein_coding",]),]

# Loading package and getting help:
library(geNetClassifier)
vignette("geNetClassifier-vignette")
?geNetClassifier

###############################################################
####               Executing geNetClassifier()
###############################################################

# Select train samples
trainSamples <- c(1:10, 13:22, 25:34, 37:46, 49:58)
summary(leukemiasEset$LeukemiaType[trainSamples])

# Execution: 
leukemiasClassifier <- geNetClassifier(eset=leukemiasEset[,trainSamples], sampleLabels="LeukemiaType", plotsName="leukemiasClassifier") # Standard
leukemiasClassifier <- geNetClassifier(eset=leukemiasEset[,trainSamples], sampleLabels="LeukemiaType", plotsName="leukemiasClassifier", estimateGError=TRUE) # With generalization error
leukemiasClassifier <- geNetClassifier(eset=leukemiasEset[,trainSamples], sampleLabels="LeukemiaType", plotsName="leukemiasClassifier", skipInteractions=TRUE, maxGenesTrain=20) # Fast

# Save / load results
save(leukemiasClassifier, file="leukemiasClassifier.RData") 
load("leukemiasClassifier.RData")	

# Load built in example :
data(leukemiasClassifier)  					
# Code used for creating this example:
leukemiasClassifier <- geNetClassifier(eset=leukemiasEset[,trainSamples], sampleLabels="LeukemiaType", plotsName="leukemiasClassifier", estimateGError=TRUE, geneLabels=geneSymbols)

# Explore resuts:
leukemiasClassifier
class(leukemiasClassifier)
names(leukemiasClassifier)

leukemiasClassifier@call

# Ranking
leukemiasClassifier@genesRanking
numGenes(leukemiasClassifier@genesRanking)
getRanking(leukemiasClassifier@genesRanking, showGeneID=TRUE)$geneID[1:5,]

miniRanking <- getTopRanking(leukemiasClassifier@genesRanking, 7)
getRanking(miniRanking)
options(width=200)
genesDetails(miniRanking)$AML

numSignificantGenes(leukemiasClassifier@genesRanking)
plot(leukemiasClassifier@genesRanking)
plot(leukemiasClassifier@genesRanking, numGenesPlot=3000, lpThreshold=0.80)

# Classifier
leukemiasClassifier@classifier

leukemiasClassifier@classificationGenes
numGenes(leukemiasClassifier@classificationGenes)
genesDetails(leukemiasClassifier@classificationGenes)$ALL

leukemiasClassifier@generalizationError
overview(leukemiasClassifier@generalizationError)
leukemiasClassifier@generalizationError@confMatrix
leukemiasClassifier@generalizationError@probMatrix
leukemiasClassifier@generalizationError@classificationGenes.stats$CLL
leukemiasClassifier@generalizationError@classificationGenes.num

# Network
leukemiasClassifier@genesNetwork
overview(leukemiasClassifier@genesNetwork$AML)
getNumEdges(leukemiasClassifier@genesNetwork$AML)
getNumNodes(leukemiasClassifier@genesNetwork$AML)
getEdges(leukemiasClassifier@genesNetwork$AML)[1:4,]
getNodes(leukemiasClassifier@genesNetwork$AML)[1:12]

network2txt(leukemiasClassifier@genesNetwork, filePrefix="leukemiasNetwork")
geneNtwsInfo <- lapply(leukemiasClassifier@genesNetwork, function(x) write.table(getEdges(x), file=paste("leukemiaNtw_",getEdges(x)[1,"class1"],".txt",sep="")))

###############################################################
####               External Validation
###############################################################

# Select samples and Query
testSamples <- c(1:60)[-trainSamples]
queryResult <- queryGeNetClassifier(leukemiasClassifier, leukemiasEset[,testSamples])
queryResult$class
queryResult$probabilities

# Explore/evaluate results
confusionMatrix <- table(leukemiasEset[,testSamples]$LeukemiaType, queryResult$class)
confusionMatrix
externalValidation.stats(confusionMatrix)
externalValidation.probMatrix(queryResult, leukemiasEset[,testSamples]$LeukemiaType, numDecimals=3)

# Assignment Conditions
queryResult_AssignAll <- queryGeNetClassifier(leukemiasClassifier, leukemiasEset[,testSamples], minProbAssignCoeff=0, minDiffAssignCoeff=0)
queryResult_AssignAll$class
which(queryResult_AssignAll$class=="NotAssigned")

queryResult_AssignLess <- queryGeNetClassifier(leukemiasClassifier, leukemiasEset[,testSamples], minProbAssignCoeff=1.5, minDiffAssignCoeff=1)
queryResult_AssignLess$class
queryResult_AssignLess$probabilities[,queryResult_AssignLess$class=="NotAssigned", drop=FALSE]
confusionMatrix2 <- table(leukemiasEset[,testSamples]$LeukemiaType, queryResult_AssignLess$class)
confusionMatrix2
externalValidation.stats(confusionMatrix2)

###############################################################
####               Classifier Query
###############################################################

# Select samples and Query
testSamples <- c(1:60)[-trainSamples]
queryResult_AsUnkown <- queryGeNetClassifier(leukemiasClassifier, leukemiasEset[,testSamples])

# Explore results
names(queryResult_AsUnkown)
queryResult_AsUnkown$class
t(queryResult_AsUnkown$probabilities[ , queryResult$class=="NotAssigned"])
querySummary(queryResult_AsUnkown, numDecimals=3)


###############################################################
####               Plots
###############################################################

# Gene's ranking
plot(leukemiasClassifier@genesRanking)
plot(leukemiasClassifier@genesRanking, numGenesPlot=3000, plotTitle="5 classes: ALL, AML, CLL, CML, NoL", lpThreshold=0.80)
ranking <- calculateGenesRanking(leukemiasEset[,trainSamples], "LeukemiaType")

# Expression profiles
myGenes <- c("ENSG00000169575", "ENSG00000078399", "ENSG00000176890", "ENSG00000121742")
plotExpressionProfiles(leukemiasEset, myGenes, sampleLabels="LeukemiaType")

plotExpressionProfiles(leukemiasEset[,trainSamples], leukemiasClassifier, sampleLabels="LeukemiaType", fileName="leukExprs_trainSamples.pdf")

classGenes <- getRanking(leukemiasClassifier@classificationGenes, showGeneID=TRUE)$geneID[,"AML"]
plotExpressionProfiles(leukemiasEset, genes=classGenes, sampleLabels="LeukemiaType", fileName="AML.pdf", geneLabels=geneSymbols)


# Discriminant Power
plotDiscriminantPower(leukemiasClassifier, classificationGenes="ENSG00000169575")
discPowerTable <- plotDiscriminantPower(leukemiasClassifier, classificationGenes=getRanking(leukemiasClassifier@classificationGenes, showGeneID=T)$geneID[1:4,"AML",drop=FALSE], returnTable=TRUE)

# Network
clGenesSubNet <- getSubNetwork(leukemiasClassifier@genesNetwork, leukemiasClassifier@classificationGenes)
clGenesInfo <- genesDetails(leukemiasClassifier@classificationGenes)
plotNetwork(genesNetwork=clGenesSubNet$ALL, genesInfo=clGenesInfo, plotAllNodesNetwork=FALSE, plotOnlyConnectedNodesNetwork=TRUE)
plotNetwork(genesNetwork=clGenesSubNet$ALL, genesInfo=clGenesInfo)

top30g <- getRanking(leukemiasClassifier@genesRanking, showGeneID=TRUE)$geneID[1:30,]
top30gSubNet <- getSubNetwork(leukemiasClassifier@genesNetwork, top30g)
top30gInfo <- lapply(genesDetails(leukemiasClassifier@genesRanking), function(x) x[1:30,]) 
plotNetwork(genesNetwork=top30gSubNet$AML, genesInfo=top30gInfo$AML)

top100gRanking <- getTopRanking(leukemiasClassifier@genesRanking, numGenes=100)
top100gSubNet <- getSubNetwork(leukemiasClassifier@genesNetwork, getRanking(top100gRanking, showGeneID=TRUE)$geneID)
plotNetwork(genesNetwork=top100gSubNet, classificationGenes=leukemiasClassifier@classificationGenes, genesRanking=top100gRanking, plotAllNodesNetwork=TRUE, plotOnlyConnectedNodesNetwork=TRUE, returniGraphs=FALSE, plotType="pdf", labelSize=0.4, fileName="leukemiasNetwork")


