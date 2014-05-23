##
## plotExpressionProfiles()
## function(eset, genes=NULL, fileName=NULL, geneLabels=NULL, sampleLabels=NULL, labelsOrder=NULL, showSampleNames=FALSE, showMean= FALSE, sameScale=TRUE, verbose=TRUE)
##
## Plots the expression of the given genes in all the given samples
# Input:  eSet, gene list...
# Output: Plot

test_plotExpressionProfiles <- function()
{
    myEset <- rbind(1:4, c(rep(1,2), rep(2,2)))
    rownames(myEset) <- c("gene1","gene2")

    # Nothing should be returned (NULL)
    checkEquals(plotExpressionProfiles(eset=myEset, genes=rownames(myEset), sampleLabels=c(rep("one", 2), rep("two", 2))), NULL)
}
