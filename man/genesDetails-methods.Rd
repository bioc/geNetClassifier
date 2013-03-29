\name{genesDetails-methods}
\docType{methods}
\alias{genesDetails}
\alias{genesDetails-methods}
\alias{genesDetails,GenesRanking-method}
\alias{genesDetails.GenesRanking}
\title{Details of the genes in the network.}
\description{
Information of the genes in the ranking (table format).
}
% \usage{
% \method{genesDetails}{GenesRanking}(object, nGenes=NULL, numDecimals=4, classes=NULL, genes=NULL)
% }
\section{Methods}{
\describe{
\code{genesDetails(object, nGenes=NULL, numDecimals=4, classes=NULL, genes=NULL)
}
}}
\arguments{
\item{object}{a GenesRanking}
\item{nGenes}{integer. Number of genes to show per class}
\item{numDecimals}{integer. Number of decimals to show in the numeric values}
\item{classes}{character. Classes of the genes to show}
\item{genes}{character. Genes to show}
 }
\value{
A list containing a dataframe with the details of the genes of each class. For each gene, the following information is provided:
 \item{ranking}{
Ranking of the gene.
}
 \item{gERankMean}{
Mean rank the gene obtained in the cross-validation loops. Only available if geNetClassifier() was called with option \code{estimateGError=TRUE} (False by default).
}
 \item{class}{
Class the gene was chosen for (the class the gene differentiates from the other classes).
}
 \item{postProb}{
Posterior probability which the gene was assigned by the expectation-maximization algorithm (emfit). Tied values are ranked based on the higher  absolute value of exprsMeanDiff. Values are rounded. Several genes may look tied at posterior probability '1' but may actually be i.e. 0.999998 and 0.999997.
}
 \item{exprsMeanDiff}{
Difference betwen the mean expression of the gene within its class and its mean expression in the other classes.
}
 \item{exprsUpDw}{
Gene repressed (DOWN) or over-expressed(UP) for the current class (compared to the other classes).
}
 \item{discriminantPower}{
 Measure calculated based on the coordinates of the support vectors. Represents the weight that the classifier gives to each gene to separate the classes. 
}
 \item{discrPwClass}{
Class for which the Discriminant Power was calculated for. 
}
 \item{isRedundant}{
Does the gene have a high correlation or mutual information with other genes in the list? The threshold to consider a gene redundant can be set through the arguments (by default: correlationsThreshold=0.8 and interactionsThreshold=0.5).
}
}

\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's class (\code{\link[geNetClassifier:GenesRanking-class]{GenesRanking}}) help page.
}
\examples{
data(leukemiasClassifier) # Sample geNetClassifier() return
options(width=200) # Optional, use in case the table rows are wrapped

genesDetails(leukemiasClassifier@classificationGenes)$CML
genesDetails(leukemiasClassifier@genesRanking, nGenes=5, numDecimals=2, 
classes="AML")
genesDetails(leukemiasClassifier@genesRanking, genes=c("ENSG00000096006", 
"ENSG00000168081","ENSG00000105699"))$CLL
}
\keyword{methods}



