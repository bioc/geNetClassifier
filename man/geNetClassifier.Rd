\name{geNetClassifier}
\alias{geNetClassifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function of the geNetClassifier package.\cr
Trains the multi-class SVM classifier based on the given gene expression data through transparent detection of gene markers and their associated networks.
}
\description{
Allows to train the classifier, calculate the genes network...
}
\usage{
geNetClassifier(eset, sampleLabels, plotsName = NULL,
buildClassifier = TRUE, estimateGError = FALSE,
calculateNetwork = TRUE, labelsOrder = NULL, geneLabels = NULL,
numGenesNetworkPlot = 100,
minGenesTrain = 1, maxGenesTrain = 100, continueZeroError = FALSE,
numIters = 6, lpThreshold = 0.95, numDecimals = 3,
removeCorrelations = FALSE, correlationsThreshold = 0.8,
correlationMethod = "pearson",
removeInteractions = FALSE, interactionsThreshold = 0.5,
minProbAssignCoeff = 1, minDiffAssignCoeff = 0.8,
IQRfilterPercentage = 0, skipInteractions = TRUE,
precalcGenesNetwork = NULL, precalcGenesRanking = NULL,
returnAllGenesRanking = TRUE, kernel="linear", verbose=TRUE, ...)
}
\arguments{
 \item{eset}{
ExpressionSet or matrix. Gene expression of the train samples (positive & non-logaritmic normalized values).
}
  \item{sampleLabels}{
Character. PhenoData variable (column name) containing the train samples class labels. \cr
Matrix or Factor. Class labels of the train samples.
}
  \item{labelsOrder}{
Vector or Factor. Order in which the labels should be shown in the returned results and plots.
}
  \item{plotsName}{
Character. File name with which the plots should be saved. If not provided, no plots will be drawn.
}
  \item{buildClassifier}{
Logical. If TRUE trains a classifier with the given samples.
}
  \item{estimateGError}{
Logical. If TRUE uses cross-validation to estimate the Generalization Error of a classiffier trained with the given samples.
}
  \item{calculateNetwork}{
Logical. If TRUE calculates the coexpression network between the best genes.
}
  \item{geneLabels}{
Vector or Matrix. Gene name, ID or label which should be shown in the returned results and plots.
}
  \item{numGenesNetworkPlot}{
Integer. Number of genes to show in the coexpression network for each class.
}
  \item{minGenesTrain}{
Integer. Maximum number of genes per class to train the classifier with.
}
  \item{maxGenesTrain}{
Integer. Maximum number of genes per class to train the classifier with.
}
  \item{continueZeroError}{
Logical. If TRUE, the program will continue testing combinations with more genes even if error 0 has been reached.
}
  \item{numIters}{
Integer. Number of iterations to determine the optimum number of genes (between \code{minGenesTrain} and \code{maxGenesTrain}).
}
  \item{lpThreshold}{
Numeric between 0 and 1. Required posterior probability value to consider a gene 'significant'.
}
  \item{removeCorrelations}{
Logical. If TRUE, no correlated genes will be chosen to train the classifier.
}
  \item{correlationsThreshold}{
Numeric between 0 and 1. Minimum Pearson's correlation coefficient to consider genes correlated.
}
  \item{correlationMethod}{
"pearson", "kendall" or "spearman". Type of correlation to calculate between genes.
}
  \item{removeInteractions}{
Logical. If TRUE, genes with Mutual Information coefficient over the threshold will not be chosen to train the classifier.
}
  \item{interactionsThreshold}{
Numeric between 0 and 1. Minimum Mutual Information coefficient to consider two genes equivalent.
}
  \item{numDecimals}{
Integer. Number of decimals to show in the statistics.
}
  \item{minProbAssignCoeff}{
Numeric. Allows modifying the required probability to assign a sample to a class in the internal crossvalidation. For details see: \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}}
}
  \item{minDiffAssignCoeff}{
Numeric. Allows modifying the difference of probabilities required between the most likely class and second most likely class to assign a sample. For details see: \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}}
}
  \item{IQRfilterPercentage}{
Integer. InterQuartile Range (IQR) filter applied to the initial data. Not recommended for more than two classes.
}
  \item{skipInteractions}{
Logical. If TRUE, the interactions between genes are not calculated (they will not appear on the genes network). Saves some execution time. Only available if \code{removeInteractions=FALSE}.
}
  \item{precalcGenesNetwork}{
\code{GenesNetwork} from a previous execution with the same expression data and parameters.
}
  \item{precalcGenesRanking}{
\code{GenesRanking} from a previous execution with the same expression data and parameters.
}
  \item{returnAllGenesRanking}{
Logical. If TRUE, returns the whole genes ranking. If FALSE the returned ranking contains only the significant genes (genes over lpThreshold).
}
  \item{verbose}{
Logical. If TRUE, messages indicating the execution progress will be shown.
}
  \item{kernel}{
Character. Type of SVM kernel. Default: "linear",
}
  \item{...}{
Other arguments to pass to the \code{\link{svm}} function.
}

}

\value{
A \code{\link[geNetClassifier:GeNetClassifierReturn-class]{GeNetClassifierReturn}} object containing the classifier and the genes chosen to train it (\code{classificationGenes}), Cross-Validation statistics, the whole \code{GenesRanking} and each class' \code{GenesNetwork} (if requested).
Several plots saved as '\code{plotsName_....pdf}' in the working directory.
}

\author{
Bioinformatics and Functional Genomics Group. Centro de Investigacion del Cancer (CIC-IBMCC, USAL-CSIC). Salamanca. Spain
}
% \note{
% TO DO:
% Any further notes?
% i.e. Explicacion de LP 1/4?
% }

\references{
% TO DO:
% HAy que poner alguna referencia?
% i.e. Paquetes que hemos usado?
Packages used by this function:\cr
\code{\link[EBarrays:emfit]{EBarrays: emfit}} (Implements EM algorithm for gene expression mixture model) and ebPatterns, for calculating the gene ranking.\cr
 Ming Yuan, Michael Newton, Deepayan Sarkar and Christina Kendziorski (2007). EBarrays: Unified Approach for Simultaneous Gene Clustering and Differential Expression Identification. R package.\cr
 \cr
\code{\link[e1071:svm]{e1071: svm}}.\cr
Evgenia Dimitriadou, Kurt Hornik, Friedrich Leisch, David Meyer and Andreas Weingessel (2011). e1071: Misc Functions of the Department of Statistics (e1071), TU Wien. R package.\cr
http://CRAN.R-project.org/package=e1071\cr
\cr
\code{\link[ipred:kfoldcv]{ipred: kfoldcv}} (computes feasible sample sizes for the k groups in k-fold cv) for the cross-validations.\cr
  Andrea Peters and Torsten Hothorn (2012). ipred: Improved Predictors. R package. http://CRAN.R-project.org/package=ipred\cr
\cr
\code{\link[minet:minet]{minet}} for the Mutual Information network.\cr
Patrick E. Meyer, Frederic Lafitte and Gianluca Bontempi (2008). MINET: An open source R/Bioconductor Package for Mutual Information based Network Inference. BMC Bioinformatics.\cr
http://www.biomedcentral.com/1471-2105/9/461\cr
\cr
\code{\link{RColorBrewer}} (\code{brewer.pal}) for palettes in some of the plots.\cr
  Erich Neuwirth (2011). RColorBrewer: ColorBrewer palettes. R package. \cr
  http://CRAN.R-project.org/package=RColorBrewer\cr
\cr
\code{\link{igraph}} for the graphical representation of the networks.\cr
  Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. http://igraph.sf.net\cr
}

\seealso{
To query the classifier: \code{\link[geNetClassifier:queryGeNetClassifier]{queryGeNetClassifier}}\cr
All functions in the package: \code{\link[geNetClassifier:geNetClassifier-package]{geNetClassifier-package}}

}
\examples{
########
# Load libraries and training data
########

# Load an expressionSet:
library(leukemiasEset)
data(leukemiasEset)

# Select the train samples:
trainSamples<- c(1:10, 13:22, 25:34, 37:46, 49:58)
# summary(leukemiasEset$LeukemiaType[trainSamples])


########
# Training
########

# NOTE: Training the classifier takes a while...
# Choose ONE of the followings, or modify to suit your needs:
\dontrun{

# "Basic" execution: All default parameters
leukemiasClassifier <- geNetClassifier(eset=leukemiasEset[,trainSamples],
    sampleLabels="LeukemiaType", plotsName="leukemiasClassifier")

# All default parameters also estimatings the classiffier's Generalization Error:
# ( by default:  buildClassifier=TRUE, calculateNetwork=TRUE)
# Takes longer time than the basic execution
leukemiasClassifier <- geNetClassifier(eset=leukemiasEset[,trainSamples],
    sampleLabels="LeukemiaType", plotsName="leukemiasClassifier",
    estimateGError=TRUE)

# Faster execution (few minutes - depending on the computer):
# By skipping the calculation of the interactions (MI) betwen the genes,
# and reducing the number of genes to explore when training the classifier
# (100 by default), the execution time can be sightly reduced
leukemiasClassifier <- geNetClassifier(eset=leukemiasEset[,trainSamples],
sampleLabels="LeukemiaType", plotsName="leukemiasClassifier",
skipInteractions= TRUE, maxGenesTrain=20)

# To any of these examples, you can add/remove the argument geneLabels,
# in order to show/remove the gene name in the rankings and plots:
# The argument labelsOrder allows showing the classes in a specific order
# i.e.: labelsOrder=c("ALL","CLL","AML",CML","NoL")

save(leukemiasClassifier, file="leukemiasClassifier.RData")  # Save execution result
# For loading the saved object in the future...
# (If it doesn't find it, use getwd() to make sure you are in the right directory)
#load("leukemiasClassifier.RData")


# To avoid having to train a classifier to continue learning to use the package,
# you can load the package's pre-executed example:
data(leukemiasClassifier)
#This example classifier was trained with the following code:
#leukemiasClassifier <- geNetClassifier(leukemiasEset[,trainSamples],
#    "LeukemiaType", plotsName="leukemiasClassifier", buildClassifier=TRUE,
#    estimateGError=TRUE, calculateNetwork=TRUE, geneLabels=geneSymbols)

########
# Explore the returned object:
########
names(leukemiasClassifier)
# More details on the class' help file:
?GeNetClassifierReturn

# Further options:
# The trained classifier can be used to find the class of new samples:
?queryGeNetClassifier

# The default plots can be modified and presonalized to fit the user needs:
?calculateGenesRanking
?plotNetwork
?plotDiscriminantPower
?plotExpressionProfiles
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}
