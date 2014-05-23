\name{geneSymbols}
\docType{data}
\alias{geneSymbols}
\title{Gene symbols associated to human Ensemble IDs.}
\description{
Gene symbols to use as gene labels in the package examples.\cr

Source: simplified version of genes.human.annotation from GATExplorer (\url{http://bioinfow.dep.usal.es/xgate/mapping/mapping.php?content=annotationfiles}).
}
\usage{
data(geneSymbols)
}
\format{
Named character vector containing the gene symbol as content, and the associated Ensemble ID as name.
}
\examples{
data(geneSymbols)
head(geneSymbols)
}
