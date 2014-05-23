\name{network2txt}
\docType{methods}
\alias{network2txt}
\alias{network2txt-methods}
\alias{network2txt,GenesNetwork-method}
\alias{network2txt,list-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
network2txt
}
\description{
Saves the GenesNetwork as text file.
}
\usage{
network2txt(network, filePrefix = NULL, nwClass = NULL)
}
\arguments{
 \item{network}{
GenesNetwork or list of GenesNetworks.
}
\item{filePrefix}{
Character. File name prefix.
}
\item{nwClass}{
Character. Network class.
}
}
\value{
Saves the networks as text (.txt) files.
The files will be saved in the current working directory as \code{filePrefix}_\emph{className}.txt.
}

\seealso{
Main package function and classifier training:
\code{\link[geNetClassifier:geNetClassifier]{geNetClassifier}} \cr
This method's class (\code{\link[geNetClassifier:GenesNetwork-class]{GenesNetwork}}) help page.
}

\examples{
## Load or calculate a network:

data(leukemiasClassifier)

## Export as text:
network2txt(leukemiasClassifier@genesNetwork, filePrefix="leukemiasNetwork")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
