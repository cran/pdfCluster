\name{groups}
\alias{groups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extracts groups}
\description{
Extracts the detected groups from objects of \code{\link{pdfCluster-class}}. 
}
\usage{
groups(obj, stage = length(obj@stages))
}

\arguments{
  \item{obj}{
An object of \code{\link{pdfCluster-class}}
}
	\item{stage}{
The stage of classification at which the clusters have to be extracted. Set this value to 0
to extract the cluster cores. Default value is the total number of classification stages, that is, 
the final partition is given. When \code{obj} contains the clusters cores only, these are given by default.
}
}
\details{
This function is an user-friendly version of command \code{obj@clusters}, now obsolete, to ease
extraction of groups from objects of \code{\link{pdfCluster-class}}.
}
\value{
A numeric vector containing the group labels. NA values are associated to points not classified at the
selected stage of the classification procedure.   
}
%\references{
%L. Hubert and P. Arabie (1985) Comparing Partitions, \emph{Journal of the Classification}, 2, 193-218. } 

\seealso{
\code{\link{pdfCluster}}}

\examples{
# load data
data(wine)

# select a subset of variables
x <- wine[, c(2, 5, 8)]

#clustering
cl <- pdfCluster(x)

groups(cl)

#equivalent to:
cl@clusters

#to extract the cluster cores
groups(cl, stage = 0)
}
\keyword{cluster}
