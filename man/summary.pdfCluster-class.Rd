\name{summary.pdfCluster-class}
\Rdversion{0.1-4}
\docType{class}
\alias{summary.pdfCluster-class}
\alias{show,summary.pdfCluster-method}

\title{Class "summary.pdfCluster"}
\description{
This class encapsulates results of the application of function {\code{\link{summary.pdfCluster}}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summary.pdfCluster", ...)} or as a result of a call to {\code{\link{summary.pdfCluster}}}.
}
\section{Slots}{
   \describe{
    \item{\code{obj.class}:}{Object of class \code{"character"} indicating the class of the object whose summary is being done.}
    \item{\code{cluster.cores}:}{Object of class \code{"numeric"} defining the cluster cores membership of allocated data.}
	\item{\code{clusters}:}{Object of class \code{"ANY"} reporting the final group labels.}
	\item{\code{tree}:}{Object of class \code{"ANY"} representing the hierarchical structure corresponding to the clusters associated to
	different sections of the density estimate. The object has class \code{"dendrogram"} 
	when the tree has more than one leaf.}  
  }
}
\section{Methods}{
\describe{
    \item{show}{\code{signature(object = "summary.pdfCluster")} 
	
S4 method for showing objects of \code{\link{summary.kepdf-class}}, printing the following elements:
		\itemize{
		\item {the class of the object;}
		\item {the frequency table of the cluster cores}
		\item {the frequency table of the final grouping}
		}	}
	 }
}
%%\references{
%%  ~~put references to the literature/web site here~~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfCluster}}, \code{\link{pdfCluster-class}}, \code{\link{show-methods}}, \code{\link{summary.pdfCluster}}.
}
\examples{
#
showClass("pdfCluster")

#
# load data
data(wine)
gr <- wine[,1]

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
summary(cl)
}

\keyword{classes}
\keyword{cluster}
\keyword{multivariate}
\keyword{smooth}
