\name{pdfSilhouette-class}
\Rdversion{1.1}
\docType{class}
\alias{pdfSilhouette-class}
%\alias{plot,pdfSilhouette,missing-method}
\alias{show,pdfSilhouette-method}

\title{Class "pdfSilhouette"}
\description{
This class encapsulates results of the application of function {\code{\link{dbs}}}.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pdfSilhouette", ...)} or as a result from calling function \code{\link{dbs}}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"matrix"} representing the clustered data points. }
    \item{\code{prior}:}{Object of class \code{"numeric"} being the prior probabilities of belonging to the groups. }
    \item{\code{dbs}:}{Object of class \code{"numeric"} reporting the density-based Silhouette information of the clustered data. }
    \item{\code{clusters}:}{Object of class \code{"numeric"} reporting the group labels of grouped data. }
    \item{\code{nc}:}{Object of class \code{"numeric"} indicating the number of clusters.}
    \item{\code{stage}:}{Object of class \code{"ANY"} corresponding to the stage of the classification at which the density-based Silhouette information
	is computed when \code{\link{dbs}} is applied to an object of \code{\link{pdfCluster-class}}. }
    \item{\code{call}:}{Object of class \code{"call"} reporting the matched call. }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "pdfSilhouette", y = "missing")}: 

S4 method for plotting objects of \code{\link{pdfSilhouette-class}}. Data are partitioned 
into the clusters, sorted in a decreasing order with respect to their dbs value and displayed
on a bar graph. See \code{\link{plot.pdfSilhouette}} for further details.}
    \item{show}{\code{signature(object = "pdfSilhouette")}: 
	
S4 method for showing objects of \code{\link{pdfSilhouette-class}}. The following elements are shown:
		\itemize{
		\item {the dbs index computed at the observed data;}
		\item {The cluster membership of each data point;}
		}	}
	 }
}
\references{
Menardi, G. (2010) Density-based Silhouette diagnostics for clustering methods. \emph{Statistics and Computing}, DOI: 10.1007/s11222-010-9169-0.}

%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}


\seealso{
\code{\link{dbs}}, \code{\link[cluster:silhouette]{silhouette}},
\code{\link{summary.pdfSilhouette}}, \code{\link{summary.pdfSilhouette-class}}. 
}
\examples{
showClass("pdfSilhouette")

#wine example
#data loading
data(wine)

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
 
dsil <- dbs(cl)
dsil
}
\keyword{classes}
\keyword{cluster}
\keyword{smooth}

