\name{pdfCluster-class}
\Rdversion{1.1}
\docType{class}
\alias{pdfCluster-class}
%\alias{dbs,pdfCluster,missing-method}
%\alias{plot,pdfCluster,missing-method}
\alias{show,pdfCluster-method}

\title{Class "pdfCluster"}
\description{
This class encapsulates results of the application of function {\code{\link{pdfCluster}}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pdfCluster", ...)} or as a result to a call to \code{\link{pdfCluster}}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"call"} representing the matched call;}
    \item{\code{x}:}{Object of class \code{"matrix"} representing the clustered data points; }
    \item{\code{estimate}:}{Object of class \code{"numeric"} reporting the kernel density estimates at data points \code{x}. }
    \item{\code{h}:}{Object of class \code{"numeric"} reporting the smoothing parameters used to estimate the density of data. }
    \item{\code{hmult}:}{Object of class \code{"numeric"} representing the multiplying factor of the smoothing parameters. }
    \item{\code{nc}:}{Object of class \code{"list"} summarizing the result of the connected components search for different sections of
		the estimated density.}
    \item{\code{cluster.cores}:}{Object of class \code{"ANY"} reporting the group labels of the data allocated to the cluster cores. }
    \item{\code{tree}:}{Object of class \code{"ANY"}, namely class \code{dendrogram} if the procedure detects more than one group, \code{list} otherwise.
		It reports the cluster tree structure associated to the different connected components for different density levels.}
    \item{\code{stages}:}{Object of class \code{"ANY"}, being NULL if the cluster cores only are detected, \code{"list"} when also the lower density data are allocated. 
							The elements of the list correspond to the group labels at the different stages of the classification procedure.
							NA values correspond to unlabeled data.  }
    \item{\code{ctrl}:}{Object of class \code{"logical"} being TRUE if the procedure detects one group only. }
	\item{\code{clusters}:}{Object of class \code{"ANY"} being NULL if the cluster cores only are detected, \code{"numeric"} when
	all the data are clustered.}
  }
See \code{\link{pdfCluster}} for further details.  
}
\section{Methods}{
  \describe{
    \item{dbs}{\code{signature(x = "pdfCluster", clusters = "missing")}
	
	Computes the density based Silhouette diagnostics of 
	clustered data. See \code{\link{dbs}} for further details.	}
    \item{plot}{\code{signature(x = "pdfCluster", y = "missing")} 
	
	Plots objects of {\code{\link{pdfCluster-class}}}. 
	{\code{\link{plot-methods}}} are available for:
	\itemize{
		\item {the mode function, gives the number of connected components when the proportion of data points with density above a threshold varies. 
				Set argument \code{which} to 1 to display this plot.}
		\item {the cluster tree: plot the hierarchical structure associated to the clusters detected by different sections of
			the density estimate. Set argument \code{which} to 2 to display this plot. Not available when \code{ctrl}=TRUE. }
		\item {the data points: scatterplot of data or of all the possible couples of coordinates reporting the label group. Set 
			argument \code{which} to 3 to display this plot. }
		}. See \code{\link{plot.pdfCluster}} for further details.}
    \item{show}{\code{signature(object = "pdfCluster")}.

	Prints the following elements:
		\itemize{
		\item {the smoothing parameters used to estimate the probability density function;} 
		\item {the density estimates at data points;}  
		\item{the groups tree (if available);}
		\item{the cluster cores;} 
		\item{the cluster labels at the different stages of the classification procedure;}
		\item{the final clustering;}
		} }
	 }
}
\references{
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation. \emph{Statistics and Computing}. vol. 17, pp. 71-80.
}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dbs}}, \code{\link{pdfCluster}}, \code{\link{plot-methods}}, \code{\link{show-methods}}, \code{\link{summary.pdfCluster-class}}  
}
\examples{
showClass("pdfCluster")
}
\keyword{classes}
\keyword{cluster}
\keyword{multivariate}
\keyword{smooth}
