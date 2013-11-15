\name{pdfCluster-class}
\Rdversion{1.1}
\docType{class}
\alias{pdfCluster-class}
%\alias{dbs,pdfCluster,missing-method}
%\alias{plot,pdfCluster,missing-method}
\alias{show,pdfCluster-method}
\alias{summary.pdfCluster}

\title{Class "pdfCluster"}
\description{
This class pertains to results of the application of function {\code{\link{pdfCluster}}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pdfCluster", ...)} or as a result to a call to \code{\link{pdfCluster}}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"call"} representing the matched call;}
    \item{\code{x}:}{Object of class \code{"matrix"} representing the clustered data points; }
    \item{\code{pdf}:}{Object of class \code{"list"} reporting details about the kernel density estimate at data points \code{x}. }
    \item{\code{nc}:}{Object of class \code{"list"} summarizing the result of the connected components search for different sections of
		the estimated density.}
	 \item{\code{graph}:}{An object of class \code{"list"} defining details about the graph built to find the connected sets 
	 of high density regions.} 
    \item{\code{cluster.cores}:}{Object of class \code{"ANY"} reporting the group labels of the data allocated to the cluster cores. }
    \item{\code{tree}:}{Object of class \code{"ANY"}, namely class \code{dendrogram} if the procedure detects more than one group, \code{list} otherwise.
		It reports the cluster tree structure associated to the different connected components for different density levels.}
    \item{\code{noc}:}{Object of class \code{"numeric"} giving the number of clusters. }
	\item{\code{stages}:}{Object of class \code{"ANY"}, being \code{NULL} if the cluster cores only are detected, \code{"list"} when also the lower density data are allocated. 
							The elements of the list correspond to the group labels at the different stages of the classification procedure.
							\code{NA} values correspond to unlabeled data.  }
    \item{\code{clusters}:}{Object of class \code{"ANY"} being \code{NULL} if the cluster cores only are detected, \code{"numeric"} when all the data are clustered. This slot is obsolete. Groups can be extracted by a call to function \code{\link{groups}}.}
  }
See \code{\link{pdfCluster}} for further details.  
}
\section{Methods}{
  \describe{
    \item{dbs}{\code{signature(x = "pdfCluster", clusters = "missing")}
	
	Computes the density based Silhouette diagnostics of 
	clustered data. See \code{\link{dbs}} for further details.	}
    
	\item{pdfCluster}{\code{signature(x="pdfCluster")} 
    
	Speeds up time for re-running the \code{\link{pdfCluster}} procedure with different values of 
	\code{tau} when \code{graphtype = "pairs"}}
	
	\item{plot}{\code{signature(x = "pdfCluster", y = "missing")} 
	
	Plots objects of {\code{\link{pdfCluster-class}}}. 
	{\code{\link{plot-methods}}} are available for:
	\itemize{
		\item {the mode function: gives the number of connected components when the proportion of data points with density above a threshold varies. 
				Set argument \code{which} to 1 to display this plot.}
		\item {the cluster tree: plot the hierarchical structure associated to the clusters detected by different sections of
			the density estimate. Set argument \code{which} to 2 to display this plot.  }
		\item {the data points: scatterplot of data or of all the possible couples of coordinates reporting the label group. Set 
			argument \code{which} to 3 to display this plot. }
		\item {the density-based Silhouette information: graphical diagnostics
		of the clustering. See \code{\link{plot,dbs-method}}. 
		   Set argument \code{which} to 4 to display this plot. Not available when \code{noc}=1.}
		} See \code{\link{plot,pdfCluster-method}} for further details.}
    \item{show}{\code{signature(object = "pdfCluster")}.

	Prints the following elements:
		\itemize{
		\item{the matched Call;}
		\item{the type of kernel estimator;} 
		\item{the type of graph built;}  
		\item{the groups tree (if available);}
		\item{the cluster cores;} 
		\item{the cluster labels at the different stages of the classification procedure;}
		\item{the final clustering.}
		} }
	
    \item{summary}{\code{signature(object = "pdfCluster")}.

	Provides a summary of \code{\link{pdfCluster-class}} objects by printing
    the following elements: 	
		\itemize{
		\item {the matched call to pdfCluster function}
		\item {the frequency table of the cluster cores;}
		\item {the frequency table of the final grouping;}
		\item {the tree of clusters.}
		} }	
	 }
}
%\references{
%Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation.
% \emph{Statistics and Computing}. vol. 17, pp. 71-80.

%Menardi G., Azzalini, A. (2013). An advancement in clustering via nonparametric
% density estimation. DOI: 10.1007/s11222-013-9400-x.
%}

\seealso{
\code{\link{pdfCluster}}, \code{\link{plot,pdfCluster-method}}, 
\code{\link{show-methods}}, \code{\link{summary-methods}}   
}
\examples{
showClass("pdfCluster")

data(wine)
x <-wine[ ,-1]
gr <- wine[ ,1]

# clustering
cl <- pdfCluster(x, graphtype="pairs", bwtype="adaptive")
summary(cl)
cl
plot(cl)
}
\keyword{classes}
\keyword{cluster}
\keyword{multivariate}
