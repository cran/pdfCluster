\name{dbs-class}
\Rdversion{1.1}
\docType{class}
\alias{dbs-class}
\alias{summary.dbs}
%\alias{plot,dbs,missing-method}
\alias{show,dbs-method}

\title{Class "dbs"}
\description{
This class pertains to results of the application of function 
{\code{\link{dbs}}}.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("dbs", ...)} or as a 
result from calling function \code{\link{dbs}}.
}
\section{Slots}{
  \describe{
 \item{\code{call}:}{Object of class \code{"call"} reporting the matched 
	call. }   
 \item{\code{x}:}{Object of class \code{"matrix"} representing the clustered
	data points. }
    \item{\code{prior}:}{Object of class \code{"numeric"} being the prior 
	probabilities of belonging to the groups. }
    \item{\code{dbs}:}{Object of class \code{"numeric"} reporting the 
	density-based silhouette information of the clustered data. }
    \item{\code{clusters}:}{Object of class \code{"numeric"} reporting the 
	group labels of grouped data. }
    \item{\code{noc}:}{Object of class \code{"numeric"} indicating the number 
	of clusters.}
    \item{\code{stage}:}{Object of class \code{"ANY"} corresponding to the 
	stage of the classification at which the density-based silhouette information
	is computed when \code{\link{dbs}} is applied to an object of 
	\code{\link{pdfCluster-class}}. }
     }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "dbs", y = "missing")}: 

S4 method for plotting objects of \code{\link{dbs-class}}. Data are partitioned
into the clusters, sorted in a decreasing order with respect to their dbs value
and displayed on a bar graph. See \code{\link{plot,dbs-method}} for further details.}
    \item{show}{\code{signature(object = "dbs")}: 
	
S4 method for showing objects of \code{\link{dbs-class}}. The following 
elements are shown:
		\itemize{
		\item {the dbs index computed at the observed data;}
		\item {The cluster membership of each data point;}
		}	}
	 \item{summary}{\code{signature(object = "dbs")}: 
	
S4 method for summarizing objects of \code{\link{dbs-class}}. The following
elements are shown:
		\itemize{
		\item {a summary (minimum, 1st quartile, median, mean, 3rd quartile, 
		maximum) of the dbs values for each cluster;}
		\item {a summary (minimum, 1st quartile, median, mean, 3rd quartile, 
		maximum) of the dbs values for all the observations.}
		}	}
	 }	 
}
%\references{
%Menardi, G. (2011) Density-based Silhouette diagnostics for clustering methods.
% \emph{Statistics and Computing}, 21:295-308.}


\seealso{
\code{\link{dbs}}, \code{\link[cluster:silhouette]{silhouette}},
\code{\link{plot,dbs-method}}, \code{\link{plot-methods}}, 
\code{\link{show-methods}}, \code{\link{summary-methods}}. 
}

\examples{
showClass("dbs")

#wine example
#data loading
data(wine)

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
 
dsil <- dbs(cl)
dsil
summary(dsil)
}

\keyword{classes}
\keyword{cluster}

