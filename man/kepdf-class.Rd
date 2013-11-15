\name{kepdf-class}
\Rdversion{1.1}
\docType{class}
\alias{kepdf-class}
\alias{show,kepdf-method}
\alias{summary.kepdf}

\title{Class "kepdf"}
\description{
This class encapsulates results of the application of function {\code{\link{kepdf}}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("kepdf", ...)} or as a result of a call to {\code{\link{kepdf}}}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"call"}, corresponding to the matched call.}
    \item{\code{x}:}{Object of class \code{"matrix"} representing the data points used to estimate the probability density function.}
    \item{\code{eval.points}:}{Object of class \code{"matrix"} representing the data points at which the density is evaluated.}
      \item{estimate }{The values of the density estimate at the evaluation points.}     
	\item{\code{kernel}: }{Object of class \code{"character"} giving the selected kernel.}
    \item{\code{bwtype}: }{Object of class \code{"character"} giving the selected type of estimator.} 
    \item{\code{par}: }{Object of class \code{"list"} providing the parameters used to estimate the density. Its elements are \code{h}, \code{hx},
	and possibly \code{alpha}.
	}
  }
See \code{\link{kepdf}} for further details.
}
\section{Methods}{
\describe{
    \item{plot}{\code{signature(x = "kepdf", y = "ANY")} 
	
	Plots objects of {\code{\link{kepdf-class}}}. 
	{\code{\link{plot-methods}}} are available for density estimates of:
	\itemize{
		\item {one-dimensional data;}
		\item {two-dimensional data: contour, image or perspective plots are available;}
		\item {multi-dimensional data: matrix of plots of all the pairs of two-dimensional marginal kernel density estimates.}
		}
		See {\code{\link{plot,kepdf-method}}} for further details.	}
    \item{show}{\code{signature(object = "kepdf")} 
	
	Prints the following elements:
		\itemize{
		\item {the class of the object;}
		\item {the selected kernel;}
		\item {the selected type of estimator;}
		\item {either the fixed smoothing parameters or the smoothing parameters of each observation;}
		\item {the density estimates at the evaluation points.}
		}	}
	\item{summary}{\code{signature(object = "kepdf")} 
	
    Provides a summary of \code{\link{kepdf-class}} object by printing
	the highest density data point and the row or index position of a 
	percentage top density data points, possibly given as optional argument
	\code{prop}.     
		}
	
	 }
}


\seealso{
\code{\link{h.norm}}, \code{\link{kepdf}}, {\code{\link{plot,kepdf-method}}},
\code{\link{plot-methods}}, \code{\link{show-methods}}, \code{\link{summary-methods}}.
}
\examples{
#
showClass("kepdf")

#
data(wine)
#select only "Barolo"-type wines
x <- wine[1:59,3] 
pdf <- kepdf(x)
pdf
summary(pdf)
summary(pdf, props = 10*seq(1, 9, by = 1))

}

\keyword{classes}
\keyword{nonparametric}
\keyword{smooth}
