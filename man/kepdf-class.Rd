\name{kepdf-class}
\Rdversion{1.1}
\docType{class}
\alias{kepdf-class}
\alias{show,kepdf-method}

\title{Class "kepdf"}
\description{
This class encapsulates results of the application of function {\code{\link{kepdf}}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("kepdf", ...)} or as a result of a call to {\code{\link{kepdf}}}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"matrix"} representing the data points used to estimate the probability density function.}
    \item{\code{eval.points}:}{Object of class \code{"matrix"} representing the data points at which the density is evaluated.}
    \item{\code{estimate}:}{Object of class \code{"numeric"} corresponding to the kernel density estimate at the \code{eval.points}.}
    \item{\code{h}:}{Object of class \code{"numeric"} with the same length as the number of columns of \code{x} and corresponding to 
	the diagonal elements of the smoothing matrix used to estimate the density.}
    \item{\code{h.weights}:}{Object of class \code{"numeric"} representing the frequencies of individual observations to be used to estimate the density. }
    \item{\code{call}:}{Object of class \code{"call"}, corresponding to the matched call.}
  }
See \code{\link{kepdf}} for further details.
}
\section{Methods}{
\describe{
    %%\item{lines}{\code{signature(x = "kepdf")} add graphical display of a {\code{\link{kepdf-class}}} 
	%%object to an existing plot. To be used for displaying kernel density estimate of one-dimensional data only.}
    \item{plot}{\code{signature(x = "kepdf", y = "ANY")} 
	
	Plots objects of {\code{\link{kepdf-class}}}. 
	{\code{\link{plot-methods}}} are available for density estimates of:
	\itemize{
		\item {one-dimensional data;}
		\item {two-dimensional data: contour, image or perspective plots are available;}
		\item {multi-dimensional data: matrix of plots of all the couples of 
		two-dimensional marginal kernel density estimates.}
		}. See {\code{\link{plot.kepdf}}} for further details.	}
    \item{show}{\code{signature(object = "kepdf")} 
	
	Prints the following elements:
		\itemize{
		\item {the class of the object;}
		\item {The diagonal elements of the smothing matrix;}
		\item {The density estimates at the evaluation points.}
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
\code{\link{h.norm}}, \code{\link{kepdf}}, \code{\link{plot-methods}}, \code{\link{show-methods}}, \code{\link{summary.kepdf-class}}.
}
\examples{
#
showClass("kepdf")

#
data(wine)
#select only "Barolo"-type wines
x <- wine[1:59,3] 
pdf <- kepdf(x)
summary(pdf)
pdf
}

\keyword{classes}
\keyword{nonparametric}
\keyword{smooth}
