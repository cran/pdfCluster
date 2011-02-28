\name{summary.kepdf-class}
\Rdversion{0.1-4}
\docType{class}
\alias{summary.kepdf-class}
\alias{show,summary.kepdf-method}

\title{Class "summary.kepdf"}
\description{
This class encapsulates results of the application of function {\code{\link{summary.kepdf}}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summary.kepdf", ...)} 
or by a call to the function {\code{\link{summary.kepdf}}}.
}
\section{Slots}{
  \describe{
    \item{\code{obj.class}:}{Object of class \code{"character"} indicating the class of the object whose summary is being done.}
    \item{\code{mode}:}{Object of class \code{"numeric"} indicating the highest density data point. }
    \item{\code{props}:}{Object of class \code{"numeric"} indicating the fraction of the data included within each density level. }
    \item{\code{indices}:}{Object of class \code{"list"} indicating the row or index position of the \code{props} percent 
	top density data points }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.kepdf")}:
	
S4 method for showing objects of \code{\link{summary.kepdf-class}}}}
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
\code{\link{kepdf}}, \code{\link{kepdf-class}}, \code{\link{summary.kepdf}}, \code{\link{show-methods}}.}

\examples{
showClass("summary.kepdf")

data(wine)
x <- wine[,c(2,3,8)]
pdf <- kepdf(x)
summary(pdf)
}

\keyword{classes}
\keyword{nonparametric}
\keyword{smooth}
