\name{summary.pdfSilhouette-class}
\Rdversion{1.1}
\docType{class}
\alias{summary.pdfSilhouette-class}
\alias{show,summary.pdfSilhouette-method}

\title{Class "summary.pdfSilhouette"}
\description{
This class encapsulates results of the application of function {\code{\link{summary.pdfSilhouette}}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summary.pdfSilhouette", ...)}
or as a result of a call to {\code{\link{summary.pdfSilhouette}}}.}

\section{Slots}{
  \describe{
    \item{\code{obj.class}:}{Object of class \code{"character"} indicating the class of the object whose summary is being done.}
    \item{\code{dbs.groups}:}{Object of class \code{"list"}. Each element of the list corresponds to a group and reports the dbs 
	indexes for each element of the group. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.pdfSilhouette")}: 

S4 method for showing objects of \code{\link{summary.pdfSilhouette-class}}}
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
\code{\link{dbs}}, \code{\link{summary.pdfSilhouette}}, \code{\link{pdfSilhouette-class}}.
}
\examples{
showClass("summary.pdfSilhouette")

data(wine)

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
 
dsil <- dbs(cl)
summary(dsil)

}
\keyword{classes}
\keyword{cluster}
\keyword{smooth}
