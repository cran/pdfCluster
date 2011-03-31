\name{summary.pdfCluster}
\alias{summary.pdfCluster}


\title{
Summary of objects of class pdfCluster
}
\description{
Provides a summary of objects having \code{\link{pdfCluster-class}}. }

\usage{
summary.pdfCluster(object, ...)
}

\arguments{
  \item{object}{
An object of \code{\link{pdfCluster-class}}
}

  \item{\dots}{
not used; for compatibility with S3 methods;
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \link{summary.pdfCluster-class} object, providing the class of the \code{object} argument, 
the highest density data point and the row or index position of the \code{props} percent top density data points. 
}

%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfCluster}}, \code{\link{pdfCluster-class}}, \code{\link{summary.pdfCluster-class}}.
}
\examples{
# load data
data(wine)
gr <- wine[,1]

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
summary(cl)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{multivariate}
\keyword{smooth}

