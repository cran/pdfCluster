\name{summary.pdfSilhouette}
\alias{summary.pdfSilhouette}

\title{
Summary of objects of class pdfSilhouette
}
\description{
Provides summurizing information about \link{pdfSilhouette-class} objects .}
\usage{
summary.pdfSilhouette(object, ...)
}

\arguments{
  \item{object}{
A \code{\link{pdfSilhouette-class}} object.}

  \item{\dots}{
Not used; for compatibility with S3 methods.
}}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \link{summary.pdfSilhouette-class} object, providing the class of the \code{object} argument and a list whose length
coincides with the number of clusters and whose elements are, for each cluster, the vector of the dbs indexes. 
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
\code{\link{dbs}}, \code{\link{pdfSilhouette-class}}, \code{\link{summary.pdfSilhouette-class}}.
}
\examples{
data(wine)

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
 
dsil <- dbs(cl)
summary(dsil)
}

\keyword{cluster}
\keyword{smooth}

