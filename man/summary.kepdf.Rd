\name{summary.kepdf}
\alias{summary.kepdf }
\title{
Summary of objects of class kepdf.
}
\description{
Provides summurizing information about \link{kepdf-class} objects .}
\usage{
summary.kepdf(object, ..., props = c(75, 50, 25))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{\link{kepdf-class}} object.}

  \item{\dots}{
not used; for compatibility with S3 methods;
}

  \item{props}{
A vector defining the fraction of the data to be included within each density level. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \link{summary.kepdf-class} object, providing the class of the \code{object} argument, 
the highest density data point and the row or index position of the \code{props} percent top density data points. 
 }

%% \references{
%% ~put references to the literature site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{kepdf}}, \code{\link{kepdf-class}}, \code{\link{summary.kepdf-class}}.
}
\examples{
data(wine)
x <- wine[,c(2,3,8)]
pdf <- kepdf(x)
summary(pdf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{smooth}

