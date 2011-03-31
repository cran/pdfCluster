\name{kepdf}
\alias{kepdf}
\title{
Kernel estimate of the probability density function %%  ~~function to do ... ~~
}
\description{
Estimates density of 1-d, 2-d, 3-d, and higher dimensional data by kernel methods. A normal kernel function is used.
}
\usage{
kepdf(x, h = h.norm(x), eval.points = x, h.weights = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector, or a matrix of data whose density should be estimated.  
}
  \item{h}{
a vector of length set to the same dimension as \code{x}, defining the smoothing parameters. It represents
the vector of standard deviations of the kernel function. If this parameter is 
omitted, an asymptotically normal optimal smoothing parameter is used. 
}
  \item{eval.points}{
a vector, or a matrix of data points at which the density estimate should be evaluated.  
}
  \item{h.weights}{
a vector of integers representing the frequencies of individual observations.  
}
}
\details{
The current version of \code{\link{pdfCluster-package}} allows for selecting a diagonal smoothing matrix only when density estimation is performed on multidimensional data.
}
\value{
An S4 object of \code{\link{kepdf-class}}, whose slots are: the data points \code{x}, the values of the density estimate at the evaluation points, the smoothing parameters and
the smoothing parameters' weights. 
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for Data Analysis: Kernel Approach with S-Plus Illustrations}. Oxford University Press, Oxford.

Wand M. P. and Jones, M. C. (1995). \emph{Kernel Smoothing}. Monographs on Statistics and Applied Probability. Chapman & Hall.
}

\seealso{
\code{\link{h.norm}}, \code{\link{kepdf-class}}, \code{\link{plot-methods}}, \code{\link{show-methods}}, 
\code{\link{summary.kepdf-class}}.
}
\examples{
## A 1-dimensional example
data(wine)
x <- wine[,3] 
pdf <- kepdf(x, eval.points=seq(0,7,by=.1))
plot(pdf, n.grid= 100, main="wine data")

## A 2-dimensional example
x <- wine[,c(2,8)] 
pdf <- kepdf(x)
plot(pdf, main="wine data", props=c(5,50,90), ylim=c(0,4))
plot(pdf, main="wine data", method="perspective", phi=30, theta=60)

### A 3-dimensional example
x <- wine[,c(2,3,8)] 
pdf <- kepdf(x)
plot(pdf, main="wine data", props=c(10,50,70), gap=0.2)
plot(pdf, main="wine data", method="perspective", gap=0.2, phi=30, theta=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{nonparametric}