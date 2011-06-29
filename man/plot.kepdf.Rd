\name{plot.kepdf}
\alias{plot.kepdf1d}
\alias{plot.kepdf2d}
\alias{plot.kepdfdd}
\alias{plot.kepdf}
\alias{plot,kepdf,missing-method}
\title{
Plot objects of class kepdf
}
\description{
Function and methods for plotting objects of \code{\link{kepdf-class}}
}

\usage{
plot.kepdf1d(x, y, eval.points=NULL, n.grid=50, data=NULL, 
	add=FALSE, main	= NULL, xlab=NULL, ylab=NULL, col=NULL,
	col.data=2, type="l", ...)

plot.kepdf2d(x, y, eval.points=NULL, n.grid=c(30,30), data=NULL,
	add=FALSE, main=NULL, xlab=NULL, ylab=NULL, zlab=NULL, 
	col=NULL, col.data=2, props=c(75,50,25), method="contour", 
	ticktype = "detailed", ...)

plot.kepdfdd(x, y, ..., main=NULL, method="contour", 
	text.diag.panel = NULL, gap = 0.5)

plot.kepdf(x, y, ...)
}
\arguments{
  \item{x}{
An object of \code{\link{kepdf-class}};
}
  \item{y}{
Not used; for compatibility with generic plot;
}
  \item{eval.points}{
A matrix of data points at which the density to be plotted has to be
evaluated; the number of columns must correspond to the dimension of the
sample data used to estimate the density. If not provided, density is
evaluated on a grid defined on the range of sample data.}
  \item{n.grid}{
A vector with length set to the number of column of sample data, defining the
length of the grid on which the density to be plotted is evaluated; this
argument is ignored when eval.points is not \code{NULL}.
  }
\item{data}{Data to be add to the plot.}
\item{add}{Logical. If \code{TRUE}, add to a current plot.}
\item{main}{An overall title for the plot}
\item{xlab}{A title for the x axis}
\item{ylab}{A title for the y axis}
\item{zlab}{A title for the z axis}
\item{col}{A specification for the plotting color}
\item{col.data}{A specification for the color of \code{data}. Ignored if \code{data} is \code{NULL}.}
\item{type}{What type of plot should be drawn. This argument applies when kernel density estimate is performed on unidimensional data only. Default value is "l".}
\item{props}{A vector defining the fraction of the data to be included within each density level. This argument applies when kernel density estimate is performed on
multidimensional data only.}
\item{method}{One of \code{c("contour", "image", "perspective")}. To be used when two or higher dimensional data 
have been used to estimate the density.}
\item{ticktype}{Character: "simple" draws just an arrow parallel to the axis to indicate direction of increase;
 "detailed" draws normal ticks; to be used if \code{method="perspective"} only.}
\item{text.diag.panel}{Text to be displayed on the diagonal panels when plotting densities estimated on higher than 2-dimensional data.}
\item{gap}{Distance between subplots, when plotting densities estimated of 2-dimensional data or higher-dimensional data.}
  \item{\dots}{
Further arguments to be passed to \code{plot}, \code{image}, \code{contour}, \code{persp}.
}
}
\details{
When density estimates is based on two or higher dimensional data, these functions make use of functions \code{\link[graphics:contour]{contour}},
\code{\link[graphics:image]{image}} and \code{\link[graphics:persp]{persp}}. 
% from package \code{\link[graphics:graphics-package]{graphics-package}}.
For densities estimated on higher than 2-d data, the paiwaise marginal estimated densities are plotted for alle the possible couple of
coordinates.
}
\value{
A list containing the following elements:
  \item{eval.points }{data points at which the plotted density has been evaluated}
  \item{estimate }{the estimated density at \code{eval.points}}
}

\section{Methods}{
\describe{
\item{\code{signature(x = "kepdf", y = "missing")}}{

S4 method for plotting objects of \code{\link{kepdf-class}}}
}
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
\code{\link{kepdf-class}}, \code{\link[graphics:plot]{plot}}, 
\code{\link[graphics:contour]{contour}}, \code{\link[graphics:image]{image}},
 \code{\link[graphics:persp]{persp}}
}
\examples{
#1-d example
set.seed(123)
x1 <- rnorm(50)

#normal optimal bandwidth
pdf1a <- kepdf(x1)
#shrink the smoothing parameter
pdf1b <- kepdf(x1, h=0.5*h.norm(x1))

plot(pdf1a, n.grid=50, data=x1, xlab="x1", ylim=c(0, max(c(pdf1a@estimate, pdf1b@estimate))))
plot(pdf1b, n.grid=50, lty=2, add=TRUE)

#2-d example
set.seed(123)
x2 <- cbind(rnorm(50),rnorm(50))

pdf2 <- kepdf(x2)

plot(pdf2, n.grid=c(50,50), data=x2)
plot(pdf2, n.grid=c(50,50), method="image")
plot(pdf2, n.grid=c(50,50), method="perspective", phi=30, theta=30)

#3-d example
set.seed(123)
x3 <- cbind(rnorm(50), rnorm(50), rnorm(50))

pdf3 <- kepdf(x3)

plot(pdf3, n.grid=c(50,50,50))
plot(pdf3, n.grid=c(50,50,50), method="image", col = terrain.colors(30))
plot(pdf3, n.grid=c(50,50,50), method="perspective", phi=30, theta=30)
}
\keyword{hplot} 

