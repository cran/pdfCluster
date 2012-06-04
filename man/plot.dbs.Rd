\name{plot,dbs-method}
\docType{methods}
\alias{plot,dbs-method}
\alias{plot.dbs}
\alias{plot,dbs,missing-method}

\title{
Plot objects of class dbs 
}
\description{
This function provides a graphical tool to display diagnostics 
of density-based cluster analysis by means of the density-based silhouette information.
}
\usage{
\S4method{plot}{dbs}(x, y , xlab = "", ylab = "", col = NULL, lwd = 3, cex = 0.9, 
    cex.axis = 0.5, main = NULL, labels = FALSE, ...)
}

\arguments{
  \item{x}{
An object of \code{\link{dbs-class};}
}
  \item{y}{
Not used; for compatibility with generic plot;
}
\item{xlab}{A title for the x axis;}
\item{ylab}{A title for the y axis;}
  \item{col}{
A specification for the plotting color. Default are colors in palette corresponding to the group labels;}
  \item{lwd}{
A specification for the width of the bars in the plot;}
  \item{cex}{
A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default; 
  }
  \item{cex.axis}{
The magnification to be used for axis annotation relative to the current setting of cex;}
  \item{main}{
An overall title for the plot; }
  \item{labels}{
Logical. Should row index of data be added to the plot? }
  \item{\dots}{
Further arguments to be passed to \code{plot};}
}
\details{
After computing the density-based silhouette index by applying \code{\link{dbs-methods}}, data are partitioned 
into the clusters, sorted in a decreasing order with respect to their dbs value and displayed
on a bar graph.}
%%\value{
%%Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}

\section{Methods}{
\describe{
\item{\code{signature(x = "dbs", y = "missing")}}{

S4 method for plotting objects of \code{\link{dbs-class}}}
}
}
\references{
Menardi, G. (2010) Density-based Silhouette diagnostics for clustering methods. \emph{Statistics and Computing}, 21, 295-308.}

%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dbs}}, \code{\link{dbs-class}}, \code{\link[cluster:silhouette]{silhouette}}.
}
\examples{
#example 1: no groups in data
#random generation of group labels
set.seed(54321)
x <- rnorm(50)
groups <- sample(1:2, 50, replace=TRUE)
groups
dsil <- dbs(x=as.matrix(x), clusters=groups)
dsil
summary(dsil)
plot(dsil, labels=TRUE, lwd=6)

#example 2: wines data
# load data
data(wine)

gr <- wine[,1]

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
 
dsil <- dbs(cl)
plot(dsil)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
 
