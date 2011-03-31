\name{pdfCluster-package}
\alias{pdfCluster-package}
\docType{package}
\title{The pdfCluster package: summary information}

\description{
This package performs cluster analysis via kernel density estimation (Azzalini and Torelli, 2007). 
Clusters are associated to the maximally connected components with estimated density above a threshold.
As the threshold varies, these clusters may be represented according to a hierarchical structure in the form of a tree. 
Detection of the connected regions is conducted by means of the Delaunay tesselation. The allocation of lower density
data is performed by following a supervised classification-like approach. 
The number of clusters, corresponding to the number of the modes of the estimated density, is automatically detected by the procedure. 
Diagnostics methods for evaluating the quality of clustering are also available (Menardi, 2010).
Moreover, the package provides a routine to estimate the probability density function by kernel methods, given a set of data 
with arbitrary dimension. 
}
\details{
%%\tabular{ll}{
%%Package: \tab pdfCluster\cr
%%Type: \tab Package\cr
%%Version: \tab 1.0\cr
%%Date: \tab 2011-02-11\cr
%%License: \tab GPL-2\cr
%%LazyLoad: \tab yes\cr
%%Depends: \tab geometry, methods\cr 
%%}
 The \code{\link{pdfCluster-package}} is written in S4 and, in its current form, is a preliminary version. 
 It includes a foreign function written in the Fortran-77 language which computes the kernel density 
 estimate of data. This function is interfaced by the R function \code{\link{kepdf}}. 
 Two other C routines are included in the package and allow for a quicker detection of the connected 
 high density regions in multivariate spaces. These functions are directly drawn from the homonymous ones
 in the \code{\link[spdep]{spdep}} package. 
}
\author{
Adelchi Azzalini, Giovanna Menardi, Tiziana Rosolin

Maintainer: Giovanna Menardi <menardi at stat.unipd.it>
}
\references{
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation. \emph{Statistics and Computing}. 17, 71-80.

Menardi G. (2010). Density based Silhouette diagnostics for clustering methods. \emph{Statistics and Computing}. DOI: 10.1007/s11222-010-9169-0.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
# load data
data(wine)
gr <- wine[, 1]

# select a subset of variables
x <- wine[, c(2, 5, 8)]

#density estimation
pdf <- kepdf(x)
summary(pdf)
plot(pdf)

#clustering
cl <- pdfCluster(x)
summary(cl)
plot(cl)

#comparison with original groups
table(cl@clusters,gr)

#density based silhouette diagnostics
dsil <- dbs(cl)
plot(dsil)

}

\keyword{cluster}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{smooth}

