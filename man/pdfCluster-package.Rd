\name{pdfCluster-package}
\alias{pdfCluster-package}
\docType{package}
\title{The pdfCluster package: summary information}

\description{
This package performs cluster analysis via kernel density estimation (Azzalini 
and Torelli, 2007; Menardi and Azzalini, 2014). 
Clusters are associated to the maximally connected components with estimated 
density above a threshold. As the threshold varies, these clusters may be 
represented according to a hierarchical structure in the form of a tree. 
Detection of the connected regions is conducted by means of the Delaunay 
tesselation when data dimensionality is low to moderate, following Azzalini
and Torelli (2007). For higher dimensional data, detection of connected regions
is performed according to the procedure described in Menardi and Azzalini (2013).
In both cases, after that a number of high-density cluster-cores is identified, 
lower density data are allocated by following a supervised classification-like
approach. The number of clusters, corresponding to the number of the modes of the 
estimated density, is automatically selected by the procedure. 
Diagnostics methods for evaluating the quality of clustering are also available
(Menardi, 2011). Moreover, the package provides a routine to estimate the 
probability density function by kernel methods, given a set of data with 
arbitrary dimension. 
The main features of the package are described and illustrated
in Azzalini and Menardi (2014).
}
\details{
 The \code{\link{pdfCluster-package}} makes use of classes and methods of the 
 S4 system.
 It includes some foreign functions written in the C language: two of them
 compute the kernel density estimate of data and are interfaced by the R 
 function \code{\link{kepdf}}. 
 Other C routines included in the package allow for a quicker detection of the
 connected components of the subgraphs associated with the level sets of the data.
 Two of them are directly drawn from the homonymous ones in the \code{\link[spdep]{spdep}} 
 package. 
 
 Starting from version 1.0-0, new features have been introduced: 
 \itemize{
	\item{kernel density estimation may be performed by using either a
	a fixed or an adaptive bandwidth; moreover, the option of selecting 
	a Student's \eqn{t} kernel has been included, for computational convenience;}
	\item{detection of connected components of the level sets is performed
	by means of the Delaunay triangulation when data dimensionality is 
	up to 6, following Azzalini and Torelli (2007); 
        for higher dimensional data a new procedure, which is less time-consuming,
        is now adopted (Menardi and Azzalini, 2014);}
    \item{the order of classification of lower density data depends now also on
	the standard error of the estimated density ratios; moreover, a cluster-specific 
	bandwidth is the default option to classify low density data.}
     }
  See examples below to understand how to set arguments of the main function
  of the package, in order to obtain the same results as the ones obtained with
  versions 0.1-x.
  }
\author{
Adelchi Azzalini, Giovanna Menardi, Tiziana Rosolin

Maintainer: Giovanna Menardi <menardi at stat.unipd.it>
}
\references{
Azzalini, A., Menardi, G. (2014). Clustering via Nonparametric Density Estimation: The R Package pdfCluster.
        \emph{Journal of Statistical Software}, 57(11), 1-26,
        URL http://www.jstatsoft.org/v57/i11/.
        
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation.
 \emph{Statistics and Computing}, 17, 71-80.

Menardi G. (2011). Density based Silhouette diagnostics for clustering methods.
 \emph{Statistics and Computing}, 21, 295-308.

Menardi G., Azzalini, A. (2014). An advancement in clustering via nonparametric density estimation.
 \emph{Statistics and Computing}, DOI: 10.1007/s11222-013-9400-x, to appear.
}
\keyword{ package }

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
table(groups(cl),gr)

#density based silhouette diagnostics
dsil <- dbs(cl)
plot(dsil)

##########
# higher dimensions

x <- wine[, -1]

#density estimation with adaptive bandwidth 
pdf <- kepdf(x, bwtype="adaptive")
summary(pdf)
#density plot is not much clear for high- dimensional data
#select a few variables
plot(pdf, indcol = c(1,4,7))

#clustering
#when dimension is >= 6, default method to find connected components is "pairs"
#density is better estimated by using an adaptive bandwidth
cl <- pdfCluster(x, bwtype="adaptive")
summary(cl)
plot(cl)

########
# this example shows how to set the arguments in function pdfCluster
# in order to obtain the same results as the ones of versions 0.1-x.
x <- wine[, c(2, 5, 8)]

# previous versions of the package 
# do not run
# old code: 
# cl <- pdfCluster(x)

# same result is obtained now obtained as follows:
cl <- pdfCluster(x, se=FALSE, hcores= TRUE, graphtype="delaunay", n.grid=50)
}

\keyword{cluster}
\keyword{multivariate}
\keyword{nonparametric}


