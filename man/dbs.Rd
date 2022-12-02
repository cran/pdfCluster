\name{dbs}
\docType{methods}
\alias{dbs}
%\alias{dbs.cluster}
%\alias{dbs.pdfCluster}
\alias{dbs-methods}
\alias{dbs,matrix,numeric-method}
\alias{dbs,pdfCluster,missing-method}
\alias{dbs,matrix-method}
\alias{dbs,pdfCluster-method}


\title{Density-based silhouette information methods}

\description{
Computes the density-based silhouette information of clustered data. Two methods are associated to this function. The first
method applies to two arguments: the matrix of data and the vector of cluster labels; the second method applies to objects of 
\code{\link{pdfCluster-class}}.
}

\usage{
\S4method{dbs}{matrix}(x, clusters, h.funct="h.norm", hmult=1, prior, ...)

\S4method{dbs}{pdfCluster}(x, h.funct="h.norm", hmult = 1, prior = 
   as.vector(table(x@cluster.cores)/sum(table(x@cluster.cores))), 
   stage=NULL, ...)
}

\arguments{
  \item{x}{
A matrix of data points partitioned by any density-based clustering method or an object of \code{\link{pdfCluster-class}}.}

  \item{clusters}{
Cluster labels of grouped data. This argument has not to be set when \code{x} is a \code{\link{pdfCluster-class}} object.
}

  \item{h.funct}{
Function to estimate the smoothing parameters. Default is \code{\link{h.norm}}. 
}

  \item{hmult}{
Shrink factor to be multiplied by the smoothing parameters. Default value is 1. 
}

  \item{prior}{
Vector of prior probabilities of belonging to the groups. When \code{x} is of \code{\link{pdfCluster-class}},
  default value is set proportional to the cluster cores cardinalities. Otherwise, equal prior probabilities are 
  given to the clusters by default. 
 }

  \item{stage}{
When \code{x} is a \code{\link{pdfCluster-class}} object, this is the stage of classification of low-density data at which the dbs has to be computed.
Default value is the number of stages of the procedure. Set it to 0 if the dbs has to be computed at cluster cores only.}
 
  \item{\dots}{
Further arguments to be passed to methods (see \code{\link{dbs-methods}}) or arguments to \code{\link{kepdf}}. See details below.}
}

\details{
This function provides diagnostics for a clustering produced by any density-based clustering method. The \code{dbs}
information is a suitable modification of the \code{\link[cluster:silhouette]{silhouette}} information aimed at evaluating 
the cluster quality in a density based framework. It is based on the estimation of data posterior probabilities of belonging to the clusters. It may be 	 
used to measure the quality of data allocation to the clusters. High values of the \eqn{\hat{dbs}} are evidence of a good quality clustering. 	 

Define
\deqn{
\hat{\tau}_m(x_i)=\frac{\pi_{m} \hat{f}(x_i|x_ \in m)}{\sum_{m=1}^M \pi_{m}\hat{f}(x_i|x_i \in m)} \quad m=1,\ldots,M, 
}

where \eqn{\pi_{m}} is a prior probability of \eqn{m} and \eqn{\hat{f}(x_i|x_i \in m)} is a density estimate at \eqn{x_i} evaluated with function \code{\link{kepdf}} by using the only data points in \eqn{m}. Density estimation is performed with fixed bandwidths \code{h}, as evaluated by function \code{h.funct}, possibly multiplied by the shrink factor \code{hmult}.

Density-based silhouette information of \eqn{x_i}, the \eqn{i^{th}} row of the data matrix \code{x}, is defined as follows:
\deqn{
\hat{dbs}_i=\frac{\log\left(\frac{\hat{\tau}_{m_{0}}(x_i)}{\hat{\tau}_{m_{1}}(x_i)}\right)}{{\max}_{x_i }\left| \log\left(\frac{\hat{\tau}_{m_{0}}(x_i)}{\hat{\tau}_{m_{1}}(x_i)}\right)\right|},
}
where \eqn{m_0} is the group where \eqn{x_i} has been allocated and \eqn{m_1} is the group for which \eqn{\tau_m} is maximum, \eqn{m\neq m_0}. 

\emph{Note:} when there exists \eqn{x_j} such that \eqn{\hat{\tau}_{m_{1}}(x_j)} is zero, \eqn{\hat{dbs}_j} is forced to 1 and \eqn{{\max}_{x_i }\left| \log\left(\frac{\hat{\tau}_{m_{0}}(x_i)}{\hat{\tau}_{m_{1}}(x_i)}\right)\right|}
is computed by excluding \eqn{x_j} from the data matrix \code{x}.

See Menardi (2011) for a detailed treatment.}


\value{
An object of class \code{"dbs"}, with slots:  
  \item{call }{The matched call.}
  \item{x }{The matrix of clustered data points.}
  \item{prior }{The vector of prior probabilities of belonging to the groups.}
  \item{dbs }{A vector reporting the density-based silhouette information of the clustered data.}
  \item{clusters }{Cluster labels of grouped data.}
  \item{noc }{Number of clusters}
  \item{stage }{If argument \code{x} of \code{dbs} is a  \code{\link{pdfCluster-class}} object, this slot provides the stage of the classification at which the dbs is computed.}

  See \code{\link{dbs-class}} for more details.
}


\section{Methods}{
\describe{
\item{\code{signature(x = "matrix", clusters = "numeric")}}{

Computes the density based silhouette information for objects partitioned 
according to any density-based clustering method.
}
\item{\code{signature(x = "pdfCluster", clusters = "missing")}}{

Computes the density based silhouette information for objects of class 
\code{"pdfCluster"}}.
}}

\references{
Menardi, G. (2011) Density-based Silhouette diagnostics for clustering methods. \emph{Statistics and Computing}, 21, 295-308.}


\seealso{
\code{\link{dbs-class}}, \code{\link{plot,dbs-method}}, \code{\link[cluster:silhouette]{silhouette}}.
}
\examples{
#example 1: no groups in data
#random generation of group labels
set.seed(54321)
x <- rnorm(50)
groups <- sample(1:2, 50, replace = TRUE)
groups
dsil <- dbs(x = as.matrix(x), clusters=groups)
dsil
summary(dsil)
plot(dsil, labels=TRUE, lwd=6)

#example 2: wines data
# load data
data(wine)

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
 
dsil <- dbs(cl)
plot(dsil)

}

\keyword{cluster}
