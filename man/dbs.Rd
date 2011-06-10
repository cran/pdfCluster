\name{dbs}
\docType{methods}
\alias{dbs.cluster}
\alias{dbs}
\alias{dbs.pdfCluster}
\alias{dbs-methods}
\alias{dbs,matrix,numeric-method}
\alias{dbs,pdfCluster,missing-method}

\title{Density-based Silhouette information methods}

\description{
Computes the Density-based Silhouette information at clustered data. Two methods are associated to these functions. The first
method applies to two arguments: the matrix of data and the vector of cluster labels; the secondo method applies to objiect of 
\code{\link{pdfCluster-class}}.
}

\usage{
dbs.cluster(x, clusters, h.funct="h.norm", hmult=0.75, prior, ...)
dbs.pdfCluster(x, h.funct="h.norm", hmult=0.75, prior, stage=NULL, ...)
}

\arguments{
  \item{x}{
A matrix of data points partitioned by any density-based clustering method or an object of \code{\link{pdfCluster-class}}}

  \item{clusters}{
Cluster labels of grouped data. This argument is missing when \code{x} is a \code{\link{pdfCluster-class}} object.
}

  \item{h.funct}{
Function to estimate the smoothing matrix. Default is \code{\link{h.norm}}. 
}

  \item{hmult}{
Shrink factor to be multiplied by the smoothing parameters. Default value is 1. 
}

  \item{prior}{
Vector of prior probabilities of belonging to the groups. If missing, a uniform distribution of the prior probabilities over the groups is chosen.  
 }

  \item{stage}{
When \code{x} is a \code{\link{pdfCluster-class}} object, this is the stage of classification of low-density data at which the dbs has to be computed.
Default value is the total number of stages of the procedure. Set it to 0 if the dbs has to be computed at cluster cores only.}
 
  \item{\dots}{
Further arguments to be passed to methods (see \code{\link{dbs-methods}}) or arguments to \code{\link{kepdf}}. See details below.}
}

\details{
This function provides diagnostics on data clustered according to any density-based clustering method. The \code{dbs}
information is a suitable modification of the \code{\link[cluster:silhouette]{silhouette}} information aimed at evaluating 
the clusters quality in a density based framework. It is based on the estimation of data posterior probabilities of belonging to the clusters. It may be 	 
used to measure the quality of data allocation to the clusters as well as to choose the best partition among different ones. 	 

Density-based Silhouette information of \eqn{x_i}, \eqn{i^{th}} row of the data matrix \code{x} is defined as follows:

\deqn{
\hat{dbs}_i=\frac{\log\left(\frac{\hat{\tau}_{m_{0}}(x_i)}{\hat{\tau}_{m_{1}}(x_i)}\right)}{{\max}_{x_i }\left| \log\left(\frac{\hat{\tau}_{m_{0}}(x_i)}{\hat{\tau}_{m_{1}}(x_i)}\right)\right|},
}
where \eqn{m_0} is the group where \eqn{x_i} has been allocated and \eqn{m_1} is the group for which \eqn{\tau_m} is maximum, \eqn{m\neq m_0}. \eqn{\tau_m} 
is estimated by
\deqn{
\hat{\tau}_m(x_i)=\frac{\pi_{m} \hat{f}(x_i|x_ \in m)}{\sum_{m=1}^M \pi_{m}\hat{f}(x_i|x_i \in m)} \quad m=1,\ldots,M 
}

\eqn{\hat{f}(x_i|x_i \in m)} is a density estimate at \eqn{x_i} evaluated with function \code{\link{kepdf}} by using the only data points in \eqn{m},
and \eqn{\pi_{m}} is a prior probability of \eqn{m}.

\emph{Note:} when there exists \eqn{x_j} such that \eqn{\hat{\tau}_{m_{1}}(x_j)} is zero, \eqn{\hat{dbs}_j} is forced to 1 and \eqn{{\max}_{x_i }\left| \log\left(\frac{\hat{\tau}_{m_{0}}(x_i)}{\hat{\tau}_{m_{1}}(x_i)}\right)\right|}
is taken by excluding \eqn{x_j} from tha data matrix \code{x}.}

\value{
An object of class \code{"pdfSilhouette"}, with slots:  
  \item{x }{The matrix of clustered data points.}
  \item{prior }{The vector of prior probabilities of belonging to the groups.}
  \item{dbs }{A vector reporting the density-based Silhouette information of the clustered data.}
  \item{clusters }{Cluster labels of grouped data.}
  \item{nc }{Number of clusters}
  \item{stage }{If argument \code{x} of \code{dbs} is a  \code{\link{pdfCluster-class}} object, this slot provides the stage of the classification at which the dbs is computed.}
  \item{call }{The matched call.}

  See \code{\link{pdfSilhouette-class}} for more details.
}


\section{Methods}{
\describe{

\item{\code{signature(x = "matrix", clusters = "numeric")}}{

Computes the density based Silhouette information for objects partitioned according to any density-based clustering method.
}

\item{\code{signature(x = "pdfCluster", clusters = "missing")}}{

Computes the density based Silhouette information for objects of class \code{"pdfCluster"}
}
}}

\references{
Menardi, G. (2010) Density-based Silhouette diagnostics for clustering methods. \emph{Statistics and Computing}, DOI: 10.1007/s11222-010-9169-0.}

%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfSilhouette-class}}, \code{\link{plot.pdfSilhouette}}, \code{\link[cluster:silhouette]{silhouette}}
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

\keyword{cluster}
\keyword{smooth}