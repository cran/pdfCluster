\name{pdfClassification}
\alias{pdfClassification}

\title{
Classification of low density data
}
\description{
Allocates low density data points in a multi-stage procedure after that cluster cores have been detected by applying \code{\link{pdfCluster}}}
\usage{
pdfClassification(obj, h.funct = "h.norm", n.stage = 5)
}

\arguments{
  \item{obj}{
An object of \code{\link{pdfCluster-class}} obtained through a call to \code{\link{pdfCluster}} by setting 
the \code{n.stage} argument greater than 0.}
  \item{h.funct}{
Function to be used to evaluate the smoothing parameters in order to estimate the probability density function conditioned 
to the group membership. Default value is \code{\link{h.norm}}. If missing, the smoothing parameter used to detect the cluster 
cores are used.}
  \item{n.stage}{
Allocation of low density data is performed by following a multi-stages procedure in \code{n.stage} stages. 
}
}
\details{
The basic idea of the classification stage of the procedure is as follows: for an unallocated data point \eqn{x_0},
compute the estimated density \eqn{\hat{f}_m(x_0)} based on the data already assigned to group \eqn{m, m = 1, 2, \ldots, M},
 and assign \eqn{x_0} to the group with highest ratio \eqn{\frac{\hat{f}_m(x_0)}{max_m \hat{f}_m(x_0)}}. 
In case \eqn{\hat{f}_m(x_0)}=0, \eqn{m, m = 1, 2, \ldots, M}, \eqn{x_0} is considered as an outlier. The procedure gives a warning 
message and the outlier remains unclassified. The cluster label of \eqn{x_0} will be set to zero.   
 The current implementation of this idea proceeds in a block-sequential manner, allocating several points at a time,
 updating the estimates \eqn{\hat{f}_m(\cdot)} based on the new members of each group and then allocating a new block of points. 
 }
\value{
An object of \code{\link{pdfCluster-class}} with slot \code{stages} non \code{NULL} and of class \code{"list"}. 
See \code{\link{pdfCluster-class}} for further details.  
}
\references{
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation. \emph{Statistics and Computing}. 17, 71-80.
}
%%\author{
%%  ~~who you are~~
%%}
\note{ 
Function \code{\link{pdfClassification}} may be called both internally, from \code{\link{pdfCluster}}, by setting the argument
\code{n.stage} greater than zero, or externally, by providing as argument an object of \code{\link{pdfCluster-class}} whose slot
\code{stages} is \code{NULL}.

When \code{\link{pdfClassification}} is internally called from \code{\link{pdfCluster}} and one group only is detected, 
the slot \code{stages} is a list with \code{n.stage} elements, each of them being a vector with length set to the number of data and set all to 1.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfCluster}}, \code{\link{pdfCluster-class}}, 
\code{\link{summary.pdfCluster-class}}  
}
\examples{
# load data
data(wine)

# select a subset of variables
x <- wine[, c(2,5,8)]

#whole procedure, included the classification phase
cl <- pdfCluster(x)
summary(cl)
table(cl@clusters)

#cluster cores only are detected, classification is performed separately 
#same output as above
cl0 <- pdfCluster(x, n.stage=0)
summary(cl0)
table(cl0@cluster.cores)

cl1 <- pdfClassification(cl0)
table(cl1@clusters)

#comparison 
table(cl1@clusters,cl@clusters)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{smooth}

