\name{pdfClassification}
\alias{pdfClassification}

\title{
Classification of low density data
}
\description{
Allocates low density data points in a multi-stage procedure after that cluster cores have been detected 
by applying \code{\link{pdfCluster}}.}

\usage{
pdfClassification(obj, n.stage = 5, se = TRUE, hcores = FALSE)
}

\arguments{
  \item{obj}{
An object of \code{\link{pdfCluster-class}}.}
  \item{n.stage}{
Allocation of low density data is performed by following a multi-stages procedure in \code{n.stage} stages. 
}
  \item{se}{
Logical. Should the standard-error of the density estimates be taken into account to define the order of allocation? Default value is TRUE. 
See details below.
  }
  \item{hcores}{
Logical. Set this value to TRUE to build cluster density estimates by selecting the same bandwidths as the ones used to form 
the cluster cores. Otherwise, bandwidths specific for the clusters are selected. Default value is FALSE. See details below.   

  }
}
\details{
The basic idea of the classification stage of the procedure is as follows: for an unallocated data point \eqn{x_0},
compute the estimated density \eqn{\hat{f}_m(x_0)} based on the data already assigned to group \eqn{m, m = 1, 2, \ldots, M},
 and assign \eqn{x_0} to the group with highest log ratio \eqn{\hat{f}_m(x_0)/\max_m \hat{f}_m(x_0)}. 

 In case \eqn{\hat{f}_m(x_0)}=0, for all \eqn{m = 1, 2, \ldots, M}, \eqn{x_0} is considered as an outlier. The procedure gives a warning 
message and the outlier remains unclassified. The cluster label of \eqn{x_0} will be set to zero.   
 
The current implementation of this idea proceeds in \code{n.stage} stages, allocating a block of points at a time,
updating the estimates \eqn{\hat{f}_m(\cdot)} based on the new members of each group and then allocating a new block of points. 
When \code{se = TRUE}, classification is performed by further weighting the log-ratios inversely with their approximated standard
error, so that points whose density estimate has highest precision are allocated first.

Each of the \eqn{\hat{f}_m(\cdot)} is built by selecting either the same bandwidths \eqn{h_0} as the ones used to form the cluster cores (when \code{hcores = TRUE}) or cluster-specific bandwidths, obtained as follows: 
\deqn{h_m^{*} = \exp [(1-a_m) \log(h_0) + a_m \log(h_m)],}
where \eqn{a_m} is the proportion of data points in the \eqn{m}-th cluster core and \eqn{h_m} are asymptotically optimal for a normal distribution of the \eqn{m}-th cluster or computed according to the Silverman (1986) approach, if the kernel estimator has fixed or adaptive bandwidth, respectively. 
 }
\value{
An object of \code{\link{pdfCluster-class}} with slot \code{stages} of class \code{"list"} having length equal to \code{n.stage}. 
See \code{\link{pdfCluster-class}} for further details.  
}
\references{
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation. \emph{Statistics and Computing}. 17, 71-80.

Silverman, B. (1986). \emph{Density estimation for statistics and data analysis}. Chapman and Hall, London.
}
%%\author{
%%  ~~who you are~~
%%}
\note{ 
Function \code{\link{pdfClassification}} is called internally, from \code{\link{pdfCluster}}, when the argument
\code{n.stage} is set to a value greater than zero. Alternatively, it may be called externally, by providing as 
argument an object of \code{\link{pdfCluster-class}}.

When \code{\link{pdfClassification}} is internally called from \code{\link{pdfCluster}} and one group only is detected, 
the slot \code{stages} is a list with \code{n.stage} elements, each of them being a vector with length equal to the number
of data points and all elements equal to 1.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfCluster}}, \code{\link{pdfCluster-class}}
}

\examples{
# load data
data(wine)

# select a subset of variables
x <- wine[, c(2,5,8)]

#whole procedure, included the classification phase
cl <- pdfCluster(x)
summary(cl)
table(groups(cl))

#use of bandwidths specific for the group 
cl1 <- pdfClassification(cl, hcores= TRUE)
table(groups(cl1))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

