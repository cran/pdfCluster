\name{pdfCluster}
\alias{pdfCluster}
\title{
Clustering via kernel density estimation
}
\description{
Cluster analysis is performed by the density-based procedure described in Azzalini and Torelli (2007). 
}
\usage{
pdfCluster(x, h = h.norm(x), hmult = 0.75,
     n.grid = min(50, nrow(as.matrix(x))), n.stage = 5)
}
\arguments{
  \item{x}{
A vector, a matrix or a data frame of numeric data to be partitioned.  
}
  \item{h}{
A vector with length equal to \code{ncol(x)} (or 1, if \code{x} is a
vector), defining the smoothing parameters to be used 
to estimate density of data by kernel methods. An asymptotically
 normal optimal smoothing parameter is used as a default argument. 
}
  \item{hmult}{
Shrink factor with length 1, applied to \code{h}. Default value is 0.75.
}
  
  \item{n.grid}{
Defines the length of the grid (plus 1) on which evaluating the connected components of the level sets, that is
the number of sections of the density estimate to be considered to find the clusters. The default value is 
set to the minimum between 50 points and the number of rows data.
  }
  
  \item{n.stage}{
Allocation of low density data is performed by following a multi-stages procedure in \code{n.stage} steps through an internal call
to \code{\link{pdfClassification}}. Default value is set to 5. Set the parameter to zero to detect the cluster cores only. See details below.
  }
}

\details{
Clusters are defined as high density regions. Any section of the density \eqn{f} underlying data, at level \eqn{k}, induces a partition of the sample space
 into two sets, one having density up to \eqn{k}, one having density greater than \eqn{k}.  The clusters corresponds to the maximum connected regions of the latter set. 
 As \eqn{k} varies, these clusters may be represented according to a hierarchical structure in the form of a tree. 
 Density estimation is here performed by kernel methods and allows for the detection of high density regions. The associated connected regions are 
 approximated by the polyhedrons formed by applying a Delaunay triangulation on data with density above \eqn{k.}. This procedure creates M groups
 of points, namely the "cluster cores", and it leaves a number of points unlabeled. 
 The allocation of lower density data is subsequently performed by following a multi-stages supervised classification-like approach.  
  }
\value{
An S4 object of \code{\link{pdfCluster-class}} with slots:
  \item{call }{The matched call.}
  \item{x }{The matrix of data input (if a vector of data is provided as input, a one-column matrix is returned as output).}
  \item{estimate }{Vector of density estimate evaluated at the data points.}
  \item{h }{Vector of smoothing parameters used to estimate the probability density function.}
  \item{hmult }{Shrink factor multiplied by the smoothing parameters.}
  \item{nc }{An object of class list defining details about the identification of the connected regions. It includes:
  \itemize{
   \item \code{nc} Number of connected for each value of section of the density estimate, corresponding to the values of \code{p}.
   \item \code{p}  Vector of proportions of data with estimated density below a threshold.
   \item \code{id} Group label at each point for different sections of the density estimate, corresponding to the values of \code{p}. 
					Negative values of \code{id} means that the estimated density is below the considered threshold.
   \item \code{pq} For each \code{p} gives the corresponding quantile \code{q} of the values of the density
  }}
  \item{cluster.cores }{Vector with the same length as the length (or number of rows) of \code{x}, defining the cluster cores membership. 
						\code{NA} values correspond to low density, unlabeled data, to be classified in the second phase of the procedure by calling \code{\link{pdfClassification}};}
  \item{tree }{Hierarchical structure in the form of a cluster tree with leaves corresponding to the connected components associated to different sections of the density estimate.
  The object has class \code{\link[stats:dendrogram]{dendrogram}}.}             
  \item{noc }{number of clusters.}
  \item{stages }{List with elements corresponding to the data allocation to groups at the different stages of the classification procedure.
   \code{NA} values correspond to unlabeled data.}
  \item{clusters }{set to \code{NULL} if the n.stages=0, that is, if data belonging to the cluster cores only have been allocated. 
  Otherwise it reports the final label groups.}	
  }
  
\references{
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation. \emph{Statistics and Computing}. 17, 71-80.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfCluster-class}}, \code{\link{pdfClassification}},
\code{\link{summary.pdfCluster-class}}  
}

\examples{
##########
#example 1
# load data
data(wine)
gr <- wine[,1]

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
summary(cl)
plot(cl)

#comparison with original groups
table(cl@clusters,gr)

###########
# The code below essentially reproduces the example in Section 4.3 of
# Azzalini and Torelli (2007), up to a minor difference, presumably due
# to a slightly different choice of the tuning parameters, whose trace
# has been lost.

#loading data
data(oliveoil)

#preparing data
olive1 <- 1 + oliveoil[, 3:10]
margin <- apply(data.matrix(olive1),1,sum)
olive1 <- olive1/margin
alr <- (-log( olive1[, -4]/olive1[, 4]))
#select the first 5 principal components
x <- princomp(alr, cor=TRUE)$scores[, 1:5]

#estimation of density
pdf <- kepdf(x)
summary(pdf)
plot(pdf, gap=0.25)

#clustering
cl <- pdfCluster(x)
summary(cl)
plot(cl)

#comparing groups with original macro-area membership
groups <- cl@clusters
table(groups, oliveoil$macro.area)

#cluster cores
table(cl@cluster.cores)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{multivariate}
\keyword{smooth}
\keyword{smooth}
