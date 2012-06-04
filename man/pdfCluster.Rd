\name{pdfCluster}
\docType{methods}
\alias{pdfCluster}
%\alias{pdfCluster.data}
%\alias{pdfCluster.graphpairs}
\alias{pdfCluster-methods}
\alias{pdfCluster,data.frame,ANY,ANY,ANY,ANY-method}
\alias{pdfCluster,matrix,ANY,ANY,ANY,ANY-method}
\alias{pdfCluster,numeric,ANY,ANY,ANY,ANY-method}
\alias{pdfCluster,pdfCluster,missing,missing,ANY,missing,ANY-method}
\alias{pdfCluster,pdfCluster,missing,missing,ANY,missing-method}
\alias{pdfCluster,data.frame-method}
\alias{pdfCluster,matrix-method}
\alias{pdfCluster,numeric-method}
\alias{pdfCluster,pdfCluster-method}



\title{ Clustering via nonparametric density estimation}

\description{
Cluster analysis is performed by the density-based procedures described in 
Azzalini and Torelli (2007) and Menardi and Azzalini (2012). Operationally, the 
kernel estimation method is used.
}

\usage{
\S4method{pdfCluster}{numeric}(x, graphtype, Q="QJ", lambda = 0.10, grid.pairs = 10, 
      n.grid=min(50, nrow(as.matrix(x))), ... ) 

\S4method{pdfCluster}{matrix}(x, graphtype, Q="QJ", lambda = 0.10, grid.pairs = 10, 
      n.grid=min(50, nrow(as.matrix(x))), ... ) 

\S4method{pdfCluster}{data.frame}(x, graphtype, Q="QJ", lambda = 0.10, grid.pairs = 10, 
      n.grid=min(50, nrow(as.matrix(x))), ... ) 
	  
\S4method{pdfCluster}{pdfCluster,missing,missing,ANY,missing,ANY}(x, graphtype, Q, 
     lambda, grid.pairs, n.grid = min(50, nrow(as.matrix(x@x))), ...) 
}
\arguments{
  \item{x}{
A vector, a matrix or a data frame of numeric data to be partitioned. 
Alternatively, an object of \code{\link{pdfCluster-class}}
itself, obtained when \code{graphtype} is set to "pairs". See details below.  
}
  \item{graphtype}{
Either "unidimensional", "delaunay" or "pairs", it defines the procedure used 
to build the graph associated with the data. If missing, a "delaunay" graph is 
built for data having dimension less than 7, otherwise a "pairs" graph is built. 
See details below. This argument has not to be set when \code{x} is of
\code{\link{pdfCluster-class}}. 
}
  \item{Q}{
Optional arguments to be given when \code{graphtype = "delaunay"}. See 
\code{\link[geometry:delaunayn]{delaunayn}} in package \code{geometry} 
for further details. This argument has not to be set when \code{x} is of
\code{\link{pdfCluster-class}}.  
}
  \item{lambda}{
Tolerance threshold to be used when \code{graphtype = "pairs"}. An edge is set 
between two observations if the density function, evaluated along the segment
linking them, does not exhibit any valley having a measure exceeding \code{lambda}.
Default value is set to 0.10.
   }
  \item{grid.pairs}{
When \code{graphtype = "pairs"}, this arguments defines the length of the grid of
points along the segment linking each pair of observations, on which the density
is evaluated. Default is 10. This argument has not to be set when \code{x} is of
\code{\link{pdfCluster-class}}. 
}
\item{n.grid}{
Defines the length of the grid (minus one) on which evaluating the connected 
components of the density level sets. The default value is set to the minimum
between 50 points and the number of rows data.
  }
%\item{mult}{
%
%  }
  \item{\dots}{
Further arguments to be passed to \code{\link{kepdf}} or to 
\code{\link{pdfClassification}}. Furthermore, a shrink factor \code{hmult} to be
multiplied by the smoothing parameters may be optionally given. If missing, it is 
taken to be 1 when data have dimension greater than 6, 0.75 otherwise. 
}
}

\details{
%Clusters are defined as high density regions. Any section of the density 
%\eqn{f} underlying data, at level \eqn{k}, induces a partition of the sample
% space into two sets, one having density up to \eqn{k}, one having density
% greater than \eqn{k}. Clusters correspond to the maximum connected regions
% of the latter set. As \eqn{k} varies, these clusters may be represented 
%according to a hierarchical structure in the form of a tree. 
%Density estimation is here performed by kernel methods and allows for the 
%detection of high density regions. M groups of points, namely the "cluster
% cores", are thus created and a number of points are left unlabeled. 
%The allocation of lower density data is subsequently performed by following
% a multi-stages supervised classification-like approach.  
Clusters are associated to the connected components of the level sets of the 
density underlying the data.
Density estimation is performed by kernel methods and the connected regions 
are approximated by the connected components of a graph built on data. 
Univariate data are handled by the "unidimensional" procedure. When data are
multidimensional, two alternative procedures to build the graph are adopted:
\describe{
\item{Delaunay triangulation}{An edge is set between two observations when 
they are contiguous in the Voronoi diagram.}
\item{Pairs evaluation}{An edge is set between two observations when the 
density function, evaluated along the segment joining them, does not exhibit 
any valley having a relative amplitude greater than a tolerance threshold
 \eqn{0 \le \lambda < 1}.} 
}
As the level set varies, the number of detected components gives rise to the
tree of clusters, where each leave corresponds to a mode of the density 
function. Observations around these modes form a number of cluster cores, 
while the lower density observations are allocated according to a 
classification procedure (see also \code{\link{pdfClassification}}). 
}
\value{
An S4 object of \code{\link{pdfCluster-class}} with slots:
  \item{call }{The matched call.}
  \item{x }{The matrix of data input (if a vector of data is provided as input,
  a one-column matrix is returned as output).}
  \item{pdf }{An object of class \code{list} defining details about the 
  density estimation. It includes: 
\itemize{
\item \code{kernel} character vector defining the kernel function used to 
estimate the density;
\item \code{bwtype} character vector defining if a fixed or an adaptive kernel
 estimator has been used;
\item \code{par} list of components defining the parameters used in density 
estimation;
\item \code{estimate} vector of density estimates evaluated at the data points.
}
See \code{\link{kepdf}} for further details.} 
  \item{nc }{An object of class \code{list} defining details about the 
  identification of the connected regions. It includes:
  \itemize{
   \item \code{nc} number of connected sets for each density level set.
   \item \code{p}  vector of level sets, giving the proportions of data with 
   estimated density below a threshold.
   \item \code{id} group label of each point at different sections of the 
   density estimate. Negative values of \code{id} mean that the estimated 
   density is below the considered threshold.
   \item \code{pq} for each \code{p} gives the corresponding quantile \code{q}
   of the values of the density.
  }}
 \item{graph}{An object of class \code{list} defining details about the graph 
 built to find the connected sets of high density regions. Its length
 depends on the value of its first element:
  \itemize{
   \item \code{type} either "unidimensional", "delaunay" or "pairs", 
   defines the procedure used to set edges among the observations. 
   In the last case only, the list includes also the following elements:
   \item \code{comp.area} a list containing the vector \code{area} and the 
   matrix \code{pairs.ord}. The element \eqn{i} of vector \code{area} is the
   measure of the maximum valley in the density function linking the 
   observations having row position as given in column \eqn{i} of 
   \code{pairs.ord}.  		 
   \item \code{lambda} tolerance threshold.    
  }
  }
  \item{cluster.cores}{Vector with the same length as \code{NROW(x)}, defining
  the cluster cores membership. \code{NA} values correspond to low density, 
  unlabeled data, to be classified in the second phase of the procedure by the
  intarnal call of \code{\link{pdfClassification}}.
  }
  \item{tree }{Cluster tree with leaves corresponding to the connected 
  components associated to different sections of the density estimate.
  The object is of class \code{\link[stats:dendrogram]{dendrogram}}.
  }             
  \item{noc }{number of clusters.
  }
  \item{stages }{List with elements corresponding to the data allocation to 
  groups at the different stages of the classification procedure.
   \code{NA} values correspond to unlabeled data.
   }
  \item{clusters }{set to \code{NULL} if \code{n.stages} = 0, that is, if data
  belonging to the cluster cores only have been allocated. 
  Otherwise it reports the final label groups.
  }	
 }


\section{Methods}{
\describe{

\item{\code{signature(x="data.frame", graphtype="ANY", Q="ANY", lambda="ANY",
                      grid.pairs="ANY")}}{

This method applies the \code{pdfCluster} procedure to objects of class 
\code{data.frame}.}
\item{\code{signature(x="matrix", graphtype="ANY", Q="ANY", lambda="ANY",
                       grid.pairs="ANY")}}{

This method applies the \code{pdfCluster} procedure to objects of class 
\code{matrix}.}
\item{\code{signature(x="numeric", graphtype="ANY", Q="ANY", lambda="ANY",
                       grid.pairs="ANY")}}{

This method applies the \code{pdfCluster} procedure to objects of class 
\code{numeric}.}

\item{\code{signature(x="pdfCluster", graphtype="missing", Q="missing", 
                      lambda="ANY", grid.pairs="missing")}}{

This method applies to objects of \code{\link{pdfCluster-class}} when the graph
 has been built according to the "pairs" procedure. It allows to save time and
 computations if the user wants to compare results of cluster analysis for 
 different values of the \code{lambda} parameter. See examples below.}
}}


\references{
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation.
 \emph{Statistics and Computing}. 17, 71-80.

Menardi, G., Azzalini, A. (2012) An advancement in clustering via nonparametric
 density estimation. \emph{Submitted}.
}


\section{Warning }{
Using \code{graphtype= 'delaunay'} when the dimensionality \eqn{d} of data is 
greater than 6 is highly time-consuming unless the number of rows \eqn{n} 
is fairly small. Computational complexity to run the Delaunay triangulation 
is \eqn{O(\frac{n^{\lfloor \frac{d}{2}\rfloor}}{\lfloor \frac{d}{2}\rfloor!})}.
Use \code{graphtype= "pairs"}, instead, whose complexity grows quadratically 
with number of observations.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfCluster-class}}, \code{\link{pdfClassification}}.  
}

\examples{
##########
#example 1
###########
# The code below essentially reproduces the example in Section 4.3 of
# Azzalini and Torelli (2007) up to a minor difference, presumably due
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

#clustering
cl <- pdfCluster(x, h = h.norm(x), hmult=0.75)
summary(cl)
plot(cl)

#comparing groups with original macro-area membership
groups <- cl@clusters
table(oliveoil$macro.area, groups)

#cluster cores
table(cl@cluster.cores)

##########
#example 2
###########

# loading data
data(wine)
x <-wine[ ,-1]
gr <- wine[ ,1]

# when data are high-dimensional, an adaptive kernel estimator is preferable 
# building the Delaunay graph entails a too high computational effort
# use option "pairs" to build the graph 
# it is the default option for dimension >6
# clustering
cl <- pdfCluster(x, graphtype="pairs", bwtype="adaptive")
summary(cl)
plot(cl)

#comparison with original groups
table(cl@clusters,gr)

# a better classification is obtained with larger value of lambda
# not necessary to run the whole procedure again
# a pdfCluster method applies on pdfCluster-class objects!

cl1 <- pdfCluster(cl, lambda=0.27)
table(gr, cl1@clusters)
}

\keyword{methods}
\keyword{cluster}
