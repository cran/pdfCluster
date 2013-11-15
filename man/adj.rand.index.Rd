\name{adj.rand.index}
\alias{adj.rand.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjusted Rand index}
\description{
Computes the adjusted Rand index to compare two altenative partitions of the same set. 
}
\usage{
adj.rand.index(cl1, cl2)
}

\arguments{
  \item{cl1}{
the vector containining the class labels of the first partition.}
	\item{cl2}{
the vector containining the class labels of the second partition.  }

}
\details{
The adjusted Rand index is a correction of the Rand index that measures the similarity between two classifications of
the same objects by the proportions of agreements between the two partitions. The correction is obtained by subtracting from the Rand index its expected value. 
}
\value{
A numeric vector of length 1.    
}
\references{
L. Hubert and P. Arabie (1985) Comparing partitions, \emph{Journal of Classification}, 2, 193-218. } 

\seealso{
\code{\link[base:table]{table}}}

\examples{
# load data
data(wine)
#actual groups
gr <- wine[, 1]

# select a subset of variables
x <- wine[, c(2, 5, 8)]

#clustering
cl <- pdfCluster(x)

#comparison with original groups
table(groups(cl), gr)
adj.rand.index(groups(cl), gr)
}
\keyword{cluster}
