\name{h.norm}
\alias{h.norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normal optimal choice of smoothing parameter in density estimation
}
\description{
This function computes the smoothing parameter to be used in kernel density estimation, as asymptotically optimal
when the underlying distribution is Normal. Unidimensional as well as multidimensional data
can be handled. When multidimensional data are provided the diagonal elements of the smoothing matrix are estimated. 
}
\usage{
h.norm(x)
}

\arguments{
  \item{x}{
vector or matrix of data
}
}
\details{
The smoothing parameter of dimension \eqn{j} of a \eqn{n\times d} data matrix is estimated as follows:
\deqn{\sigma_j{\left(\frac{4}{(d+2)n }\right)}^{\frac{1}{d+4}}}

where \eqn{\sigma_j} is the estimated standard deviation of dimension \eqn{i}.
See Section 2.4.2 of the reference below.
}
\value{
Returns a vector with the same length as the number of columns of \code{x} or with length one if \code{x} is a vector. When \code{x} is a matrix, the diagonal elements of the 
asymptotically optimal smoothing matrix for Normal data are provided.    
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for Data Analysis: Kernel Approach with S-Plus Illustrations}. Oxford University Press, Oxford.
} 

\seealso{
\code{\link[sm:hnorm]{hnorm}}}
\examples{
set.seed(123)
x <- rnorm(30)
sm.par <- h.norm(x)
pdf <- kepdf(x,h=sm.par)
plot(pdf,eval.points=seq(-4,4,by=.2))
}
\keyword{nonparametric}
\keyword{smooth}