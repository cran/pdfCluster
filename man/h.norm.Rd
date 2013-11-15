\name{h.norm}
\alias{h.norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normal optimal choice of smoothing parameter in density estimation
}
\description{
This function computes the smoothing parameter to be used in kernel density estimation, as asymptotically optimal
when the underlying distribution is Normal. Unidimensional as well as multidimensional data
can be handled. When multidimensional data are supplied, a vector of smoothing parameters is computed having one element for each 
component. 
}
\usage{
h.norm(x)
}

\arguments{
  \item{x}{
vector, matrix or data-frame of data.
}
}
\details{
The smoothing parameter of component \eqn{j} of a \eqn{n\times d} data matrix is estimated as follows:
\deqn{\sigma_j{\left(\frac{4}{(d+2)n }\right)}^{\frac{1}{d+4}}}

where \eqn{\sigma_j} is the estimated standard deviation of component \eqn{j}.
See Section 2.4.2 of the reference below.
}
\value{
Returns a numeric vector with the same length as the number of columns of \code{x} or with length one if \code{x} is a vector. When \code{x} is a matrix,
 a vector of smoothing parameters is provided having one element for each component.    
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied smoothing techniques for data analysis: the kernel approach with S-Plus illustrations}. Oxford University Press, Oxford.
} 

\seealso{
\code{\link[sm:hnorm]{hnorm}}}
\examples{
set.seed(123)
x <- rnorm(30)
sm.par <- h.norm(x)
pdf <- kepdf(x, bwtype= "fixed", h = sm.par)
plot(pdf, eval.points=seq(-4,4,by=.2))
}
\keyword{nonparametric}
\keyword{smooth}
