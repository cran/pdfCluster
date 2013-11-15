\name{hprop2f}
\alias{hprop2f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample smoothing parameters in adaptive density estimation
}
\description{
This function computes the sample smoothing parameters to be used in adaptive kernel density estimation, according to Silverman (1986).   
}
\usage{
hprop2f(x, h = h.norm(x), alpha = 1/2, kernel = "gaussian")
}

\arguments{
	\item{x}{
Vector or matrix of data.}
	\item{h}{
Vector of smoothing parameters to be used to get a pilot estimate of the density function. It has length equal to \code{NCOL(x)}.
}
	\item{alpha}{
Sensitivity parameter satysfying \eqn{0 \leq \alpha \leq 1}, giving the power to which raise the pilot density. Default value is 1/2.
See details.	
}
	\item{kernel}{
Kernel to be used to compute the pilot density estimate. It should be one of 
"gaussian" or "t7". See \code{\link{kepdf}} for further details.	
}

}
\details{
A vector of smoothing parameters \eqn{h_{i}} is chosen for each sample point \eqn{x_i}, as follows:
\deqn{h_i = h \left(\frac{\hat{f}_h(x_i)}{g}\right)^{- \alpha }}

where \eqn{\hat{f}_h} is a pilot kernel density estimate of the density function \eqn{f}, with vector of bandwidths \code{h},
 and \eqn{g} is the geometric mean of \eqn{\hat{f}_h(x_i)},
 \eqn{i=1, ..., n}.
See Section 5.3.1 of the reference below.
}
\value{
Returns a matrix with the same dimensions of \code{x} where row \eqn{i} provides
the vector of smoothing parameters for sample point \eqn{x_i}.    
}

\references{
Silverman, B. (1986). \emph{Density estimation for statistics and data analysis}. Chapman and Hall, London.
} 

\seealso{
\code{h.norm}
}
\examples{
set.seed(123)
x <- rnorm(10)

sm.par <- hprop2f(x)
pdf <- kepdf(x, bwtype= "adaptive")

pdf@par$hx
sm.par

plot(pdf,eval.points=seq(-4,4,by=.2))
}
\keyword{nonparametric}
\keyword{smooth}
