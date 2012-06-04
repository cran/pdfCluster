\name{kepdf}
\alias{kepdf}

\title{
Kernel estimate of a probability density function. 
}
\description{
Estimates density of 1-d, 2-d, 3-d, and higher dimensional data by kernel methods. 
}
\usage{
kepdf(x, eval.points = x, kernel = "gaussian", 
      bwtype = "fixed", h = h.norm(x), hx = NULL, alpha = 1/2)}

\arguments{
  \item{x}{
A vector, a matrix or data-frame of data whose density should be estimated.  
}
  \item{eval.points}{
A vector, a matrix or a data-frame of data points at which the density estimate should be evaluated.  
}
  \item{kernel}{
Either 'gaussian' or 't7', it defines the kernel function to be used. See details below.
}
  \item{bwtype}{
Either 'fixed' or 'adaptive', corresponding to a kernel estimator with fixed or adaptive bandwidths respectively. 
See details below.
}
  \item{h}{
A vector of length set to \code{NCOL(x)}, defining the smoothing parameters to be used either to estimate
the density in kernel estimation with fixed bandwidth or to estimate the pilot density in kernel estimation 
with adaptive bandwidths. Default value is the result of \code{\link{h.norm}} applied to \code{x}. 
}
  \item{hx}{
A matrix with the same number of rows and columns as \code{x}, where each row defines the vector of smoothing parameters 
specific for each sample point. To be used when \code{bwtype = "adaptive"}. Default value is the result of \code{\link{hprop2f}} 
applied to \code{x}. Set to NULL when \code{bwtype= "fixed"}. 
  }
  \item{alpha}{
Sensitivity parameter to be given to \code{\link{hprop2f}} when \code{bwtype= "adaptive"} and the vectors of smoothing parameters are computed according to the Silverman's (1986) approach.
%when satysfying \eqn{0 \leq \alpha \leq 1}, that van be given   
%giving the power to which raise the pilot 
%density  . See details.
}
}

\details{
The current version of \code{\link{pdfCluster-package}} allows for computing estimates by a kernel product 
estimator of the form: 

\deqn{\hat{f}(y)= \sum_{i=1}^n \frac{1}{n h_{i,j} \cdots h_{i,d}} \prod_{j=1}^d K\left(\frac{y_{j} - x_{i,j}}{h_{i,j}}\right).
}

The kernel function \eqn{K} can either be a Gaussian density (if \code{kernel = "gaussian"}) or a \eqn{t_\nu} density, with \eqn{\nu = 7} degrees of freedom (when \code{kernel = "t7"}).
Although uncommon, the opportunity of selecting a \eqn{t} kernel is motivated by computational efficiency reasons. Hence, its use is suggested when either \code{x} or \code{eval.points} have a huge number of rows.
 
The vectors of bandwidths \eqn{h_{i} = (h_{i,1} \cdots h_{i,d})'} are defined as follows: 
\describe{
\item{Fixed bandwidth}{ When \code{bwtype='fixed'}, \eqn{h_{i} = h} that is, a constant smoothing 
vector is used for all the observations \eqn{x_i}. Default values are set as asymptotically optimal for a multivariate Normal distribution (e.g., Bowman and Azzalini, 1997).
See \code{\link{h.norm}} for further details.}
\item{Adaptive bandwidth}{When \code{bwtype='adaptive'}, a vector of bandwidths \eqn{h_i} is  
specified for each observation \eqn{x_i}. Default values are selected according to Silverman 
(1986, Section 5.3.1). See \code{\link{hprop2f}}.
%, as follows:
%\deqn{h_i = h \left(\frac{\hat{f}(x_i)}{g}\right)^{- \alpha }}
%where \eqn{\hat{f}} is a pilot kernel density estimate of the density function \eqn{f}, with vector of bandwidths \code{h}
%and \eqn{g} is the geometric mean of \eqn{\hat{f}_h(x_i)} \eqn{i=1, ..., n}.
%See Section 5.3.1 of Silverman (1986).
}
}

}

\value{
An S4 object of \code{\link{kepdf-class}} with slots:

  \item{call }{The matched call.}   
  \item{x }{The data input, coerced to be a matrix.} 
  \item{eval.points }{The data points at which the density is evaluated.}
  \item{estimate }{The values of the density estimate at the evaluation points.}
  \item{kernel }{The selected kernel.}
  \item{bwtype }{The type of estimator.} 
  \item{par }{A list of parameters used to estimate the density, with elements:
\itemize{
\item \code{h} the smoothing parameters used to estimate either the density or the pilot density;
\item \code{hx}  the matrix of sample smoothing parameters, when \code{bwtype='adaptive'};
\item \code{alpha}  sensitivity parameter used if \code{bwtype='adaptive'}.
}
}
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for Data Analysis: Kernel Approach with S-Plus Illustrations}. Oxford University Press, Oxford.

Silverman, B. (1986). \emph{Density estimation for statistics and data analysis}. Chapman and Hall, New York.
}

\seealso{
\code{\link{h.norm}}, \code{\link{hprop2f}}, \code{\link{kepdf-class}}.
}
\examples{
## A 1-dimensional example
data(wine)
x <- wine[,3] 
pdf <- kepdf(x, eval.points=seq(0,7,by=.1))
plot(pdf, n.grid= 100, main="wine data")

## A 2-dimensional example
x <- wine[,c(2,8)] 
pdf <- kepdf(x)
plot(pdf, main="wine data", props=c(5,50,90), ylim=c(0,4))
plot(pdf, main="wine data", method="perspective", phi=30, theta=60)

### A 3-dimensional example
x <- wine[,c(2,3,8)] 
pdf <- kepdf(x)
plot(pdf, main="wine data", props=c(10,50,70), gap=0.2)
plot(pdf, main="wine data", method="perspective", gap=0.2, phi=30, theta=10)

### A 6-dimensional example
### adaptive kernel density estimate is preferable in high-dimensions
x <- wine[,c(2,3,5,7,8,10)]
pdf <- kepdf(x, bwtype="adaptive")
plot(pdf, main="wine data", props=c(10,50,70), gap=0.2)
plot(pdf, main="wine data", method="perspective", gap=0.2, phi=30, theta=10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{nonparametric}
