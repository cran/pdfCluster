\name{plot,pdfCluster-method}
\docType{methods}
\alias{plot.pdfCluster}
\alias{plot,pdfCluster-method}
\alias{plot,pdfCluster,missing-method}

\title{
Plot objects of class pdfCluster
}
\description{
Functions and methods for plotting objects of \code{\link{pdfCluster-class}}.
}
\usage{
\S4method{plot}{pdfCluster}(x, y, which = 1:4, stage = Inf, pch = NULL, col = NULL, ...)
}

\arguments{
  \item{x}{
An object of \code{\link{pdfCluster-class}};
}
  \item{y}{
Not used; for compatibility with generic plot;
}
  \item{which}{
To be used to select the type of plot: \itemize{
		\item when \code{which = 1} plots the mode function, corresponding to the number of modes for different proportions of data points with density above a threshold. 
		\item when \code{which = 2} plots the cluster tree associated to different sections of
			the density estimate.
		\item when \code{which = 3} displays the scatterplot of data or of all the possible pairs of coordinates reporting the label group. 
		\item when \code{which = 4} the \code{\link{plot.dbs}} is displayed. 
		}
Multiple choices are possible.  
  }
  \item{stage}{
Plots the data points at the indicated \code{stage} of the classification procedure. Unallocated data are indicated by 0s. This argument applies if \code{which=3} only. 
}
  \item{pch}{
Either an integer specifying a symbol or a single character to be used in plotting points. If a vector of the same length as the number of groups is given, 
different symbols or characters are used for different groups. The default value denotes points as their 
group label. This argument applies if \code{which=3} only.
}
  \item{col}{
Colors to be used in plotting points. If a vector of the same length as the number of groups is given, 
different colors or characters are used for different groups. The default value use colors in palette corresponding to the the 
group labels of the data. This argument applies if \code{which=3} only.}
  \item{\dots}{
Further arguments to be passed to \code{\link{plot-methods}}.}
}
%\details{
%If necessary, more details than the description above
%}
%\value{
%Value here
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}

\section{Methods}{
\describe{
\item{\code{signature(x = "pdfCluster", y = "missing")}}{

S4 method for plotting objects of \code{\link{pdfCluster-class}}}
}
}
\references{
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation. \emph{Statistics and Computing}. vol. 17, pp. 71-80.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfCluster-class}}, \code{\link[graphics:plot]{plot}}, \code{\link{plot-methods}}.
}
\examples{
data(wine)
gr <- wine[,1]

# select a subset of variables
x <- wine[, c(2,5,8)]

#clustering
cl <- pdfCluster(x)
plot(cl, which=3, stage=2 )

table(cl@clusters, gr)
#set "B" for Barolo, "G" for Grignolino, "A" for Barbera
plot(cl,  pch=c("B", "G", "A"), col=c(3,4,5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
