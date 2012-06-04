\name{plot-methods}
\docType{methods}
\alias{plot-methods}
%\alias{plot,ANY,ANY-method}
%\alias{plot,kepdf,ANY-method}
%\alias{plot,pdfCluster,missing-method}
%\alias{plot,dbs,missing-method}
\title{Methods for Function plot}
\description{
Methods for functions \code{plot} aimed at graphically displaying the S4 classes
included in the \code{\link{pdfCluster-package}}.  
}

\section{Methods}{
\describe{

\item{\code{signature(x = "kepdf", y = "ANY")}}{

S4 method for plotting objects of \code{\link{kepdf-class}}. See 
\code{\link{plot,kepdf-method}} for further details.}

\item{\code{signature(x = "dbs", y = "missing")}}{

S4 method for plotting objects of \code{\link{kepdf-class}}. See 
\code{\link{plot,dbs-method}} for further details.}

\item{\code{signature(x = "pdfCluster", y = "missing")}}{

S4 method for plotting objects of \code{\link{pdfCluster-class}}. See
 \code{\link{plot,pdfCluster-method}} for further details.}

}}

\seealso{
\code{\link{plot,dbs-method}}, \code{\link{plot,kepdf-method}}, 
\code{\link{plot,pdfCluster-method}} 
}
\keyword{methods}
%%\keyword{ ~~ other possible keyword(s) ~~ }
