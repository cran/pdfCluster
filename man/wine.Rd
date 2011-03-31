\name{wine}
\alias{wine}
\title{Wine data}
\description{
These data are the results of a chemical analysis of wines grown in the same region in Italy but derived from three
different cultivars. The analysis determined the quantities of 13 constituents found in each of the three types of wines: Barolo, Grignolino, Barbera. 
The data set is used to evaluate the \code{pdfCluster} ability of understanding the type of wine, given the chemical measurement. 
}

\usage{data(wine)}

\format{
This data frame contains 178 rows, each corresponding to a different cultivar of wine produced in Piedmont (Italy), and 14 columns. 
The first column is the type of wine, a factor variable with the following levels: Barolo, Grignolino, Barbera.
The attributes of the three types of wines are the following: Alcohol, Malic acid, Ash, Alcalinity, Magnesium, Phenols, Flavanoids,
Nonflavanoids, Proanthocyanins, Color intensity, Hue, Dilution, Proline. All attributes are continuous.}
            
\details{
The original data set had around 28 variables, here a subset of variables only has been considered.
}

\source{
Forina, M. et al, PARVUS - \emph{An Extendible Package for Data
Exploration, Classification and Correlation}. Institute of Pharmaceutical
and Food Analysis and Technologies, Via Brigata Salerno, 
16147 Genoa, Italy.
}
 
\references{
Forina M., Armanino, C., Castino, M. Ubigli, M. (1986) 
Multivariate data analysis as a discriminating method of the origin of wines. \emph{Vitis} 25, 189-201. }
  
\keyword{datasets}