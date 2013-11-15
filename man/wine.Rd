\encoding{UTF-8}
\name{wine}
\alias{wine}
\title{Wine data}
\description{
These data are the results of a chemical analysis of wines grown in the same region in Italy but derived from three
different cultivars. The analysis determined the quantities of 13 constituents found in each of the three types of wine: Barolo, Grignolino, Barbera. 
The data set is used to evaluate the \code{pdfCluster} ability of understanding the type of wine, given the chemical measurement. 
}

\usage{data(wine)}

\format{
This data frame contains 178 rows, each corresponding to a different cultivar of wine produced in Piedmont (Italy), and 14 columns. 
The first column is the type of wine, a factor variable with the following levels: Barolo, Grignolino, Barbera.
The variables measured on the three types of wines are the following: Alcohol, Malic acid, Ash, Alcalinity, Magnesium, Phenols, Flavanoids,
Nonflavanoids, Proanthocyanins, Color intensity, Hue, OD280.OD315Dilution, Proline. All variables but the label class are continuous.}
            
\details{
The original data set comprises 27 variables. Here a subset of 14 variables only has been included.
}

\source{
Forina, M., Lanteri, S. Armanino, C., Casolino, C., Casale, M., Oliveri, P. (2008). V-PARVUS. \emph{An Extendible Pachage of programs for 
esplorative data analysis, classification and regression analysis}. Dip. Chimica e Tecnologie Farmaceutiche ed Alimentari, \enc{Università}{Universita} di Genova. Freely available at
\url{http://www.parvus.unige.it}. 
}
 
%\references{
%Forina M., Armanino, C., Castino, M. Ubigli, M. (1986) 
%Multivariate data analysis as a discriminating method of the origin of wines. \emph{Vitis} 25, 189-201. }
  
\keyword{datasets}
