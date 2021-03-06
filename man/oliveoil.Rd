\encoding{UTF-8}
\name{oliveoil}
\alias{oliveoil}
\title{Olive oil data}
\description{
This data set represents eight chemical measurements on different specimen of olive oil produced in various regions in Italy (northern Apulia, 
southern Apulia, Calabria, Sicily, inland Sardinia and coast Sardinia, eastern and western Liguria, Umbria) and further classifiable 
into three macro-areas: Centre-North, South, Sardinia. 
The data set is used to evaluate the \code{pdfCluster} ability of recunstructing the macro-area membership. 
}

\usage{data(oliveoil)}

\format{
This data frame contains 572 rows, each corresponding to a different specimen of olive oil, and 10 columns. 
The first and the second column correspond to the macro-area and the region of origin of the olive oils respectively; here, the term "region" refers to a geographical area and only partially to administrative borders. 
Columns 3-10 represent the following eight chemical measurements on the acid components for the oil specimens:
palmitic, palmitoleic, stearic, oleic, linoleic, linolenic, arachidic, eicosenoic.}

\details{
Since the raw data are of compositional nature, ideally totalling 10000, some preliminary transformations of data are advisable. In particular, Azzalini and Torelli (2007) 
adopt an additive log-ratio transformation (ALR). If \eqn{x_j} denotes the \eqn{j^{th}} chemical measurement \eqn{(j=1,\ldots,8)}, the ALR transformation is \eqn{y_j= \log x_j/x_k, j\neq k},
where \eqn{k} is an arbitrary but fixed variable.  
However, in this data set, the raw data do not always sum up exactly to 10000, because of measurement errors. Moreover, some 0's are present in the data, corresponding 
to measurements below the instrument sensitivity level. Therefore, it is suggested to add 1 to all raw data and normalize them by dividing each entry by the
corresponding row sum \eqn{\sum_j (x_j+1)}.  
}

\source{
Forina, M., Lanteri, S. Armanino, C., Casolino, C., Casale, M., Oliveri, P. (2008). V-PARVUS. \emph{An Extendible Package of programs for explorative data analysis, classification and regression analysis}. Dip. Chimica e Tecnologie Farmaceutiche ed Alimentari, \enc{Università}{Universita} di Genova. 
}
 
\references{
Azzalini A., Torelli N. (2007). Clustering via nonparametric density estimation.
 \emph{Statistics and Computing}, 17, 71-80.
%Forina M., Armanino, C., Lanteri, S., Tiscornia, E. (1983) \emph{Classification of olive oils from their fatty acid %composition}. 
%In: Martens, M., Russwurm, H. Jr. (eds). Food Research and data Analysis, Applied Science Publishers, London, 189-214.
}
  

 



\keyword{datasets}
