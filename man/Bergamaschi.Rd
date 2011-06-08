\name{Bergamaschi}
\alias{Bergamaschi}
\docType{data}
\title{
Gene expression data of 82 patients with 10 genes as covariates
}
\description{
A dataframe:\\
''Bergamaschi:" A dataframe of 10 covariates (in the columns) for 82 patients (in the rows)\\
}

\usage{data(Bergamaschi)}
\format{
  The format is:
 num [1:82, 1:10] 0.654 0.701 0.126 0.899 0.267 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:82] "1" "2" "3" "4" ...
  ..$ : chr [1:10] "IMAGE:753234" "IMAGE:50794" "IMAGE:302190" "IMAGE:51408" ...
}
\details{
A subset of the data set used in Bergamaschi, A., & al., e. (2006). Distinct Patterns of DNA Copy NumberAlteration are 
associated with different clinicopathological features and gene expression subtypes of breast cancer. Genes, Chromosomes
and cancer , 45, 1033-1040.
}
\source{
Bergamaschi, A., & al., e. (2006). Distinct Patterns of DNA Copy NumberAlteration are associated with different 
clinicopathological features and gene expression subtypes of breast cancer. Genes, Chromosomes and cancer , 45, 1033-1040.
}
\references{
Bergamaschi, A., & al., e. (2006). Distinct Patterns of DNA Copy NumberAlteration are associated with different 
clinicopathological features and gene expression subtypes of breast cancer. Genes, Chromosomes and cancer , 45, 1033-1040.
}
\examples{
data(Bergamaschi)
colnames(Bergamaschi)
}
\keyword{datasets}
