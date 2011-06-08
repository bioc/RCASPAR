\name{survData}
\alias{survData}
\docType{data}
\title{
Survial data of 82 patients
}
\description{
A dataframe ''survData:" A dataframe containing the survival information of the 82 patients the relevant columns of which are
"True_STs", "censored" indicating the recorded survival times of the patients, and their censorship status simultaneously. In the
"censored" column "0" indicates a non-censored patient, and "1" indicates a censored patient\\
}
\usage{data(survData)}
\format{
  A data frame with 82 observations on the following 10 variables.
  \describe{
    \item{\code{Array.ID}}{The array IDs of the patients}
    \item{\code{Sample.ID.Bauke.}}{The array IDs of the patients as annotated in one source}
    \item{\code{Experiment.Jon.}}{The array IDs of the patients as annotated in a third source (in essence the same as the previous two)}
    \item{\code{ExptID}}{The Experiment ID}
    \item{\code{Order}}{The number/order of the patient}
    \item{\code{Overall.survival..mons.undivided}}{the overall survival time of the patients in months}
    \item{\code{Relapse.free.survival..mons.}}{the relapse free survival time of the patients in months}
    \item{\code{Status.0.A..1.AWD..2.DOD..3.DOC}}{the status of the tumour (as abbreviated in the publication)}
    \item{\code{censored}}{the status of censorship with ''0" corresponding to a ''non-censored" patient and ''1" to a ''censored patient" }
    \item{\code{True_STs}}{The survival/censorship times of the patient in years}
  }
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
data(survData)
colnames(survData)
}
\keyword{datasets}
