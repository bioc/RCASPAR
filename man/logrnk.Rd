\name{logrnk}
\alias{logrnk}
\title{ Performs Log Rank test on the long and short patient sets
}
\description{
This function performs a Chi-square test on the long and short subject sets to determine if their is a significant difference between the survival
times in both sets. It returns the p-value.
}
\usage{
logrnk(dataL, dataS)
}
\arguments{
  \item{dataL}{
The set of subjects predicted to fall into the long-survivor set.
A data frame containing at least the following columns: ``PatientOrderValidation'' (the number/order of the subject); ``group'' (the group into which
the patient falls L (for long) or S(for short)); ``censored'' (the censorship status of the patient $1$ for uncensored and $0$ for censored).
}
  \item{dataS}{
Same as ``dataL'' but for the set of short survivors.
}
}
\details{
Note that the typical arguments to be passed are the results of the ``STpredict'' functions ``long\_survivors'' and ``long\_survivors''
}
\value{
The estimated p-value is returned
}
\references{
Bland JM, Altman DG. Survival probabilities (the Kaplan-Meier method). \emph{BMJ} 2004;328;1073
\url{http://www.bmj.com/statsbk/12.dtl}
}
\author{
Douaa AS Mugahid
}

\seealso{
\code{\link{survivAURC}}
}
\examples{
PatientOrderValidation_L <- c(1, 2, 3, 5, 7)
PatientOrderValidation_S <- c(4, 6, 8)
group_L <- rep("L", 5)
group_S <- rep("S", 3)
censored_L <- c(0, 0, 1, 1, 0)
censored_S <- c(0, 0, 1)
True_STs_L <- c(5, 6, 6 ,7, 8)
True_STs_S <- c(2, 3, 2)
short <- as.data.frame(cbind(PatientOrderValidation_S, group_S, censored_S, True_STs_S))
long <- as.data.frame(cbind(PatientOrderValidation_L, group_L, censored_L, True_STs_L))
names(short) <- c("PatientOrderValidation", "group", "censored", "True_STs")
names(long) <- c("PatientOrderValidation", "group", "censored", "True_STs")
logrnk(dataL=long, dataS=short)
}

\keyword{ log rank test }
