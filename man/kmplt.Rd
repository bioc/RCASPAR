\name{ kmplt }
\alias{ kmplt }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Kaplan Meier curve }
\description{
This function plots the survival curve for the provided data set as a Kaplan Meier plot. It can only be used for visualization and returns no numeric values.
}
\usage{
kmplt(data, title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing \emph{at least} the two columns ``censored'' and ``True_STs''. Where ``censored'' contains the censorship status of the subject as either ``0/F'' for 
\sQuote{uncensored subjects} or ``$1$/T'' for \sQuote{censored subjects}.
This information is essential to be able to plot the KM curve.
}
  \item{title}{
A string of characters denoting the title of the plot produced.
}
}
\details{
Note that this function was intended only for visualization and does not return any numerical values \emph{as such}.
}
\value{
A plot of the survival function ``S(t)'' against the survival times (unit-less) as calculated for the provided data set.
}
\references{
Bland JM, Altman DG. Survival probabilities (the Kaplan-Meier method). \emph{BMJ} 1998;317;1572\\
\url{http://www.bmj.com/statsbk/12.dtl}
}
\author{
Douaa Mugahid
}

\examples{
censored <- c(1, 0, 1, 1, 1, 0, 1, 0, 0, 0)
True_STs <- c(1, 4, 5, 4, 6, 3, 2, 1, 3, 4)
dat <- as.data.frame(cbind(True_STs, censored))
kmplt(dat, "test")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kaplan Meier }
\keyword{ Survival function }
