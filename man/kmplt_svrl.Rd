\name{kmplt_svrl}
\alias{kmplt_svrl}
\title{
A function that plots the KM curves of $2-3$ patient sets in one graph.
}
\description{
This function can plot the KM curves estimated for $2-3$ patients simultaneously for sake of easier comparison.
}
\usage{
kmplt_svrl(all, long, short, title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all}{
Data for the first set of patients; usually the complete set of patients, but could be any other.It is a data frame containing \emph{at least} the two columns ``censored'' and 
``True_STs''. Where ``censored'' contains the censorship status of the subject as either ``0/F'' for \bold{uncensored subjects} or ``1/T'' for \bold{censored subjects}.
This information is essential to be able to plot the KM curve.
}
  \item{long}{
Data for the second set of patients; in our case the group of patients who survived at least up to the cut off value passed to the predictor. It has essentially the same 
structure as \dQuote{all}
}
  \item{short}{
Data for the third and last set of patients; in our case the group of patients who survived less than the cut off value passed on to the predictor. It essentially has the same 
structure as the two other arguments.
}
  \item{title}{
The main title for the plot.
}
}
\details{
This function essentially is the same as \code{kmplt} but does the same for up to $3$ plots simultaneously.}
\value{
A plot with all $2-3$ KM curves.
}
\references{
\url{http://www.bmj.com/statsbk/12.dtl}
}
\author{
Douaa Mugahid
}


\seealso{
\code{\link{kmplt}}
}
\examples{
censored <- c(1, 0, 1, 1, 1, 0, 1, 0, 0, 0)
True_STs <- c(1, 4, 5, 4, 6, 3, 2, 1, 3, 4)
dat1 <- as.data.frame(cbind(True_STs, censored))
censored <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 1)
True_STs <- c(7, 7, 8, 5, 9, 11, 8, 11, 10, 6)
dat2 <- as.data.frame(cbind(True_STs, censored))
censored <- c(1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1)
True_STs <- c(1, 4, 5, 4, 6, 3, 2, 1, 3, 4, 7, 7, 8, 5, 9, 11, 8, 11, 10, 6)
dat3 <- as.data.frame(cbind(True_STs, censored))
kmplt_svrl(all=dat3, long=dat2, short=dat1, title="KM of predictions")           
}
\keyword{ Kaplan-Meier estimator}
