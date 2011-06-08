\name{pltgamma}
\alias{pltgamma}
\title{ Plotting the gamma distribution of shape parameter }
\description{ This function allows the user to visualize the shape of the gamma distribution used as a prior on the baseline hazards in the 
functions ending with "_BLH"
}
\usage{
pltgamma(a, b)
}
\arguments{
  \item{a}{Is the shape parameter for the gamma distribution.
}
  \item{b}{Is the scale parameter for the gamma distribution
}
}
\value{
A plot of the gamma distribution with the given parameters.
}

\author{
Douaa Mugahid
}

\seealso{
\code{\link{pltprior}}
}
\examples{
pltgamma(a=2,b=2)
}
\keyword{gamma distribution}
