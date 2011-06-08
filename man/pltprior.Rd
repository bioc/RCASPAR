\name{pltprior}
\alias{pltprior}
\title{
A function to visualize the shape of the prior on the weights with the chosen q and s parameters.
}
\description{
This function helps visualize the effect of the parameters \emph{q} and \emph{s} on the prior distribution used on the weights.
}
\usage{
pltprior(q, s)
}
\arguments{
  \item{q}{
One of the two parameters used to determine the prior on the weights.
}
  \item{s}{
The second of two parameters used to determine the prior on the weights.
}
}
\details{
The prior assumed on the weights in the objective function takes on the following formulation:
\deqn{[L_q(\beta; q, s)= \frac{\frac{{q^{(q-1)}} {q}} {2s\Gamma(\frac{{1}{q}})}} exp(\frac{{-1}{qs^q}} |\emph{weights}|^q)}
}
\value{
A 3D plot of the value of the prior taking into consideration only two weights}
\references{
Mazur J., et al. Reconstructing nonlinear dynamic models of gene regulation using stochastic sampling. \emph{BMC Bioinformatics} 2009;10;448
}
\author{
Douaa Mugahid
}

\seealso{
\code{\link{pltgamma}}
}
\examples{
pltprior(q=1,s=1)}

\keyword{Lq-Norm prior}