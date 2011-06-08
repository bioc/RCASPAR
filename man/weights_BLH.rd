\name{weights_BLH}
\alias{weights_BLH}
\title{
Optimization for the regression coefficients and baseline hazards that maximize the partial likelihood in our PW Cox PH regression model.
}
\description{
This function is a wrapper around the optimization function \code{optim} to allow the optimization for the regression coefficients and baseline
hazards appropriate for the data set at hand. It is where the functions \code{weight_estimator_BLH, weight_estimator_BLH_noprior,
deriv_weight_estimator_BLH, deriv_weight_estimator_BLH_noprior} are required.
}
\usage{
weights_BLH(geDataT, survDataT, q, s, a, b, groups, par, method = c("Nelder-Mead", "L-BFGS-B", "CG", "BFGS", "SANN"), noprior = 1, extras = list(),
dist = NULL)
}
\arguments{
  \item{geDataT}{
A matrix with the co-variate in the columns and the subjects in the rows.Each cell corresponds to that row\emph{th} subject's column\emph{th} 
co-variate's value.
}
  \item{survDataT}{
A data frame with the survival data of the set of subjects at hand. It should at least have the following columns \dQuote{True_STs} and
\dQuote{censored}, corresponding to the observed survival times and the censoring status of the subjects consecutively. Censored patients are assigned a \dQuote{1} while patients who experience an event are assigned
\dQuote{1}.
}
  \item{q}{
One of the two parameters on the prior distribution used on the weights (regression coefficients) in the model.
}
  \item{s}{
The second of the two parameters on the prior distribution used on the weights (regression coefficients) in the model.
}
  \item{a}{
The shape parameter for the gamma distribution used as a prior on the baseline hazards.
}
  \item{b}{
The scale parameter for the gamma distribution used as a prior on the baseline hazards.
}
  \item{groups}{
The number of partitions along the time axis for which a different baseline hazard is to be assigned. This number should be the same as the number
of initial values passed for the baseline hazards in the beginning of the \dQuote{weights_baselineH} argument.
}
  \item{par}{
A single vector with the initial values of the baseline hazards followed by the weights(regression coefficients) for the co-variates.
}
  \item{method}{
The preferred optimization method. It can be one of the following:
\code{"Nelder-Mead":} for the Nelder-Mead simplex algorithm.
\code{"L-BFGS-B":} for the L-BFGS-B quasi-Newtonian method.
\code{"BFGS":} for the BFGS quasi-Newtonian method.
\code{"CG":} for the Conjugate Gradient decent method
\code{"SANN":} for the simulated annealing algorithm.
}
  \item{noprior}{
An integer indicating the number of iterations to be done without assuming a prior on the regression coefficients.
}
  \item{extras}{
The extra arguments to passed to the optimization function optim. For further details on them, see the documentation for the \code{optim} function.
}
  \item{dist}{
The distribution function to be passed to the optimization algorithm in case of using SANN to generate a new candidate point.
}
}
\value{
The same value as the \code{optim} function. See it's documentation for details.
}
\references{
\url{http://sekhon.berkeley.edu/stats/html/optim.html}
}
\author{
Douaa Mugahid
}
\note{
Note that this function is just a wrapper around the \code{optim} function to serve our purpose, and it's main purpose is to be called within the main functions of this package 
\code{STpredictor_BLH} and \code{weights_xvBLH}
}
\examples{
data(Bergamaschi)
data(survData)
weights_BLH(geDataT=Bergamaschi[1:10,1:2], survDataT=survData[1:10, 9:10], q=1, s=1, a=1.56, b=0.17, groups=3, par=c(0.1,0.2,0.3,rep(0,ncol(Bergamaschi))), method = "CG", noprior = 1, extras =
list(reltol=1), dist = NULL)
}
\keyword{optimization}
\keyword{likelihood function}
\keyword{Piecewise baseline hazard Cox regression model}
