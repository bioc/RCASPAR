\name{weight_estimator_BLH}
\alias{weight_estimator_BLH}
\title{
Returns the value of the objective function used for optimizing for the regression parameters and baseline hazards in the model.
}
\description{
Given the arguments, it can evaluate the value of the objective function used by the optimization algorithms for determining the optimal regression parameters and baseline 
hazard values.
}
\usage{
weight_estimator_BLH(survDataT, geDataT, weights_baselineH, q, s, a, b, groups)
}

\arguments{
  \item{survDataT}{
The survival data of the patient set passed on by the user. It takes on the form of a data frame with at least have the following columns \dQuote{True_STs} and 
\dQuote{censored}, corresponding to the observed survival times and the censoring status of the subjects consecutively. Censored patients are assigned a \dQuote{1} while 
patients who experience an event are assigned \dQuote{1}.
}
  \item{geDataT}{
The co-variate data (gene expression or aCGH, etc...) of the patient set passed on by the user. It is a matrix with the co-variates in the columns and the subjects in the rows. Each cell 
corresponds to that row\emph{th} subject's column\emph{th} co-variate's value.
}
  \item{weights_baselineH}{
A single vector with the initial values of the baseline hazards followed by the weights(regression coefficients) for the co-variates.
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
The number of partitions along the time axis for which a different baseline hazard is to be assigned. This number should be the same as the number of initial values passed for 
the baseline hazards in the beginning of the \dQuote{weights_baselineH} argument.
}
}

\value{
A single numerical value corresponding to the value of the objective function with the given regression coefficients and baseline hazard values.
}
\references{
The basic model is based on the Cox regression model as first introduced by Sir David Cox in: Cox,D.(1972).Regression models & life tables. \emph{Journal of the Royal Society 
of Statistics}, 34(2), 187-220.
The extension of the Cox model to its stepwise form was adapted from: Ibrahim, J.G, Chen, M.-H. & Sinha, D. (2005). \emph{Bayesian Survival Analysis (second ed.)}.
NY: Springer.
as well as Kaderali, Lars.(2006) A Hierarchial Bayesian Approach to Regression and its Application to Predicting Survival Times in Cancer Patients. Aachen: Shaker 
The prior on the regression coefficients was adopted from: Mazur, J., Ritter,D.,Reinelt, G. & Kaderali, L. (2009). Reconstructing Non-Linear dynamic Models of Gene Regulation
using Stochastic Sampling. \emph{BMC Bioinformatics}, 10(448).
}
\author{
Douaa Mugahid
}
\note{
This function is in itself not useful to the user, but is used within the function \code{weights.BLH}
}

\seealso{
\code{\link{weight_estimator_BLH_noprior}}, \code{\link{deriv_weight_estimator_BLH_noprior}}
}
\examples{
\donttest{
data(Bergamaschi)
data(survData)
weight_estimator_BLH(survDataT=survData[1:10, 9:10], geDataT=Bergamaschi[1:10 , 1:2],weights_baselineH=c(0.1,0.2,0.3,rep(0,ncol(Bergamaschi))), q=1, s=1, a=1.5, b=0.3, groups=3)
}
}

\keyword{ Cox regression model}
\keyword{ piecewise baseline hazard cox regression model}