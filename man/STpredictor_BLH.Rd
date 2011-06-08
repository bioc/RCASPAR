\name{STpredictor_BLH}
\alias{STpredictor_BLH}
\title{Predicts the survival times of the validation set based on the regression coefficients and baseline hazards determined according to the Piecewise baseline hazard Cox 
regression model.}
\description{
This function uses the training set to estimated the best regression coefficients, and baseline hazards describing the data according to the piecewise baseline hazard Cox 
regression model. It then takes them and uses them to predict the survival times of the validation set, which are determined as the mean value of the p.d.f. of the survival 
time, as a continuous random variable, given the co-variate values of that subject.
}
\usage{
STpredictor_BLH(geDataS, survDataS, cut.off, file = paste(getwd(), "STpredictor_results", sep = "/"), geDataT, survDataT, groups = NULL, a = 2, b = 2, q = 1, s = 1, BLHs = 
NULL, geneweights = NULL, method = "BFGS", noprior = 1, extras = list())
}

\arguments{
  \item{geDataS}{
The co-variate data of the validation set passed on by the user. It is a matrix with the co-variates in the columns and the subjects in the rows. Each cell 
corresponds to that row\emph{th} subject's column\emph{th} co-variate's value.
}
  \item{survDataS}{
The survival data of the validation set passed on by the user. It takes on the form of a data frame with at least have the following columns \dQuote{True\_STs} and 
\dQuote{censored}, corresponding to the observed survival times and the censoring status of the subjects consecutively. Censored patients are assigned a \dQuote{1} while 
patients who experience an event are assigned \dQuote{1}.
}
  \item{cut.off}{
The value of the separator around which the patients are grouped according to their predicted survival times.
}
  \item{file}{
The path of the file to which the log file of this session is saved.
}
  \item{geDataT}{
The co-variate data of the kth training set passed on by the user.
}
  \item{survDataT}{
The survival data of the kth training set passed on by the user.
}
  \item{groups}{
The number of partitions along the time axis for which a different baseline hazard is to be assigned. This number should be the same as the number of initial values passed for 
the baseline hazards in the beginning of the \dQuote{weights\_baselineH} argument.
}
  \item{a}{
The shape parameter for the gamma distribution used as a prior on the baseline hazards.
}
  \item{b}{
The scale parameter for the gamma distribution used as a prior on the baseline hazards.
}
  \item{q}{
One of the two parameters on the prior distribution used on the weights (regression coefficients) in the model.
}
  \item{s}{
The second of the two parameters on the prior distribution used on the weights (regression coefficients) in the model.
}
  \item{BLHs}{
A vector with the initial values for the baseline hazards. Should be of length \emph{groups}. The default is NULL, in which case a vector of length \code{groups} with values
corresponding to the maximum of the gamma distributions with the given parameters is created.
}
  \item{geneweights}{
A vector with the initial values of the weights(regression coefficients) for the co-variates. The default is NULL, in which case a vector of zeros the same length as 
\code{ncol(geData)} is created as the initial starting value.
}
  \item{method}{
The preferred optimization method. It can be one of the following:\\
\code{"Nelder-Mead":} for the Nelder-Mead simplex algorithm.\\
\code{"L-BFGS-B":} for the L-BFGS-B quasi-Newtonian method.\\
\code{"BFGS":} for the BFGS quasi-Newtonian method.\\
\code{"CG":} for the Conjugate Gradient decent method\\
\code{"SANN":} for the simulated annealing algorithm.\\
}
  \item{noprior}{
An integer indicating the number of iterations to be done without assuming a prior on the regression coefficients.
}
  \item{extras}{
The extra arguments to passed to the optimization function optim. For further details on them, see the documentation for the \code{optim} function.
}
}
\value{
\item{log_optimization}{The result of the optimization performed on the training set as is described in the documentation for the \code{optim} function}
\item{short_survivors}{A data frame of results for the patients living less than the cut off value; with the columns True\_STs (the observed survival times), 
Predicted\_STs (the predicted survival times), censored(the censoring status of the patient, absolute\_error(the sign-less difference between the predicted and observed 
survival times), PatientOrderValidation (The patient's number)}
\item{long_survivors}{A data frame with the results for the patients living at least as long as the cut off value; with columns True\_STs (the observed survival times), Predicted\_STs (the predicted survival times), censored(the censoring status of the 
patient,absolute\_error(the sign-less difference between the predicted and observed survival times), PatientOrderValidation (The patient's number)}
}
\references{
The basic model is based on the Cox regression model as first introduced by Sir David Cox in: Cox,D.(1972).Regression models & life tables. \emph{Journal of the Royal 
Society of Statistics}, 34(2), 187-220.
The extension of the Cox model to its stepwise form was adapted from: Ibrahim, J.G, Chen, M.-H. & Sinha, D. (2005). \emph{Bayesian Survival Analysis (second ed.)}. NY: Springer.
as well as Kaderali, Lars.(2006) A Hierarchical Bayesian Approach to Regression and its Application to Predicting Survival Times in Cancer Patients. Aachen: Shaker 
The prior on the regression coefficients was adopted from: Mazur, J., Ritter,D.,Reinelt, G. & Kaderali, L. (2009). Reconstructing Non-Linear dynamic Models of Gene Regulation 
using Stochastic Sampling. \emph{BMC Bioinformatics}, 10(448).
}
\author{
Douaa Mugahid
}

\seealso{
\code{\link{STpredictor_xvBLH}}
}
\examples{
data(Bergamaschi)
data(survData)
result <- STpredictor_BLH(geDataS=Bergamaschi[1:27, 1:2], survDataS=survData[1:27, 9:10], geDataT=Bergamaschi[28:82, 1:2], survDataT=survData[28:82, 9:10], q = 1,
s = 1, a = 1.558, b = 0.179, cut.off=3, groups = 3, method = "CG", noprior = 1, extras = list(reltol=1))
}

\keyword{survival time prediction}
\keyword{piecewise baseline hazard Cox PH regression model }
