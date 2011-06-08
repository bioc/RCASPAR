\name{STpredictor_xvBLH}
\alias{STpredictor_xvBLH}
\title{
This function performs a cross validation on the full data set to help predict the survival times of the patients using the piecewise baseline hazard PH Cox model.
}
\description{
Using the full data provided by the user, this function splits the data set k times, into a smaller validation set, and a much bigger training set. The regression coefficients 
of the model are estimated from the training set and used to predict the survival times of the validation set. The patients can then be split into patients two groups around a 
cut off value also specified by the user.
}
\usage{
STpredictor_xvBLH(geData, survData, k = 10, cut.off, file = paste(getwd(), "STpredictor.xv.BLH_results", sep = "/"), q = 1, s = 1, a = 2, b = 2, groups = 3, geneweights = NULL
, BLHs = NULL, method = "BFGS", noprior = 1, extras = list())
}
\arguments{
  \item{geData}{
A matrix with the co-variate data of the full set of subjects. It is constructed with the co-variate in the columns and the subjects in the rows.Each cell corresponds to that 
row\emph{th} subject's column \emph{th} co-variate's value.
}
  \item{survData}{
The survival data of the entire set of subjects. It takes on the form of a data frame with at least have the following columns \dQuote{True_STs} and \dQuote{censored}, 
corresponding to the
observed survival times and the censoring status of the subjects consecutively. Censored patients are assigned a \dQuote{1} while patients who experience an event are assigned 
\dQuote{1}.
}
  \item{k}{
The number of times the cross-validation is.
}
  \item{cut.off}{
The value of the separator around which the patients are grouped according to their predicted survival times.
}
  \item{file}{
The path of the file to which the log file of this session is saved.
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
the baseline hazards in the beginning of the \dQuote{weights_baselineH} argument
}
  \item{geneweights}{
A vector with the initial values of the weights(regression coefficients) for the co-variates. The default is NULL, in which case a vector of zeros the same length as 
\code{ncol(geData)} is created as the initial starting value.
}
  \item{BLHs}{
A vector with the initial values for the baseline hazards. Should be of length \emph{groups}. The default is NULL, in which case a vector of length \code{groups} with values
corresponding to the maximum of the gamma distributions with the given parameters is created.
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
The extra arguments to passed to the optimization function optim. For further details on them, see the documentation for the \kbd{optim} function.
}
}
\value{
\item{predicted_STs}{A data frame of the results for all patients, with the columns True_STs (the observed survival times), Predicted_STs (the predicted survival times), 
censored(the censoring status of the 
patient,absolute_error(the signless difference between the predicted and oberved survival times), PatientOrderValidation (The patient's number)}
\item{short_survivors}{A data frame of results for the patients living less than the cut off value; with the columns True_STs (the observed survival times), 
Predicted_STs (the predicted survival times), censored(the censoring status of the patient,absolute_error(the signless difference between the predicted and oberved survival 
times), PatientOrderValidation (The patient's number)}
\item{long_survivors}{A data frame with the results for the patients living at least as long as the cut off value; with columns True_STs (the observed survival times), 
Predicted_STs (the predicted survival times), censored(the censoring status of the patient,absolute_error(the sign-less difference between the predicted and observed survival 
times), PatientOrderValidation (The patient's number)}
\item{weights}{A vector with the mean value of the regression coefficients obtained from the k training sets}
\item{baselineHs}{A vector with the mean value of the baseline hazards returned from the k training sets}
}
\references{
The basic model is based on the Cox regression model as first introduced by Sir David Cox in: Cox,D.(1972).Regression models & life tables. \emph{Journal of the Royal 
Society of Statistics}, 34(2), 187-220.
The extension of the Cox model to its stepwise form was adapted from: Ibrahim, J.G, Chen, M.-H. & Sinha, D. (2005). \emph{Bayesian Survival Analysis (second ed.)}. NY: Springer.//
as well as Kaderali, Lars.(2006) A Hierarchial Bayesian Approach to Regression and its Application to Predicting Survival Times in Cancer Patients. Aachen: Shaker 
The prior on the regression coefficients was adopted from: Mazur, J., Ritter,D.,Reinelt, G. & Kaderali, L. (2009). Reconstructing Non-Linear dynamic Models of Gene 
Regulation using Stochastic Sampling. \emph{BMC Bioinformatics}, 10(448).
}
\author{
Douaa Mugahid
}

\seealso{
\code{\link{STpredictor_BLH}}
}
\examples{
data(Bergamaschi)
data(survData)
STpredictor_xvBLH(geData=Bergamaschi[1:20, 1:2], survData=survData[1:20, 9:10], k = 10, cut.off=3, file = paste(getwd(), "STpredictor.xv.BLH_results", sep = "/"), q = 1, s = 1, a = 2, b = 2, 
groups = 3, geneweights = NULL, BLHs = NULL, method = "CG", noprior = 1, extras = list(reltol=1))
}
\keyword{ cross validation}
\keyword{ piecewise baseline hazard Cox PH model}
