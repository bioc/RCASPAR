\name{weights_xvBLH}
\alias{weights_xvBLH}
\title{
A special version of STpredictor.BLH used within k-xv to predict the survival times of the kth validation group in the cross validation step.
}
\description{
This function is an \dQuote{incomplete} version of STpredictor.BLH used within the cross validation function \code{STpredictor_xvBLH} to predicted the survival times of the
subset of patients in the kth partitioning. It is not meant for use outside that function.
}
\usage{
weights_xvBLH(geDataS, survDataS, geDataT, survDataT, q = 1, s = 1, a = 2, b = 2, groups = 3, par, method = "BFGS", noprior = 1, extras = list())
}

\arguments{
  \item{geDataS}{
The co-variate data of the kth validation set passed on by \code{STpredictor.xv.BLH}. It is a matrix with the co-variates in the columns and the
subjects in the rows. Each cell corresponds to that row\emph{th} subject's column\emph{th} co-variate's value.
}
  \item{survDataS}{
The survival data of the kth validation set passed on by \code{STpredictor_xvBLH}. It takes on the form of a data frame with at least have the
following columns \dQuote{True_STs} and \dQuote{censored}, corresponding to the observed survival times and the censoring status of the subjects 
consecutively. Censored patients are assigned a \dQuote{1} while patients who experience an event are assigned \dQuote{1}.
}
  \item{geDataT}{
The co-variate data of the kth training set passed on by \code{STpredictor_xvBLH}.
}
  \item{survDataT}{
The survival data of the kth training set passed on by \code{STpredictor_xvBLH}.
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
  \item{par}{
A single vector with the initial values of the baseline hazards followed by the weights(regression coefficients) for the co-variates.
}
  \item{method}{
The preferred optimization method. It can be one of the following:
\code{"Nelder-Mead":} for the Nelder-Mead simplex algorithm.
\code{"L-BFGS-B"} for the L-BFGS-B quasi-Newtonian method.
\code{"BFGS"} for the BFGS quasi-Newtonian method.
\code{"CG"} for the Conjugate Gradient decent method.
\code{"SANN":} for the simulated annealing algorithm.
}
  \item{noprior}{
An integer indicating the number of iterations to be done without assuming a prior on the regression coefficients.
}
  \item{extras}{
The extra arguments to passed to the optimization function optim. For further details on them, see the documentation for the \code{optim} function.
}
}
\value{
\item{prediction }{A data frame with the columns True_STs (the observed survival times), Predicted_STs (the predicted survival times), censored(the censoring status of the 
patient,absolute_error(the sign-less difference between the predicted and observed survival times), PatientOrderValidation (The patient's number)}
\item{est.geneweight }{The estimated regression coefficients from the kth training set (geDataT,survDataT)}
\item{est.baselineH}{The estimated baseline hazards from the kth training set (geDataT, survDataT)}
}
\author{
Douaa Mugahid
}
\note{
This function is not meant to be used outside its wrapper.
}
\seealso{
\code{\link{STpredictor_BLH}}
}
\examples{
data(Bergamaschi)
data(survData)
weights_xvBLH(geDataS=Bergamaschi[21:31, 1:2], survDataS=survData[21:31, 9:10],geDataT=Bergamaschi[1:20, 1:2], 
survDataT=survData[1:20, 9:10], q = 1, s = 1, a = 2, b = 2, groups = 3, par = c(0.1, 0.1, 0.1,rep(0,2)), 
method = "CG", noprior = 1, extras = list(reltol=1))
}

\keyword{survival time prediction}
