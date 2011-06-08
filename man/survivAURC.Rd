\name{survivAURC}
\alias{survivAURC}
\title{
A function that calculates the area under a curve constructed from plotting the area under a ROC curve at the corresponding time point at which it was generated.
}
\description{
This function plots the individual area under the ROC curves of different time points against the time at which they were evaluated, and calculates the area under that curve.
}
\usage{
survivAURC(Stime, status, marker, entry = NULL, cut.values = NULL, time.max = 20, by = 1)
}
\arguments{
  \item{Stime}{
The observed survival times of the patients.
}
  \item{status}{
The censoring status of the patient. 1 for a censored patient, and 0 for a patient who has an event.
}
  \item{marker}{
The predicted survival time of the patients.
}
  \item{entry}{
The time of entry of the patients, set to NULL by default.
}
  \item{cut.values}{
The cut off values for which the ROC curves are to be constructed.
}
  \item{time.max}{
The maximum time point for which the area under a ROC curve is to be plotted.
}
  \item{by}{
The step size between every ROC curve estimated and the next.
}
}
\details{
The calculations for the Area under the ROC curves at each time point is done according to Patrick Heagerty's survivalROC function in the survivalROC package in R. It is a 
value indicating the performance of the predictor when used on the results of both functions STpredictor.BLH and STpredictor.xv.BLH.
}
\value{
\item{ AUC }{The value of the area under the curve generated}
\item{AUeachROC}{A vector with the values of the area under the individual ROC curves}
A plot of the Area under the ROC curves against their corresponding time points.
}
\references{
Heagerty,P., Lumely T. & Pepe M.(2000). Time-dependent ROC curves for censored survival data & a diagnostic marker. Biometrics, 56(2), 337-344.
}
\author{
Douaa Mugahid
}

\seealso{
\code{\link{survivROC}}
}
\examples{
True_STs <- c(1.416667,2.75,2.416667,2.583333,2.166667,2.5,2.5,1.833333,1.25,0.6666667,1,6.583333,6.5,6.666667,2.75,1.666667,1.166667,2.833333,3.583333,6.166667,6.166667,
3.416667,6.083333,1.833333,5.583333,0.75,5.75,5.5,0.5833333,7.666667,5,2.833333,1.333333,5.083333,0.8333333,1.5,4.75,3.416667,4.666667,1.916667,4.666667,7.416667,0.9166667,
1.083333,3.75,3.25,3,2.416667,2.75,2.5,2.666667,4.5,4.416667,1.5,0.8333333,3.166667,3.833333,3.833333,0.4166667,3.333333,2.75,3.083333,0.3333333,0.25,0.6666667,1.833333,
2.333333,3.416667,3.416667,3,0.6666667,0.75,2.166667,1,1.416667,1.333333,1.166667,1.166667,0.4166667,1.25,1.166667,1.083333)
Predicted_STs <- c(6.030591,6.014457,3.545584,5.414229,6.41576,9.393992,5.542331,6.890859,8.090213,4.98545,2.77357,6.275699,9.163978,7.511511,9.531218,7.63715,10.08977,
11.12364,3.982502,5.441881,12.61404,12.21851,17.05850,12.78141,16.22795,21.48544,6.281354,13.83925,8.859929,6.104142,8.255909,2.335526,6.564962,2.335761,9.33772,12.62540,
10.97276,15.63089,8.01967,5.817267,5.59897,4.340784,32.40319,33.74123,27.45024,26.31024,26.88833,24.34707,32.06541,38.90473,17.37102,15.11059,8.772035,14.24816,7.852889,
7.79996,5.601459,2.802408,35.77047,24.34717,30.65796,25.93927,20.64544,22.04807,19.15037,23.83430,1.876557,3.937208,6.526354,5.886377,9.301074,12.4657,14.49783,15.41502,
2.860931,2.541947,4.543111,4.525553,4.148272,3.986912,6.246755,6.89523)
censored <- c(0,1,1,1,1,0,1,1,0,1,0,1,1,1,1,0,0,0,0,1,1,1,1,0,1,0,1,1,0,1,1,0,0,1,0,0,1,0,1,0,1,1,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,
1)
survivAURC(Stime=True_STs,status=censored, marker=Predicted_STs,time.max=5)
}

\keyword{ ROC curves }
\keyword{Area under ROC curves}
