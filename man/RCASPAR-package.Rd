\name{RCASPAR-package}
\alias{RCASPAR-package}
\alias{RCASPAR}
\docType{package}
\title{
A package for survival time prediction based on a piecewise baseline hazard Cox regression model.
}
\description{
The package is the R-version of the C-based software \bold{CASPAR} (Kaderali,2006). It is meant to help predict survival times in the presence of high-dimensional explanatory 
co-variates. The model is a piecewise baseline hazard Cox regression model with an Lq-norm based prior that selects for the most important regression coefficients, and in turn 
the most relevant co-variates for survival analysis. It was primarily tried on gene expression and aCGH data, but can be used on any other type of high-dimensional data and in 
disciplines other than biology and medicine.  
}
\details{
\tabular{ll}{
Package: \tab RCASPAR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-08-23\cr
License: GPL(>=3)
LazyLoad: \tab yes\cr
}
}
\author{
Douaa Mugahid

Maintainer: Douaa Mugahid <mugahid@stud.uni-heidelberg.de>,  Lars Kaderali <lars.kaderali@bioquant.uni-heidelberg.de>
}
\references{
The basic model is based on the Cox regression model as first introduced by Sir David Cox in: Cox,D.(1972).Regression models & life tables. \emph{Journal of the Royal Society
of Statistics}, 34(2), 187-220.
The extension of the Cox model to its stepwise form was adapted from: Ibrahim, J.G, Chen, M.-H. & Sinha, D. (2005). \emph{Bayesian Survival Analysis (second ed.)}.
NY: Springer.
as well as Kaderali, Lars.(2006) A Heirarchial Bayesian Approach to Regression and its Application to Predicting Survival Times in Cancer Patients. Aachen: Shaker 
The prior on the regression coefficients was adopted from: Mazur, J., Ritter,D.,Reinelt, G. & Kaderali, L. (2009). Reconstructing Non-Linear dynamic Models of Gene 
Regulation using Stochastic Sampling. \emph{BMC Bioinformatics}, 10(448).
}

\keyword{Piecewise baseline hazard Cox regression model}
\keyword{survival analysis}

\examples{
## Eg.(1): A simple example performed with a training and validation set:
data(Bergamaschi)
data(survData)
  ## Generate prediction:
result <- STpredictor_BLH(geDataS=Bergamaschi[1:27, 1:2], survDataS=survData[1:27, 9:10], geDataT=Bergamaschi[28:82, 1:2], survDataT=survData[28:82, 9:10], q = 1, s = 1, a = 1.558, b = 0.179
, cut.off=15, groups = 3, method = "CG", noprior = 1, extras = list(reltol=1))
 ## Plot a KM plot with both long and short survivors:
kmplt_svrl(long=result$long_survivors, short=result$short_survivors,title="KM plot of long and short survivors")
  ## Determine the area under the curve of AUROC curves vs. time to see the performance of the predictor given the chosen parameters and the current partitioning into training
  ## and validation sets:
survivAURC(Stime=result$predicted_STs$True_STs,status=result$predicted_STs$censored, marker=result$predicted_STs$Predicted_STs, time.max=20)
 ## Perform a log-rank test to see if the difference between the long and short survivors is significant:
logrnk(dataL=result$long_survivors, dataS=result$short_survivors)

## Eg.(2): A simple example performed with cross validation:
data(Bergamaschi)
data(survData)
  ## Generate prediction:
STpredictor_xvBLH(geData=Bergamaschi[1:40,1:2], survData=survData[1:40,9:10], k = 10, cut.off = 10, q = 1, s = 1, a = 1.558, b = 0.179, groups = 3, method = "BFGS", noprior = 1, extras = list(reltol=1))
 ## Plot a KM plot with both long and short survivors:
kmplt_svrl(long=result$long_survivors, short=result$short_survivors,title="KM plot of long and short survivors")
  ## Determine the area under the curve of AUROC curves vs. time to see the performance of the predictor given the chosen parameters and the current partitioning into training
  ## and validation sets:
survivAURC(Stime=result$predicted_STs$True_STs,status=result$predicted_STs$censored, marker=result$predicted_STs$Predicted_STs, time.max=20)
 ## Perform a log-rank test to see if the difference between the long and short survivors is significant:
logrnk(dataL=result$long_survivors, dataS=result$short_survivors)
}