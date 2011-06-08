#### requires weights.estimator and weights_fun_extended_ST
#### determines the predicted_STs for the given data in a k-xv run
weights_xvBLH <- function (geDataS, survDataS, geDataT, survDataT, q=1, s=1, a=2, b=2, groups=3, par, method="BFGS", noprior=1, extras=list()) {
     #### Some variable definition:
     #### Determine the groups into which the survival times are marked
        #### Reorder the patients and their data according to the ascending True_STs
        ooo <- order(survDataT$True_STs)
        True_STs <- survDataT$True_STs[ooo]
        survDataT <- survDataT[ooo, ]
        geDataT <- geDataT[ooo, ]
    #### Determine the groups into which the survival times are marked
         max_ST <- max(True_STs) + 1
         limits <- seq(0, max_ST, length.out=(groups + 1))
         diff.limits <- limits[-1] - limits[-(groups + 1)]
     #### If the estimates for the BLHs passed are not determined by the user, then assign the max. of the pdf of gamma as the initial values
    est.weights <- weights_BLH (geDataT=geDataT, survDataT=survDataT, q=q, s=s, a=a, b=b, groups=groups, par=par, method=method, noprior=noprior, extras=extras)$par
    est.baselineH <- est.weights[1:groups]
    prods <- est.baselineH * diff.limits
   #### function for the first step( t1-t2). t1=0
        #### l is the baselineH, t is the 1st time point and m is exp(<weights,gene expression values>)
    func1 <- function (l,t,m) {
            -exp(-m*l*t) * (m*l*t + 1) / (m*l)
             }
         
    #### function of the higher limit for the steps between (t2 - tn), tn != Inf
      ##### sum is the integral up to the lower time point t l is the baselineH, t is the time point and m is exp(<weights,gene expression values>), tt is the upper limit on the intergral(time point) 
    func2 <- function (sum, t, tt, l, m) {
             -exp(-m*(sum -t*l)) * -exp(-m*l*tt) * (m*l*tt + 1) / (m*l)
             }
         
    #### function of the lower limit for the steps between (t2- tn), tn != Inf
    func3 <- function (l, ll, t, m, sum) {
             l*exp(-m*(sum - t*ll)) * -exp(-m*ll*t) * (m*ll*t + 1) / (m*ll^2)
             }
         
    #### function for the last step between (tn - Inf)
    func4 <- function (sum, l, m, t) {
             -exp(-m*sum) * l*m*t^2 / 2
             }
    predicted_STs <-rep(NA,nrow(geDataS))
    for (j in 1:nrow(geDataS)) { 
    #### j is the patient's number
    m <- exp(sum(est.weights[-c(1:groups)]*geDataS[j, ]))
    #### initiate the calculation of the predicted ST to be equal to func1 for m, the first baselineH and the first non-zero limit on the time intervals for the baselineH groups
    expected <- func1(l=est.baselineH[1],t=limits[2],m=m)
    #### for the intervals between the first and last add the values of func2 and func3 with the appropriate baselineH values, integrals of the BLH function and so on.
    for (n in 2:(groups-1)) {
        expected <- expected + func2(sum=sum(prods[1:(n-1)]), t=limits[n], l=est.baselineH[n], tt=limits[n+1], m=m) - func3(sum=sum(prods[1:(n-1)]), t=limits[n], l=est.baselineH[n-1], ll= est.baselineH[n], m=m)
    }
    #### in the end add the the value of func4 with the last BLH, limit and integral as arguments.
    expected <- expected + func4(m=m, l=est.baselineH[groups], t=limits[groups+1], sum=sum(prods[1:groups]))
    predicted_STs[j] <- expected
    }
    #### put the results into a dataframe with the appropriate column names
    result <- as.data.frame (cbind (c(1:nrow (survDataS)), survDataS$True_STs, predicted_STs, abs (survDataS$True_STs - predicted_STs), survDataS$censored))
    colnames(result) <- c("PatientOrderValidation","True_STs","Predicted_STs","Absolute_Error","censored")
    #rm(limits, True_STs, survDataT, geDataT, no.patients, no.genes, pos=1)
    return (list(prediction=result, est.geneweights=est.weights[-c(1:groups)], est.baselineH=est.baselineH))
             }