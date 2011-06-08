#### requires weights.estimator and weights_fun_extended_ST
STpredictor_BLH <- function (geDataS, survDataS, cut.off, file=paste(getwd(),"STpredictor_results",sep="/"),geDataT, survDataT, groups=NULL, a=2, b=2, q=1, s=1, BLHs=NULL, geneweights=NULL, method="BFGS", noprior=1, extras=list()) {
    args <- match.call (expand.dots=TRUE)
    #### Assign some variables:
    no.patients <- nrow (geDataT)
    no.genes <- ncol (geDataT)
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
    if(length(BLHs) == 0)
        BLHs <- rep((b * (a - 1)),groups)
    if(length(geneweights) == 0)
        geneweights <- rep(0, length(geDataS))
    par <- c(BLHs, geneweights)
    est.weights <- weights_BLH (geDataT=geDataT, survDataT=survDataT, q=q, s=s, a=a, b=b, groups=groups, par=par, method=method, noprior=noprior, extras=extras) 
    est.baselineH <- est.weights$par[1:groups] 
    prods <- est.baselineH * diff.limits
# function for the first step( t1-t2). t1=0
    func1 <- function (l,t,m) {
            -exp(-m*l*t) * (m*l*t + 1) / (m*l)
             }
         
# function of the higher limit for the steps between (t2 - tn), tn != Inf
    func2 <- function (sum, t, tt, l, m) {
             -exp(-m*(sum -t*l)) * -exp(-m*l*tt) * (m*l*tt + 1) / (m*l)
             }
         
# function of the lower limit for the steps between (t2- tn), tn != Inf
    func3 <- function (l, ll, t, m, sum) {
             l*exp(-m*(sum - t*ll)) * -exp(-m*ll*t) * (m*ll*t + 1) / (m*ll^2)
             }
         
# function for the last step between (tn - Inf)
    func4 <- function (sum, l, m, t) {
             -exp(-m*sum) * l*m*t^2 / 2
             }
    predicted_STs <-rep(NA,nrow(geDataS))
    for (j in 1:nrow(geDataS)) { #### j is the patient's number
    m <- exp(sum(est.weights$par[-c(1:groups)]*geDataS[j, ]))
    expected <- func1(l=est.baselineH[1],t=limits[2],m=m)
    for (n in 2:(groups-1)) {
        expected <- expected + func2(sum=sum(prods[1:(n-1)]), t=limits[n], l=est.baselineH[n], tt=limits[n+1], m=m) - func3(sum=sum(prods[1:(n-1)]), t=limits[n], l=est.baselineH[n-1], ll= est.baselineH[n], m=m)
        }
    expected <- expected + func4(m=m, l=est.baselineH[groups], t=limits[groups+1], sum=sum(prods[1:groups]))
    predicted_STs[j] <- expected
    }
    result <- as.data.frame (cbind (c(1:nrow (survDataS)), survDataS$True_STs, predicted_STs, abs (survDataS$True_STs - predicted_STs), survDataS$censored))
    colnames(result) <- c("PatientOrderValidation","True_STs","Predicted_STs","Absolute_Error","censored")
    long_survivors <- result[result$Predicted_STs >= cut.off, ]
    long_survivors <- cbind(long_survivors, rep ("L", nrow (long_survivors)))
    colnames (long_survivors)[ncol (long_survivors)] <- "group"
    short_survivors <- result[result$Predicted_STs < cut.off, ]
    short_survivors <- cbind (short_survivors, rep ("S",nrow (short_survivors)))
    colnames (short_survivors)[ncol (short_survivors)] <- "group"
    cat (date(), "\n", "*Function called:", "\n", file=file, sep="")
    sink (file=file, append=TRUE)
    print (args)
    sink ()
    cat ("\n", "**Baseline Hazard(s):", "\n", file=file, append=TRUE)
    sink(file=file, append=TRUE)
    print(est.weights[1:groups])
    sink()
    cat ("\n", "***Weights:", "\n", file=file, append=TRUE)
    sink (file=file, append=TRUE)
    print (est.weights[-c(1:groups)])
    sink ()
    cat ("\n", "****Predicted survival times:", "\n", file=file, append=TRUE)
    sink (file=file, append=TRUE)
    print (result)
    sink ()
    cat ("\n", "*****Long survivors:", "\n", file=file, append=TRUE)
    sink (file=file, append=TRUE)
    print (long_survivors)
    sink ()
    cat ("\n", "******Short survivors:", "\n", file=file, append=TRUE)
    sink (file=file, append=TRUE)
    print (short_survivors)
    sink ()
    #rm(limits, True_STs, survDataT, geDataT, no.patients, no.genes, pos=1)
    return (list (log_optimization=est.weights, predicted_STs=result, short_survivors=short_survivors, long_survivors=long_survivors))         
             }
 