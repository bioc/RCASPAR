#### Function that runs a k-times cross validation for the estimation of weights and prediction of patient survival times
#### Requires the following arguments:
    #### 1) geData = the gene expression data for all patients (rows) corresponding to each gene (cols)
    #### 2) survData = the survival dat for all patients including censorship status and recorded censorship/survival times
    #### 3) k = the number of subsets into which the data is split into for validation
    #### 4) baselineH = the probablitiy that patients would die at any time when all the gene expression values are zero.(some basal value)
    #### 5) cut.off = the value at which patients are seperated into long and short survival groups
    #### 6) file = the path for the directory in which the generated log files are to be stored
    #### 7) q, s = the parameters for the prior on the wieghts
    #### 8) par = the initial values for the parameters, must be specified by the user as a vector of values of length equal to the ncols(geData)
    #### 9) method = the optimization algorithm used (can be any of BFGS, L-BFGS-B, nlm, SANN, Nelder-Mead)
    #### 10) hessian = should a hessian be returned? T/F
    #### 11) fraction = the fraction of initial paramaters equal to zero at which an optimization iteration without a prior is run before the  whole optimization with a prior
    #### 12) extras = additional parameters to be passed to the optimization algorithm
STpredictor_xvBLH <- function (geData, survData, k=10, cut.off, file=paste(getwd(),"STpredictor.xv.BLH_results",sep="/"), q=1, s=1, a=2, b=2, groups= 3, geneweights=NULL, BLHs=NULL, method="BFGS", noprior=1, extras=list()) {
    no.patients <- nrow(geData)
    no.genes <- ncol(geData)
     #### If the estimates for the BLHs passed are not determined by the user, then assign the max. of the pdf of gamma as the initial values
    if(length(BLHs) == 0)
        BLHs <- rep((b * (a - 1)),groups) #;print(BLHs)
    if(length(geneweights) == 0)
        geneweights <- rep(0, no.genes) #; print(geneweights[1:5])
    par <- c(BLHs, geneweights)
    #### Find if the number of patients is divisible by the number of data sets we need or not
    round <- no.patients %% k
    sep <- no.patients %/% k
    ### The boundaries of the groups into which the patients are divided if the number of patients is divisible by k
    if (round == 0) {
        groups_upper <- seq(1, no.patients, sep)
        groups_lower <- seq(sep, no.patients, sep)
    }
    #### The boundaries of the group into which the patients are divided if the number of patients is not divisible by k.
    #### Done in a way such that the one of each of the extra patients is distributed among the first "extras" groups 
    else {
        xtrs <- no.patients - sep * k
        frst <- (sep + 1) * xtrs
        groups_upper <- c(seq(1, frst, sep + 1), seq(frst + 1, no.patients, sep))
        upp <- c(seq(1, frst, sep+1), seq(frst+sep, no.patients, sep))
        groups_lower <- c(seq(sep + 1, frst, sep + 1), seq(frst + sep, no.patients, sep))
    }
    addn <- 0 #### The values added to the the patients' number to keep them updated according to their original order
    final_result <- NULL #### A table with the final predicted_STs and true_STs among other things
    final_long <- NULL #### A table with the patients predicted to be long survivors
    final_short <- NULL #### A table with the patients predicted to be short survivors
    weights <- NULL #### A matrix with the weights predicted from each optimization on the subset used for training
    baselineHs <- NULL #### A matrix with the baselineHs predicted from each optimization on the subset used for training
    mean.weights <-c() #### The average of the weights predicted for each subset used for training
    mean.baselineHs <- NULL #### The average of the baselineHs predicted for each subset used for training
    names_weights <- paste("weights", c(1:k), sep="_") #### The set number from which the weights were estimated
    names_baselineHs <- paste("weights", c(1:k), sep="_")
    #### Now for each new set of validation and training data run an optimization and determine the predicted_STs
    for (group in 1:k) { 
    cat("Progress for group",group,"\n",sep=" ")
        validation_ge <- geData[groups_upper[group]:groups_lower[group], ] #### The gene expression data for the subset used for validation
        validation_surv <- survData[groups_upper[group]:groups_lower[group], ] #### The survival data for the subset used for validation
        train_ge <- geData[- c(groups_upper[group]:groups_lower[group]), ] #### The gene expression data for the subset used for training
        train_surv <- survData[- c(groups_upper[group]:groups_lower[group]), ] #### The survival data for the subset used for training
        result <- weights_xvBLH(geDataS=validation_ge, survDataS=validation_surv, geDataT=train_ge, survDataT=train_surv, q=q, s=s, a=a, b=b, groups=groups, par=par, method=method, noprior=noprior, extras=extras)
        est.weights <- result$est.geneweights #### Extract the weights determined in this run
        est.baselineH <- result$est.baselineH
        result <- result$prediction #### reassign result as the predicion table returned by ST.predictor.xv
        #### Update the values of the patients' number to match that in which they were ordered in geData and survData
        result$PatientOrderValidation <- addn + result$PatientOrderValidation[1:nrow(result)]  
        #### Append the results from this subset to the final result
        final_result <- rbind(final_result, result)
        #### Append the weights estimated from this subset as a row to the final matrix of weights
        weights <- rbind(weights, est.weights)
        baselineHs <- rbind(est.baselineH)
        addn <- addn + nrow(result) #### Update the value of addition for the next subset
    }
    #### Calculate the mean for each of the weights estimated from the different training data subsets
    for(index in 1:ncol(weights)) {
        meani <- mean(weights[ ,index])
        mean.weights <- c(mean.weights, meani)
    }
    #### Calculate the mean for each of the baseline Hazards estimated from the different trainnig subsets
    for (indexj in 1:groups) {
        meanj <- mean(baselineHs[ ,indexj])
        mean.baselineHs <- c(mean.baselineHs, meanj)
   }
    #### Assign the long and short survivors
    long_survivors <- final_result[result$Predicted_STs >= cut.off, ]
    long_survivors <- cbind(long_survivors, rep("L", nrow (long_survivors)))
    colnames (long_survivors)[ncol(long_survivors)] <- "group"
    short_survivors <- final_result[result$Predicted_STs < cut.off, ]
    short_survivors <- cbind(short_survivors, rep("S",nrow (short_survivors)))
    colnames (short_survivors)[ncol(short_survivors)] <- "group"
    #### Write the log file for this xv run
    cat(date(), "\n", "* Function called:", "\n", file=file, sep="")
    sink(file=file, append=TRUE)
    print(args)
    sink()
    cat("\n", "** Gene Weights:", "\n", file=file, append=TRUE)
    sink(file=file, append=TRUE)
    for(i in 1:k) {
        cat(paste("   Group_",i,":","\n",sep=""))
        print(weights[i,])
    }
    cat("\n", "   Average weights", "\n", sep="")
    print(mean.weights)
    sink()
    cat("\n", "*** Baseline Hazards:", "\n", file=file, append=TRUE)
    sink(file=file, append=TRUE)
    for(i in 1:k) {
        cat(paste("   Group_",i,":","\n",sep=""))
        print(weights[i,])
    }
    cat("\n", "   Average baseline Hazards", "\n", sep="")
    print(mean.weights)
    sink()
    cat("\n", "**** Predicted survival times:", "\n", file=file, append=TRUE)
    sink(file=file, append=TRUE)
    print(result)
    sink()
    cat("\n", "***** Long survivors:", "\n", file=file, append=TRUE)
    sink(file=file, append=TRUE)
    print(long_survivors)
    sink()
    cat("\n", "****** Short survivors:", "\n", file=file, append=TRUE)
    sink(file=file, append=TRUE)
    print(short_survivors)
    sink()
    #### Remove the global variables created during the run
    #rm(a.prior, c.prior, e.prior, g.prior, no.patients, no.genes, True_STs, survDataT, geDataT, limits, diff.limits, baselineH_j, pos=1)
    return(list(predicted_STs = final_result, short_survivors = short_survivors, long_survivors = long_survivors, weights=mean.weights, baselineHs=mean.baselineHs))
              }