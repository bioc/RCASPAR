#### Formulation of the posterior function for the determination of gene expression weights according to the training data
weight_estimator_BLH <- function (survDataT, geDataT, weights_baselineH, q, s, a, b, groups) { 
cat(".")
    #### baselineH = the selected value for baseline hazard value (the posibility of death within the next infinitsimaly small time step in case gene expression data is 0)
    #### weights = a vector of values to be multiplied by the corresponding gene expression values per patient that determines the extent of influence a gene has on the patients predicted survival time if any
    #### geDataT = a dataframe/matrix of values of gene expression values for each gene on the microarray chip, for each patient. colnames = Gene Name/Clone ID, rownames = patient number
    #### survDataT = a dataframe with the observed survival time of the patients and the censorship status, 0(true),1(false). Header=Truse_STs|censored, at least.
    #### q, s = the parameters affecting the distribution of the prior on the weights
    
    #### first level of calculations: Independant of weights or gene expression data; done only once in weights function and assigned as global variable
            #### Determine the groups into which the survival times are marked
        #### Reorder the patients and their data according to the ascending True_STs
             ooo <- order(survDataT$True_STs) 
             True_STs <- survDataT$True_STs[ooo]
             survDataT <- survDataT[ooo, ]
             geDataT <- geDataT[ooo, ]    
             max_ST <- max(True_STs) + 1
             limits <- seq(0, max_ST, length.out=(groups + 1))
             no.patients <- nrow(geDataT) 
             a.prior <- 1/(q*s^q) 
             e.prior <- (1-a) 
             baselineH_j <- rep(NA, no.patients)
                 for(pat in 1:no.patients) {
                     for(lim in 1:groups) {
                         if(True_STs[pat] <= limits[lim + 1] && True_STs[pat] > limits[lim]) {
                            baselineH_j[pat] <- lim
                            break
                         }
                     }
                 }
     
    
    #### second level of calculations: Dependant on weights, but independant of gene expression data of patient
        #### Seperate the vector of parameters into those designated as lambdas and those meant to be the weights
        #### The first "group" parameters are the baseline hazards, and the rest are the weights
        weights <- weights_baselineH[-c(1:groups)]
        baselineH <- exp(weights_baselineH[1:groups])   ###### CHANGED DM-LK
    #### for a given weight vector,caluclate |weight(i)|^q and sum up for all weights in the vector
        b.prior <- sum (abs (weights)^q)
        prior <- a.prior*b.prior
        #### for the given baseline hazards, calculate the parts of -ln(the gamma distributions pdf)
        ee.prior <- sum(e.prior * weights_baselineH[1:groups])
        f.prior <- sum(baselineH/b)
    #### third level of calculations: Dependant on weights and gene expression data of the individual patients
    #### for each patient's gene expression data, calculate the value of e^<weights,geDataT of patient>
         prob.weights.total <- 0 #### the probability that the chosen weights fits the data of all the different patients' geDataT to their survival times. Reset to zero for each new chosen vector of weights
         sums <- (c(0, (limits[-1] - limits[-(groups + 1)]))) * (c(0, baselineH)) #### the product of the limits and the corresponding baseline hazards
         for (j in 1:no.patients) { #### j is the patient's number
         #### for each patient's gene expression data, calculate the value of e^<weights,geDataT of patient>
             #### gedPatient <- geDataT[j,]
             weight.ge <- sum (weights * geDataT[j,]) 
             exp.weight.ge <- exp (weight.ge) 
             #### for each patient calculate the contribution of that patient's GEdata to the objective and add it to the older value
             for (grp in 1:groups) { 
                 if ( True_STs[j] > limits[grp] && True_STs[j] <= limits[(grp + 1)]) {
                     integral_BLH <- sum(c(sums[1:(grp)], ((True_STs[j] - limits[grp])  * baselineH[grp])))
                 }    
             }
             if (survDataT$censored[j] == 1)
                 prob.weights.total <- prob.weights.total + integral_BLH * exp.weight.ge
                 # was prob.weights.total <- prob.weights.total + integral_BLH[j] * exp.weight.ge
             else
                 prob.weights.total <- prob.weights.total - log(baselineH[baselineH_j][j]) - weight.ge + integral_BLH * exp.weight.ge
                 # was prob.weights.total <- prob.weights.total - ln_baselineH_j[j] - weight.ge + integral_BLH[j] * exp.weight.ge
         }
         prob.weights.total <- prob.weights.total + prior + ee.prior + f.prior 
         #print(c(baselineH,prob.weights.total))
         return (prob.weights.total)
                  }

#### A formulation for the determination of the derivative of the posterior on the estimated weights
deriv_weight_estimator_BLH <- function (geDataT, survDataT, weights_baselineH, q, s, a, b, groups) 
{
    cat(".")
    #### baselineH = the selected value for baseline hazard value (the posibility of death within the next infinitsimaly small time step in case gene expression data is 0)
    #### weights = a vector of values to be multiplied by the corresponding gene expression values per patient that determines the extent of influence a gene has on the patients predicted survival time if any
    #### geDataT = a dataframe/matrix of values of gene expression values for each gene on the microarray chip, for each patient. colnames = Gene Name/Clone ID|Data, rownames = patient number
    #### survDataT = a dataframe with the observed survival time of the patients and the censorship status, 1(T) or 0(F). Header = patient|True_STs|censored
    #### q, s = the parameters affecting the distribution of the prior on the weights
     # Assign variables:
     ooo <- order(survDataT$True_STs) 
             True_STs <- survDataT$True_STs[ooo]
             survDataT <- survDataT[ooo, ]
             geDataT <- geDataT[ooo, ]    
             max_ST <- max(True_STs) + 1
             limits <- seq(0, max_ST, length.out=(groups + 1))
     no.patients <- length(weights_baselineH[-c(1:groups)])
    # allocate return vector
    deriv <- rep(0,length(weights_baselineH))
    # compute length of blh interval
    delta <- limits[length(limits)] / groups
    # compute position of weights
    startw <- groups + 1
    endw <- length(weights_baselineH)
    # exp-transform baseline hazard
    blh <- exp(weights_baselineH[1:groups])
    # compute derivative with respect to prior on blh
    deriv[1:groups] <- (1-a) + blh/b
    # derivative with respect to prior on weights
    fac <- s^q
    deriv[startw:endw] <- (abs(weights_baselineH[startw:endw])^(q-1))/fac*sign(weights_baselineH[startw:endw])
    deriv[is.na(deriv[startw:endw])] <- 0
    # add the derivatives of the likelihood
    summa <- rep(0,(groups+1))
    for (i in 1:groups)
      summa[i+1] <- summa[i] + delta * blh[i]
    # loop over all patients
    for (i in 1:no.patients)
      for (j in 1:groups)
      {
        # is patient i censored in interval j?
        # if not, there is nothing to do
        a <- limits[j]
        b <- limits[j+1]
        t <- True_STs[i]
        lambda <- blh[j]
        fx <- sum(weights_baselineH[startw:endw]*geDataT[i,])
        if ((t>a)&&(t<=b)){
          if (survDataT$censored[i]==0)
            deriv[startw:endw] <- deriv[startw:endw] - geDataT[i,]
          deriv[startw:endw] <- deriv[startw:endw] + geDataT[i,] * exp(fx) * (lambda * (t-a) + summa[j])
        }
        # deriv. w.r.t. blh
        if ((t>a)&&(t<=b))
        {
          if (survDataT$censored[i]==0)
            deriv[j] <- deriv[j] -1;
          deriv[j] <- deriv[j] + exp(fx) * (t-a) * lambda
        }
        else if (t>b)
          deriv[j] <- deriv[j] + exp(fx) * delta * lambda
      }
    return(deriv)
}


#### Formulation of the posterior function for the determination of gene expression weights according to the training data with no prior
   weight_estimator_BLH_noprior <- function (geDataT, survDataT, weights_baselineH, a, b, groups) { 
    cat(".")
    #### baselineH = the selected value for baseline hazard value (the posibility of death within the next infinitsimaly small time step in case gene expression data is 0)
    #### weights = a vector of values to be multiplied by the corresponding gene expression values per patient that determines the extent of influence a gene has on the patients predicted survival time if any
    #### geDataT = a dataframe/matrix of values of gene expression values for each gene on the microarray chip, for each patient. colnames = Gene Name/Clone ID, rownames = patient number
    #### survDataT = a dataframe with the observed survival time of the patients and the censorship status, 0(true),1(false). Header=Truse_STs|censored, at least.
    #### q, s = the parameters affecting the distribution of the prior on the weights
    
    #### first level of calculations: Independant of weights or gene expression data; done only once in weights function and assigned as global variable
            #### Determine the groups into which the survival times are marked
        #### Reorder the patients and their data according to the ascending True_STs
             ooo <- order(survDataT$True_STs) 
             True_STs <- survDataT$True_STs[ooo] 
             survDataT <- survDataT[ooo, ]
             geDataT <- geDataT[ooo, ]    
             max_ST <- max(True_STs) + 1  
             limits <- seq(0, max_ST, length.out=(groups + 1))
             no.patients <- nrow(geDataT) 
             e.prior <- (1-a) 
             baselineH_j <- rep(NA, no.patients)
                 for(pat in 1:no.patients) {
                     for(lim in 1:groups) {
                         if(True_STs[pat] <= limits[lim + 1] && True_STs[pat] > limits[lim]) {
                            baselineH_j[pat] <- lim
                            break
                         }
                     }
                 } 
       
    #### second level of calculations: Dependant on weights, but independant of gene expression data of patient
        #### Seperate the vector of parameters into those designated as lambdas and those meant to be the weights
        #### The first "group" parameters are the baseline hazards, and the rest are the weights
        weights <- weights_baselineH[-c(1:groups)]
        baselineH <- exp(weights_baselineH[1:groups])
        #### for the given baseline hazards, calculate the parts of -ln(the gamma distributions pdf)
        ee.prior <- sum(e.prior * weights_baselineH[1:groups])
        f.prior <- sum(baselineH/b)
    #### third level of calculations: Dependant on weights and gene expression data of the individual patients
         prob.weights.total <- 0 #### The value of the objective func. Reset to zero for each new chosen vector or weights
         #### for each patient calculcate the contribution to the objective function and add it to the previous value
         for (j in 1:no.patients) { #### j is the patient's number
         #### for each patient's gene expression data, calculate the value of e^<weights,geDataT of patient>
             weight.ge <- sum (weights * geDataT[j, ]) 
             exp.weight.ge <- exp (weight.ge) 
             sums <- (c(0, (limits[-1] - limits[-(groups + 1)]))) * (c(0, baselineH))
             for (grp in 1:groups) { 
                 if ( True_STs[j] > limits[grp] && True_STs[j] <= limits[(grp + 1)]) {
                     integral_BLH <- sum(c(sums[1:(grp)], ((True_STs[j] - limits[grp])  * baselineH[grp])))
                 }    
             }
             if (survDataT$censored[j] == 1)
                 prob.weights.total <- prob.weights.total + integral_BLH * exp.weight.ge
             else
                 prob.weights.total <- prob.weights.total - log(baselineH[baselineH_j][j]) - weight.ge + integral_BLH * exp.weight.ge
         }
         prob.weights.total <- prob.weights.total + ee.prior + f.prior
         #print(c(baselineH,prob.weights.total))
         return (prob.weights.total)
                  }

#### A formulation for the determination of the derivative of the posterior on the estimated weights with no prior
deriv_weight_estimator_BLH_noprior <- function (survDataT, geDataT, weights_baselineH, a, b, groups) 
{
    cat(".")
    #### baselineH = the selected value for baseline hazard value (the posibility of death within the next infinitsimaly small time step in case gene expression data is 0)
    #### weights = a vector of values to be multiplied by the corresponding gene expression values per patient that determines the extent of influence a gene has on the patients predicted survival time if any
    #### geDataT = a dataframe/matrix of values of gene expression values for each gene on the microarray chip, for each patient. colnames = Gene Name/Clone ID|Data, rownames = patient number
    #### survDataT = a dataframe with the observed survival time of the patients and the censorship status, 1(T) or 0(F). Header = patient|True_STs|censored
    #### q, s = the parameters affecting the distribution of the prior on the weights

    # Assign variables:
     ooo <- order(survDataT$True_STs) 
             True_STs <- survDataT$True_STs[ooo]
             survDataT <- survDataT[ooo, ]
             geDataT <- geDataT[ooo, ]    
             max_ST <- max(True_STs) + 1
             limits <- seq(0, max_ST, length.out=(groups + 1))
     no.patients <- length(weights_baselineH[-c(1:groups)])
    # allocate return vector
    deriv <- rep(0,length(weights_baselineH))
    # compute length of blh interval
    delta <- limits[length(limits)] / groups
    # compute position of weights
    startw <- groups + 1
    endw <- length(weights_baselineH)
    # exp-transform baseline hazard
    blh <- exp(weights_baselineH[1:groups])
    # compute derivative with respect to prior on blh
    deriv[1:groups] <- (1-a) + blh/b
    # derivative with respect to prior on weights
 #   fac <- s^q
 #   deriv[startw:endw] <- (abs(weights_baselineH[startw:endw])^(q-1))/fac*sign(weights_baselineH[startw:endw])
 #   deriv[is.na(deriv[startw:endw])] <- 0
    # add the derivatives of the likelihood
    summa <- rep(0,(groups+1))
    for (i in 1:groups)
      summa[i+1] <- summa[i] + delta * blh[i]
    # loop over all patients
    for (i in 1:no.patients)
      for (j in 1:groups)
      {
        # is patient i censored in interval j?
        # if not, there is nothing to do
        a <- limits[j]
        b <- limits[j+1]
        t <- True_STs[i]
        lambda <- blh[j]
        fx <- 0.0
        fx <- sum(weights_baselineH[startw:endw]*geDataT[i,])
        if ((t>a)&&(t<=b)){
          if (survDataT$censored[i]==0)
            deriv[startw:endw] <- deriv[startw:endw] - geDataT[i,]
          deriv[startw:endw] <- deriv[startw:endw] + geDataT[i,] * exp(fx) * (lambda * (t-a) + summa[j])
        }
        # deriv. w.r.t. blh
        if ((t>a)&&(t<=b))
        {
          if (survDataT$censored[i]==0)
            deriv[j] <- deriv[j] -1;
          deriv[j] <- deriv[j] + exp(fx) * (t-a) * lambda
        }
        else if (t>b)
          deriv[j] <- deriv[j] + exp(fx) * delta * lambda
      }
    return(deriv)
}