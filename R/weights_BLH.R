####Function to be called to perform the optimization and determine the best set of weights that help give the observed survival times.
####Takes as an argument the chosen Baseline Hazard,gene expression data as returned from function collect.data,survival data (a dataframe with censorship status "censored"|observed survival times"survival.times"|patients in the correct order (that correrponding to the order in which they were fed into the geDataT),values of q,s as parameters on the prior of weights,the initial values of the parameters, the method to be used,if a hessian is to be returned  and other values to be passed to the optimization function

weights_BLH <- function (geDataT, survDataT, q, s, a, b, groups, par, method=c("Nelder-Mead", "L-BFGS-B", "CG", "BFGS", "SANN"), noprior=1, extras=list(), dist=NULL) {
#### in the three dots pass the gradient (gr) argument if the required method is SANN.
    #### no.patients <<- nrow (geDataT)
    #### no.genes <<- ncol (geDataT)
    #### Determine the groups into which the survival times are marked
        #### Reorder the patients and their data according to the ascending True_STs
    ####ooo <- order(survDataT$True_STs)
    ####True_STs <<- survDataT$True_STs[ooo]
    ####survDataT <<- survDataT[ooo, ]
    ####geDataT <<- geDataT[ooo, ]
    #### Determine the groups into which the survival times are marked
    #### max_ST <- max(True_STs) + 1
    #### limits <<- seq(0, max_ST, length.out=(groups + 1))
    #### diff.limits <<- limits[-1] - limits[-(groups + 1)]
    #### Determine the groups into which each patient falls w.r.t. the baselineHs
    ####baselineH_j <<- rep(NA, no.patients)
        #### for(pat in 1:no.patients) {
            #### for(lim in 1:groups) {
                #### if(True_STs[pat] <= limits[lim + 1] && True_STs[pat] > limits[lim]) {
                    #### baselineH_j[pat] <<- lim
                    #### break
                #### }
            #### }
        #### }
    # log-transform starting values of baseline hazard (par[0:groups]
    par[1:groups] <- log(par[1:groups])
    #### defining variables needed for weights function evaluation
        #### first level of calculations: Independent of weights or gene expression data
        #### a.prior <<- 1/(q*s^q)
        #### c.prior <<- 1/(s^q)
        #### e.prior <<- (1 - a)
        #### g.prior <<- 1/b
        cat("---------------Optimizing------------------ \n")
    #### Performing the actual optimization depending on the method specified by the user
    #### If the method specified by the user is "L-BFGS-B", check whether the user defined any boundary values, then if more than a cetain percent of parameters are zero, then run one optimization iteration without the prior.
        if(method == "SANN") {
            gr <- dist 
            noprior <- 0
        }
        if (noprior > 0) { #### make the percentage a variable the user could define
            maxit.defined <- extras$maxit
            extras$maxit <- noprior
            w <- optim (fn=weight_estimator_BLH_noprior, gr=deriv_weight_estimator_BLH_noprior,geDataT=geDataT, survDataT=survDataT, a=a, b=b, groups=groups, par=par, method=method, control=extras)
            w <- w$par
            extras$maxit <- maxit.defined
            required <- optim (fn=weight_estimator_BLH, gr=deriv_weight_estimator_BLH, q=q, s=s, geDataT=geDataT, survDataT=survDataT, a=a, b=b, groups=groups, par=w, method=method, control=extras)
        }
        else
            required <- optim (fn=weight_estimator_BLH, gr=deriv_weight_estimator_BLH, q=q, s=s, survDataT=survDataT, geDataT=geDataT, a=a, b=b, groups=groups, par=par, method=method, control=extras)
        required$par[1:groups] <- exp(required$par[1:groups])
        cat("\n")
        #rm(a.prior, c.prior, e.prior, g.prior, no.patients, no.genes, True_STs, survDataT, geDataT, limits, diff.limits, baselineH_j, pos=1)
        return (required)
        }
        