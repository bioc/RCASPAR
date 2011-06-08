#### Slightly modified version of P.Heagerty's survivalROC function.
#### Stime=True survival time
#### status=1/T if patient died, 0/F if patient is censored
#### marker=Predicted survival time
#### entry=time of entry into the study
#### predicted.time= time point at which the ROC analysis is to be performed
#### cut.values=the values at which the sensitivity and specificity are evaluated
survivROC <- function (Stime, status, marker, entry= NULL, predict.time, cut.values=NULL, plot=TRUE) {
    #### Replace the 0s with 1s and vice versa since the Heagerty code deals with censored patients opposite to R-CASPAR
    status[status == 0] <- "F"
    status[status == 1] <- 0
    status[status == "F"] <- 1
    method = "KM"
    times = Stime                                                                                                                                                         
    x = marker     
    if (is.null(entry))                                                                                                                                                   
        entry <- rep(0, length(times))                                                                                                                                    
    bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)                                                                                                         
    entry <- entry[!bad]                                                                                                                                                  
    times <- times[!bad]                                                                                                                                                  
    status <- status[!bad]                                                                                                                                                
    x <- x[!bad]                                                                                                                                                          
    if (sum(bad) > 0)                                                                                                                                                     
        cat(paste("\n", sum(bad), "records with missing values dropped. \n"))                                                                                             
    if (is.null(cut.values)) {                                                                                                                                              
        cut.values <- c()
        PredictedSTs <- unique(x)
        PredictedSTs <- c(PredictedSTs, (max (PredictedSTs) + 1))
        PredictedSTs <- PredictedSTs[order (PredictedSTs)]
        for (ST in 1:length (x)) {
            cut.values <- c(cut.values, mean(PredictedSTs[ST:(ST+1)]))
        }
    }                                                                                                                           
    ncuts <- length(cut.values)                                                                                                   
    ooo <- order(times)                                                                                                                                                   
    times <- times[ooo]                                                                                                                                                   
    status <- status[ooo]                                                                                                                                                 
    x <- x[ooo]                                                                                                                                                           
    s0 <- 1                                                                                                                                                               
    unique.t0 <- unique(times)                                                                                                                                            
    unique.t0 <- unique.t0[order(unique.t0)]                                                                                                                              
    n.times <- sum(unique.t0 <= predict.time)                                                                                                                            
    for (j in 1:n.times) {                                                                                                                                                
        n <- sum(entry <= unique.t0[j] & times >= unique.t0[j])                                                                                                           
        d <- sum((entry <= unique.t0[j]) & (times == unique.t0[j]) &                                                                                                      
            (status == 1))                                                                                                                                                
        if (n > 0)                                                                                                                                                        
            s0 <- s0 * (1 - d/n)                                                                                                                                          
    }                                                                                                                                                                     
    s.pooled <- s0                                                                                                                                                        
    roc.matrix <- matrix(NA, ncuts, 2)                                                                                                                                    
    roc.matrix[ncuts, 1] <- 0                                                                                                                                             
    roc.matrix[ncuts, 2] <- 1                                                                                                                                             
    if (method == "KM") {                                                                                                                                                 
        for (c in 1:(ncuts - 1)) {                                                                                                                                        
            s0 <- 1                                                                                                                                                       
            subset <- as.logical(x > cut.values[c])
            e0 <- entry[subset]                                                                                                                                           
            t0 <- times[subset]                                                                                                                                           
            c0 <- status[subset]                                                                                                                                          
            if (!is.null(t0)) {                                                                                                                                           
                unique.t0 <- unique(t0)                                                                                                                                
                unique.t0 <- unique.t0[order(unique.t0)]                                                                                                                  
                n.times <- sum(unique.t0 <= predict.time)                                                                                                              
                if (n.times > 0) {                                                                                                                                        
                  for (j in 1:n.times) {                                                                                                                                  
                    n <- sum(e0 <= unique.t0[j] & t0 >= unique.t0[j])                                                                                                     
                    d <- sum((e0 <= unique.t0[j]) & (t0 == unique.t0[j]) &                                                                                                
                      (c0 == 1))                                                                                                                                          
                    if (n > 0)                                                                                                                                            
                      s0 <- s0 * (1 - d/n)                                                                                                                                
                  }                                                                                                                                                       
                }                                                                                                                                                         
            }                                                                                                                                                             
            p0 <- mean(subset)                                                                                                                                        
            roc.matrix[c, 1] <- (1 - s0) * p0/(1 - s.pooled)                                                                                                              
            roc.matrix[c, 2] <- 1 - s0 * p0/s.pooled                                                                                                                      
        }                                                                                                                                                                 
    }                                                                                                                                                                     
    sensitivity = roc.matrix[, 1]
    specificity = roc.matrix[, 2]
    y <- 1 - c(0, specificity)
    x <- c(1, sensitivity)
    n <- length(x)
    dx <- x[-n] - x[-1]
    mid.y <- (y[-n] + y[-1])/2
    area <- sum((dx) *( mid.y))
    if(plot == TRUE)
       plot(x, y, type="s", xlab=paste("1-Specificity", "\n", "AUC= ", signif(area, 3)), ylab="Sensitivity", main=paste("ROC at t=",predict.time,"years"))
    list(cut.values = c(-Inf, cut.values), Comp.Specificity = y, Sensitivity = x, predict.time = predict.time,
        Survival = s.pooled, AUC = area)
}
