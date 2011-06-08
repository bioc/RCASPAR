####Using SurvivROC, calculate AU(all)ROC curves.
survivAURC <- function(Stime, status, marker, entry=NULL, cut.values=NULL, time.max=20, by=1) {
                  #### Replace the 0s with 1s and vice versa since the Heagerty code deals with censored patients opposite to R-CASPAR
                  #status[status == 0] <- "F"
                  #status[status == 1] <- 0
                  #status[status == "F"] <- 1
                  
                  time.points <- seq(0, time.max, by)
                  AURC <- 0
                  AURC_p <- c()
                  for (t in time.points) {
                      AURC_t <- survivROC(Stime=Stime, status=status, marker=marker, cut.values=cut.values, predict.time=t, plot=FALSE)
                      AURC_p <- c(AURC_p, AURC_t$AUC)
                  }
                  index <- match("NaN", AURC_p, nomatch=0)
                  if (length(index) == 1 & index[1] == 0) {
                      x_axis <- time.points
                      y_axis <- AURC_p 
                  }
                  else {
                      x_axis <- time.points[-index]
                      y_axis <- AURC_p[-index]
                  }
                  AURC <- trapezoid(x_axis, y_axis)
                  plot(x_axis, y_axis, type="s", xlab=paste("time", "\n", "AUC= ", signif(AURC, 3)), ylab="AUeachROC", main="Area Under ROC curves")
                  #return(list(AUC=AURC, AUeachROC=AURC_p))
                  return(list(AUC=AURC, AUeachROC=y_axis))
            }