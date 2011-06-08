#### Plotting a Kaplan Meier plot
#### Input is expeceted to be a table with at least the domains Patient, survival.time, censored (censorship status), title=title for the generated graph
#### Censorship status can be either False or True.
kmplt_svrl <- function (all=NULL, long, short, title) {
                km.svrl<-function(data) {
                         no.subjects <- nrow(data)
                         new.row <- rep(NA,ncol(data))
                         data.new <- rbind(new.row, data)
                         data.new$True_STs[1] <- 0
                         data.new$censored[1] <- TRUE
                         data.new <- data.new[order(data.new$True_STs), ]
                         data.new$censored[data.new$censored == TRUE] <- 1
                         #### ri=patients at risk and mi=patients who die
                         risk <- c(0,no.subjects:1)
                         deaths <- c()
                         pt<-c(NA) #### proportion of patinets surviving at time t
                         St<-c(1) #### cumulitave proporiton of patients surviving at time t;S(t)
                         for(death in 1:(no.subjects+1)){
                             if(data.new$censored[death] == 1) 
                                 mi <- 0
                             else
                                 mi <- 1
                             deaths <- c(deaths, mi)
                         }
                         data.new <- cbind(data.new, risk, deaths)
                         for(patient in 2:(no.subjects + 1)) {
                             pi <- (risk[patient]-deaths[patient])/risk[patient]
                             St <- c(St, St[patient - 1]*pi)
                             pt <- c(pt, pi)
                         }
                         data.new <- cbind(data.new, pt, St)
                         y_axis <- c(data.new$St[length(data.new$True_STs)])
                         x_axis <- data.new$True_STs
                         for(survival in c(no.subjects:1)) {
                             if(data.new$True_STs[survival] != data.new$True_STs[survival+1]) 
                                 y_axis <- c(y_axis, data.new$St[survival])
                             else y_axis <- c(y_axis, y_axis[length(y_axis)])
                         }
                         y_axis <- y_axis[length(y_axis):1]
                         return(list(x=x_axis, y=y_axis))#, type="S", xlab="Time", ylab="Probability of survival S(t)", main=title))
                }
                if(length(all)!=0) {
                    required1 <-km.svrl(all)
                    plot(x=required1$x, y=required1$y,type="s",xlab="time",ylab="S(t)",main=title, ylim=0:1)
                    required2 <- km.svrl(long)
                    lines(x=required2$x, y=required2$y, type="s", col="forestgreen")
                    required3 <- km.svrl(short)
                    lines(x=required3$x,y=required3$y, type="s", col="red")
                    legend("topright", legend=c("all","long","short"), cex=1, col=c("black","forestgreen","red"), lwd=2, bty="n")
                }
                else {
                    required1 <-km.svrl(long)
                    plot(x=required1$x, y=required1$y,type="s",xlab="time",ylab="S(t)",main=title, ylim=0:1, col="forestgreen")
                    required2 <- km.svrl(short)
                    lines(x=required2$x,y=required2$y, type="s", col="red")
                    legend("topright", legend=c("long","short"), cex=1, col=c("forestgreen","red"), lwd=2, bty="n")
                    }
}