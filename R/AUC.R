#### For the calcualtion of the area under the curve
    #### Trapezoidal rule
    trapezoid <- function(x, y) {
                     sum(diff(x)*(y[-1]+y[-length(y)]))/2
                 }
    #### Simpson's rule: more accurate
    simpson <- function(x, y){
                   n <- length(y)
                   odd <- n %% 2
                   if (odd)
                       area <- 1/3*sum( y[1] + 2*sum(y[seq(3,(n-2),by=2)]) + 4*sum(y[seq(2,(n-1),by=2)]) + y[n])
                   if (!odd)
                       area <- 1/3*sum( y[1] + 2*sum(y[seq(3,(n-3),by=2)]) + 4*sum(y[seq(2,(n-2),by=2)]) + y[n-1]) + 1/2*(y[n-1] + y[n])
                   dx <- x[2] - x[1]
                   return(area * dx)
               }


