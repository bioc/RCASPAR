#### Visualization of the prior with respect to 2 weights depending on the prior parameters s,q
pltprior <- function(q, s) {
    theta1 <- seq(-2, 2)
    theta2 <- theta1
    prior.weights <- function(theta1, theta2) {
        prob <- (q^((q-1)/q)) * exp(-(abs(theta1) + abs(theta2))^q / (q * s^q)) / (2 * s * gamma(1/q))
        return(prob)
                     }
    z <- outer(theta1,theta2, prior.weights)
    persp(theta1 , theta2, z, col="lightgreen", theta=30, phi=20, r=50, d=0.1, expand=0.5, ltheta=90, lphi=180, shade=0.75, ticktype="detailed", nticks=5, zlab="p(theta1, theta2)")
              }