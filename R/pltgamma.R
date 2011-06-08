#### Plot the gamma distribution with the chosen parameters
pltgamma <- function (a, b) {
    x <- seq(0,20,000.1)
    y <- x ^ (a - 1) * exp((-x/b)) / ((b ^ a) * gamma(a)) 
    plot (x, y, type="s")
    }