\name{simpson}
\alias{simpson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function that calculates the area under a curve based on the Simposon algorithm
}
\description{
A function that calculates the approximate value of the definite integral of a continuous function. In other words, it can help plot the area under the curve of the plotted 
function between two limits.
}
\usage{
simpson(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of the values at which the function is to be plotted. 
}
  \item{y}{
A vector with the values of the function at the corresponding x-values.
}
}

\value{
A single numerical value of the approximate area under the curve generated with the x and y values.
}
\references{
Hennion, P.E.(1962). Algorithm 84: Simpson's integration. Communications of ACM. 5(4), 208
}
\author{
Douaa Mugahid
}
\note{
Compared to the trapezoidal algorithm, this is usually more accurate.
}


\seealso{
\code{\link{trapezoid}}
}
\examples{
x <- seq(0:20)
y <- seq(0, 100, 1)
simpson(x,y)}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{integration}
\keyword{area under the curve}
