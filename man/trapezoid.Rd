\name{trapezoid}
\alias{trapezoid}
\title{
A function that calculates the area under a curve based on the Simposon algorithm
}
\description{
A function that calculates the approximate value of the definite integral of a continuous function. In other words, it can help plot the area under the curve of the plotted 
function between two limits.
}
\usage{
trapezoid(x, y)
}
\arguments{
  \item{x}{
The values to be used along the x-axis while plotting the curve of the function. The \emph{x} in f(x)=y
}
  \item{y}{
The values to be used along the y-axis while plotting the curve of the function. The \emph{y} in the f(x)=y.
}
}
\value{
The area under the curve plotted with the x and y values provided as arguments.
}
\references{
Weisstein, Eric W. "Trapezoidal Rule." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/TrapezoidalRule.html
}
\author{
Douaa Mugahid
}
\note{
Using this method is slightly less accurate than using the simpson integration method
}
\seealso{
\code{\link{simpson}}
}
\examples{
x <- seq(0:20)
y <- seq(0, 100, 1)
trapezoid(x,y)
}

\keyword{ integral }
\keyword{ area under the curve }
