\name{construct.s}
\alias{construct.s}
\title{Function to construct the s(gamma)-Matrix }
\description{
  This function constructs the matrix of the first derivatives of the extrapoltion function for the value -1. 
  In Carrol, Küchenhoff, Lombard and Stefanski the resulting matrix is called s(gamma). 
  
  This function is for internal use only!
}
\usage{
construct.s(ncoef, lambda, fitting.method, extrapolation = NULL)
}
\arguments{
  \item{ncoef}{ number of coefficients in the naive model }
  \item{lambda}{ lambda (including 0) }
  \item{fitting.method}{ fitting method used in the extrapolation step }
  \item{extrapolation}{ the extrapolation model, only necessary if nonlinear or log-linear extrapolation is used}
}
\value{
	returns the matrix s(gamma)
}
\references{ Carrol, Küchenhoff, Lombard and Stefanski, 1996, \emph{Asymptotics for the SIMEX Estimator in nonlinear Measurement Error Models, Journal of the American Statistical Association}, \bold{91}, 242--250 }
\author{ Wolfgang Lederer }
\seealso{ \code{\link{simex}},\code{\link{mcsimex}}}
\examples{
construct.s(2,c(0,0.5,1,1.5,2), "quad")
}
\keyword{regression}

