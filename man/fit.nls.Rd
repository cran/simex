\name{fit.nls}
\alias{fit.nls}
\title{Fits nls-models for the extrapolation step of SIMEX models}
\description{Fits a nonlinear model of the form a+b/(c+lambda) to the data. Starting values are obtained 
by unsing quadratic extrapolation and an analytic way to obtain a three point fit to the model.

This function is for internal use only!}
\usage{fit.nls(lambda, p.names, estimates)}

\arguments{
  \item{lambda}{ vector of lambdas, for which the extrapolation should be made}
  \item{p.names}{ vector of variable names }
  \item{estimates}{ matrix of Simex estimates, as supplied by the functions \code{simex}}
}
\details{
Uses \code{nls()} with the obtained starting values as described by Carrol, Ruppert, Stefanski.
}
\value{a list of \code{nls}-objects, names are \code{p.names}
}
\author{Wolfgang Lederer}
\note{numerical not stable, especially for constant functions}
\seealso{\code{\link{simex}}, \code{\link[stats]{nls}}}
\keyword{models}
