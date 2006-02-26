\name{fit.logl}
\alias{fit.logl}
\title{ Fits a log-linear}
\description{
Fits a log-linear model to the data via \code{nls()}
}
\usage{
fit.logl(lambda, p.names, estimates)
}
\arguments{
  \item{lambda}{the vector of lambdas used in the function \code{mcsimex}}
  \item{p.names}{ parameter names of the naive model}
  \item{estimates}{ results of the simulation step done by function \code{mcsimex}}
}
\details{
Uses the function \code{nls} for direct fitting of the model y~exp(gamma.0 + gamma.1 * lambda). 
As starting values the results of a LS-fit to a linear model with a log transformed response are used. 
If \code{nls} does not converge, the model with the staring values is returned
}
\value{
  \item{extrapolation}{list of individual models for the single coefficients}
 }
\references{}
\author{Wolfgang Lederer}

\seealso{ \code{\link{mcsimex}}}
\keyword{models}% at least one, from doc/KEYWORDS
