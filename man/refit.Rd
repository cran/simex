\name{refit}
\alias{refit}
\title{Changing the fitting method}
\description{
  Refits the object of class SIMEX or MCSIMEX with a different extrapolation function and returns an object of the same class with the new estimates.
}
\usage{
refit(object, fitting.method = "quad", asymptotic = TRUE, jackknife = "quad")
}
\arguments{
  \item{object}{object of class SIMEX or MCSIMEX}
  \item{fitting.method}{the new fitting method for for the extrapolation step}
  \item{asymptotic}{logical, indicating whether asymptotic variance estimation should be done. (Only possible if done in the original model)}
  \item{jackknife}{specifying the extrapolation method for jackknife variance estimation. Can be set to FALSE if it should not be performed (Only possible if done in the original model)}
}
\value{
  object of class SIMEX or MCSIMEX
}
\author{Wolfgang Lederer}
\seealso{\code{\link{simex}},\code{\link{mcsimex}}}
\examples{
 # continues example from mcsimex()
 # refit(simex.model,"line")
}
\keyword{regression}

