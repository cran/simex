\name{summary.SIMEX}
\alias{summary.SIMEX}
\alias{print.summary.SIMEX}
\title{Summarizing Model Fits for the SIMEX method}
\description{
summary method for class SIMEX.
}
\usage{
summary.SIMEX(object,...)
}
\arguments{
  \item{object}{Object of class SIMEX}
  \item{\dots}{forther arguments}
}
\value{
  \item{coefficients}{a p x 4 matrix with columns for the estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
  \item{residuals}{residuals of the corrected model}
  \item{call}{the function call}
  \item{B}{number of iterations}
  \item{naive.model}{the naive model}
  \item{SIMEXvariable}{character vector of the SIMEXvariable}
  \item{measurement.error}{vector of the variances of the knowen measurement error}
}
\author{Wolfgang Lederer}
\seealso{\code{\link{simex}}}
\examples{
  # continues example from simex()
  # summary(simex.model)
}
\keyword{regression}
