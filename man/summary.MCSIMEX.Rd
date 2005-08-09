\name{summary.MCSIMEX}
\alias{summary.MCSIMEX}
\alias{print.summary.MCSIMEX}
\title{Summarizing Model Fits, corrected by the MCSIMEX-method}
\description{
summary method for class MCSIMEX.
}
\usage{
summary.MCSIMEX(object,...)
}
\arguments{
  \item{object}{object of class MCSIMEX}
  \item{\dots}{further arguments}
}
\value{
  \item{coefficients}{a p x 4 matrix with columns for the estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
  \item{residuals}{residuals of the corrected model}
  \item{call}{the function call}
  \item{B}{number of iterations}
  \item{naive.model}{the naive model}
  \item{SIMEXvariable}{character vector of the SIMEXvariable}
  \item{mc.matrix}{The misclassification matrices}
  }
\author{Wolfgang Lederer}
\seealso{\code{\link{mcsimex}}}
\examples{
  # continues example from mcsimex()
  # summary(simex.model)
}
\keyword{regression}

