\name{build.mc.matrix}
\alias{build.mc.matrix}
\title{Build misclassification matrices from empirical estimations}
\description{
Empirical misclassification matrices to the power of lambda my not exist for small 
values of lambda. This function provides methods to estimate the nearst version of the
misclassification matrix that satisfies the conditions a misclassification matrix has to fulfill
}
\usage{
build.mc.matrix(mc.matrix, method = "series", tuning = sqrt(.Machine$double.eps), diag.cor = FALSE, tol = .Machine$double.eps, max.iter = 100)
}
\arguments{
  \item{mc.matrix}{an empirical misclassification matrix}
  \item{method}{method used to estimate the  generator for the misclassification matrix (see Details)}
  \item{tuning}{security parameter for numerical reasons}
  \item{diag.cor}{should corrections substracted from the diagonal or from all values corresponding to the size}
  \item{tol}{tollerance level for series method for convegence}
  \item{max.iter}{maximal number of iterations for the serie method to converge}
}
\details{
Method "series" constructs a generator via the series 

$(Pi-I) - (Pi-I)^2/2 + (Pi-I)^3/3 - \dots$

Method "log" constructs the generator via taking the log of the misclassification matrix

small negative off diagonal values are corrected and set to (0 + tuning).

The amount used to correct for negative values is added to the diagonal element if \code{diag.cor =T} 
and distributed among all values if \code{diag.cor =F}

Method "jlt" uses the method described by Jarro et al. (see Israel et al.)


}
\value{
 returns a misclassification matrix that is the closest estimate for a working misclassification matrix
}
\references{Israel, Robert B., Rosenthal, Jeffrey S., Wei, Jason Z., Finding Generators for Markov Chains via Empirical Transition Matrices, with Applications to Credit Ratings, Mathematical Finance, 11, 245--265}
\author{Wolfgang Lederer}
\note{Does not always work! So check properly.}
\seealso{\code{\link{check.mc.matrix}},\code{\link{mcsimex}}}
\examples{
Pi <- matrix(data = c(0.989,0.01,0.001,0.17,0.829,0.001,0.001,0.18,0.819), nrow =3, byrow =FALSE) 
check.mc.matrix(list(Pi))
check.mc.matrix(list(build.mc.matrix(Pi)))
build.mc.matrix(Pi)


pi3 <- matrix(c(0.8,0.2,0,0,0,0.8,0.1,0.1,0,0.1,0.8,0.1,0,0,0.3,0.7), nrow =4)
check.mc.matrix(list(pi3))
build.mc.matrix(pi3)
check.mc.matrix(list(build.mc.matrix(pi3)))

}
\keyword{regression}
