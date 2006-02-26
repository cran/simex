\name{check.mc.matrix}
\alias{check.mc.matrix}
\title{ Check the misclassification matrix }
\description{
  Function to check the misclassification matrix for its existance for exponents smaller than 1.
}
\usage{
check.mc.matrix(mc.matrix, tol = .Machine$double.eps)
}
\arguments{
  \item{mc.matrix}{ a list of matrices to be checked}
  \item{tol}{tolerance for very small values, (this valus will be ignored)}
}
\value{
  vector of logicals
}
\author{ Wolfgang Lederer }
\seealso{ \code{\link{mcsimex}},\code{\link{build.mc.matrix}}}
\examples{
p1 <- matrix(c(1,0,0,1), nrow = 2)
p2 <- matrix(c(0.8,0.15,0,0.2,0.7,0.2,0,0.15,0.8), nrow = 3, byrow =TRUE)
p3 <- matrix(c(0.4,0.6,0.6,0.4), nrow = 2)

mc.matrix <- list(p1,p2,p3)

check.mc.matrix(mc.matrix) # T F F
}
\keyword{ regression }
