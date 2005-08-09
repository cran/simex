\name{extract.covmat}
\alias{extract.covmat}
\title{Extracts the covariance matrix}
\description{
  Extracts the covariance matrix from various models 
  and returns them in a way, so that the \code{simex} and \code{mcsimex} 
  functions can use them.
}
\usage{
extract.covmat(model)
}
\arguments{
  \item{model}{object of a fitting function}
}
\details{
  So far objects of class \code{lm}, \code{glm}, \code{gam}, \code{nls}, \code{lme}, \code{nlme} are supported.
}
\value{
  \item{covmat}{ covariance matrix of the supplied model}
  }
\author{Wolfgang Lederer}
\note{Easily extendable for other objects. Is required for jackknife variance estimation.

This function is for internal use only!}
\seealso{\code{\link{simex}},\code{\link{mcsimex}}}
\examples{
x <- rnorm(100,0,100)
u <- rnorm(100,0,25)
w <- x+u
y <- x +rnorm(100,0,9)
model <- lm(y~w, x = TRUE, y =TRUE)
extract.covmat(model)
summary(model)$cov.unscaled*summary(model)$sigma^2
}
\keyword{models}
