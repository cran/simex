\name{predict.SIMEX}
\alias{predict.SIMEX}
\title{Predict Method}
\description{
 Predicts values for objects of class SIMEX
}
\usage{
predict.SIMEX(object, newdata, ...)
}
\arguments{
  \item{object}{Object of class SIMEX}
  \item{newdata}{optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.}
  \item{\dots}{additional arguments passed the original predict method}
}
\details{
  Just a wrapper for various predict methods
}
\value{
    data.frame of predicted values
}
\author{Wolfgang Lederer}
\seealso{\code{\link[stats]{predict}},\code{\link{simex}}}
\keyword{regression}
