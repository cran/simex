\name{predict.MCSIMEX}
\alias{predict.MCSIMEX}
\title{ Predict method for the MCSIMEX class }
\description{
  Predicts values for objects of class MCSIMEX
  }
\usage{
predict.MCSIMEX(object, newdata, ...)
}
\arguments{
  \item{object}{Object of class MCSIMEX}
  \item{newdata}{optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.}
  \item{\dots}{additional arguments passed the original predict method}
}
\details{
  Just a wrapper for various predict methods.
}
\value{
  data.frame of predicted values
}
\author{Wolfgang Lederer}
\seealso{\code{\link[stats]{predict}},\code{\link{mcsimex}}}
\keyword{models}
