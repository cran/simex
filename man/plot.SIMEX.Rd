\name{plot.SIMEX}
\alias{plot.SIMEX}
\alias{plot.MCSIMEX}
\title{Plot method for the SIMEX- and MCSIMEX-class}
\description{
 Creates plots for the SIMEX- or MCSIMEX-class. For each variable a plot is produced showing the results of the extrapolation step.
}
\usage{
plot.SIMEX(x, xlab = expression((1 + lambda)), ylab = colnames(b[, -1])
	, ask = FALSE, show = rep(TRUE, NCOL(b) - 1), ...)
}
\arguments{
  \item{x}{object of class SIMEX}
  \item{xlab}{optional name for the X-Axis}
  \item{ylab}{vector containing the names for the Y-Axis}
  \item{ask}{logical. If \code{TRUE}, the user is asked for input, before a new figure is drawn.}
  \item{show}{vector of logicals indicating for wich variables a plot should be produced}
  \item{\dots}{optoinal, passed to \code{par()}}
}
\author{Wolfgang Lederer}
\seealso{ \code{\link{simex}},\code{\link{mcsimex}},\code{\link[graphics]{par}}}
\examples{
a <- rnorm(100)
b <- a*20 + rnorm(100,100,10)
e <- a + rnorm(100,1,5)
naive.model <- lm(a~ b + e, x= TRUE,y =TRUE)

w<- simex(model = naive.model
	, SIMEXvariable = c("b","e")
	, measurement.error = c(10,5)
	, lambda=c(0.5,1,1.5,1.75,2.25,2.5,3))

plot(w, ask = FALSE, mfrow = c(2,2))
}
\keyword{models}
