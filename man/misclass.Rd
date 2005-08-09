\name{misclass}
\alias{misclass}
\title{Generates Misclassified Data }
\description{
  Takes a \code{data.frame} and produces misclassifed data. 
  Probabilities for the missclassification are given in the \code{mc.matrix}.
}
\usage{
misclass(data.org, mc.matrix, k)
}
\arguments{
  \item{data.org}{ \code{data.frame} containing the factor variabels. Must be \code{factors}.}
  \item{mc.matrix}{ a \code{list} of matrices giving the probabilities for the misclassification. 
 Names of the \code{list} must correspond to the variable names in \code{data.org}.
 The \code{colnames} must be named according to the factor levels}
  \item{k}{ the exponent for the misclassification matrix }
}
\value{
	\code{data.frame} containing the misclassified variables
}
\author{Wolfgang Lederer}
\seealso{\code{\link{mcsimex}}, \code{\link{check.mc.matrix}}}
\examples{
x1 <- factor(rbinom(100,1,0.5))
x2 <- factor(rbinom(100,2,0.5))

p1 <- matrix(c(1,0,0,1), nrow = 2)
p2 <- matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8), nrow = 3)

colnames(p1) <- levels(x1)
colnames(p2) <- levels(x2)

x <- data.frame(x1 = x1, x2 = x2)
mc.matrix <- list(x1 = p1, x2 = p2)

x.mc <- misclass(data.org = x, mc.matrix = mc.matrix,k = 1)

identical(x[,1],x.mc[,1]) # T
identical(x[,2],x.mc[,2]) # F
}
\keyword{datagen}

