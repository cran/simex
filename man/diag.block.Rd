\name{diag.block}
\Rdversion{1.1}
\alias{diag.block}

\title{Constructs a block diagonal matrix}

\description{
  The function takes a \code{list} and constructs a block diagonal matrix with
  the elements of the list on the diagonal. If \code{d} is not a list then
  \code{d} will be repeated \code{n} times and written on the diagonal (a wrapper
  for \code{kronecker()})
}

\usage{
diag.block(d, n)
}
\arguments{
  \item{d}{a \code{list} of matrices or vectors, or a matrix or vector}
  \item{n}{number of repetitions}
}

\value{
  returns a matrix with the elements of the list or the repetitions of the
  supplied matrix or vector on the diagonal.
}

\author{ Wolfgang Lederer, \email{wolfgang.lederer@googlemail.com}}

\seealso{\code{\link[base]{diag}}, \code{\link[base]{kronecker}}}

\examples{
a <- matrix(rep(1, 4), nrow = 2)
b <- matrix(rep(2, 6), nrow = 2)
e <- c(3, 3, 3, 3)
f <- t(e)
d <- list(a, b, e, f)
diag.block(d)
diag.block(a, 3)
}

\keyword{datagen}
