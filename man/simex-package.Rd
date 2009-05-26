\name{simex-package}
\Rdversion{1.1}
\encoding{latin1}
\alias{simex-package}
\docType{package}

\title{Error or misclassification correction in models using (MC)SIMEX}
\description{
  Package \code{simex} is an implementation of the SIMEX--algorithm by Cook
  and Stephanski and the MCSIMEX--Algorithm by Küchenhof, Mwalili and Lesaffre.

}

\details{
  \tabular{ll}{
    Package: \tab simex\cr
    Type: \tab Package\cr
    Version: \tab 1.4\cr
    Date: \tab 2009-04-19\cr
    License: \tab GPL 2 or above\cr
    LazyLoad: \tab yes\cr
  }
  The package includes first of all the implementation for the SIMEX-- and
  MCSIMEX--Algorithms. Jackknife and asymptotic variance estimation are
  implemented. Various methods and analytic tools are provided for a simple and
  fast access to the SIMEX-- and MCSIMEX--Algorithm.

  Functions \code{simex()} and \code{mcsimex()} can be used on models issued
  from \code{lm()}, \code{glm()} with asymtotic estimation.
  Models from \code{nls()}, \code{gam()} (package \pkg{mgcv}),
  \code{lme()} and \code{nlme()} (package \pkg{nlme}) can also be corrected with
  these algorithms, but without asymptotic estimations.
}

\author{
Wolfgang Lederer

Maintainer: Wolfgang Lederer,\email{wolfgang.lederer@googlemail.com}
}

\references{
Lederer, W. and Küchenhoof, H. (2006) A short introduction to the SIMEX and
  MCSIMEX. \emph{R News}, \bold{6/4}, 26 -- 31
}

\seealso{
  \code{\link[simex]{simex}}, \code{\link[simex]{mcsimex}}, \code{\link[simex]{misclass}}

  and for functions generating the initial naive models:
  \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{nls},
  \code{gam}, \code{lme}, \code{nlme}
}

\examples{
# See example(simex) and example(mcsimex)
}

\keyword{package}
\keyword{models}
\keyword{regression}
