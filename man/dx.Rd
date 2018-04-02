\name{dx}
\alias{dx}
\title{Simplest diagnostic plot for minimization result}
\description{
It performs simplest diagnostic plot from the result of \code{nlr}.
}
\usage{
dx(r)
}
\arguments{
  \item{r}{a result from \code{nlr} or \code{wnl5}}
}
\details{
This plots 'Observation vs. Prediction' and 'Normalized Redisual vs. Prediction' only. Normalized residual are meant to be distributed as standard normal distribution, N(0, 1).
}
\value{
This just draws plot.
}
\author{Kyun-Seop Bae <k@acr.kr>}

\keyword{Diagnostic Plot}
