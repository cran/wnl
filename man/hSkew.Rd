\name{hSkew}
\alias{hSkew}
\title{Hougaard Measure of Skewness}
\description{
  Hougaard measure of skewness with nonlinear regression
}
\usage{
  hSkew(rx)
}
\arguments{
  \item{rx}{a result of nls function}
}
\details{
  Hougaard measure of skewness can be used to check if the parameters of nonlinear regression behavior in linear fashion, i.e. symmetric confidence interval. Be cautious on the variable name conflict. All the variables in the nonlinear function should be able to be accessed by the function.
}
\value{
  Hougaard estimate of skewness for each parameter
  \item{(0, 0.1]}{The estimate is very close-to-linear in behavior.}
  \item{(0.1, 0.25]}{The estimate is reasonably close-to-linear in behavior.}
  \item{(0.25, 1]}{The skweness is apparent.}
  \item{>1}{The estimate is considerably nonlinear in behavior.}
}
\examples{
  r1 = nls(density ~ b1*conc/(b2 + conc), DNase[DNase$Run == 1, ], start=list(b1=3, b2=5))
  hSkew(r1)
}
\references{
  EL-Shehawy SA. On Calculating the Hougaard Measure of Skewness in a Nonlinear Regression Model with Two Parameters. J Math & Stat. 2009;5(4):360-364.
}
\author{Kyun-Seop Bae k@acr.kr}
