\name{SolComp2}
\alias{SolComp2}
\title{Get Lambdas and Coefficients of Two-compartment Model}
\description{
It calculates lambdas and coefficients for two-compartment model from K10, K12, and K21.
}
\usage{
SolComp2(K10, K12, K21)
}
\arguments{
  \item{K10}{Ke, Elimination rate constant from central compartment}
  \item{K12}{Rate constant from the central to the peripheral compartment}
  \item{K21}{Rate constant from the peripheral to the central compartment}
}
\details{
It calculates lambdas and coefficients of two-compartment model from K10, K12, and K21. Lambdas should have no identical values.
}
\value{
This returns a list of lambdas and coefficients.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
DAT
DAT2 = ExpandDH(DAT)
Sol = SolComp2(K10=0.1, K12=3, K21=1)
X2 = nComp(Sol, Ka=1, DAT2)
X2
matplot(DAT2[, "TIME"], X2, type="l")
}
\keyword{Two-compartment}
