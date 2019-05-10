\name{SolComp3}
\alias{SolComp3}
\title{Analytical Solution for Three-compartment Model}
\description{
It calculates lambdas and coefficients for three-compartment model from K10, K12, K21, K13, and K31.
}
\usage{
SolComp3(K10, K12, K21, K13, K31)
}
\arguments{
  \item{K10}{Ke, Elimination rate constant from central compartment}
  \item{K12}{Rate constant from the central to the first peripheral compartment}
  \item{K21}{Rate constant from the first peripheral to the central compartment}
  \item{K13}{Rate constant from the central to the second peripheral compartment}
  \item{K31}{Rate constant from the second peripheral to the central compartment}
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
Sol = SolComp3(K10=0.1, K12=3, K21=1, K13=2, K31=0.5)
X3 = nComp(Sol, Ka=1, DAT2)
X3
matplot(DAT2[, "TIME"], X3, type="l")
}
\keyword{Three-compartment}
