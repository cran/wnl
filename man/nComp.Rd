\name{nComp}
\alias{nComp}
\title{Get Amounts of Each Compartments using Lambdas and Coefficients of Multi-compartment Model}
\description{
It calculates using multi-compartment model.
}
\usage{
nComp(Sol, Ka=0, DH)
}
\arguments{
  \item{Sol}{Solution list of lambdas and coefficients}
  \item{Ka}{Absorption rate constant}
  \item{DH}{Expanded dosing history table}
}
\details{
First compartment is the gut compartment for oral dosing. IV bolus and infusion dosing should be done at the second compartment.
If a bolus dose was given at time T, it is reflected at times of larger than T. This is more close to real observation. ADAPT does like this, but NONMEM does not.
}
\value{
This returns a table with the gut and the other compartment columns
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
\keyword{Multi-compartment}
