\name{Comp1}
\alias{Comp1}
\title{One compartment model - analytical}
\description{
It calculates using one compartment model.
}
\usage{
Comp1(Ke, Ka=0, DH)
}
\arguments{
  \item{Ke}{Elimination rate constant}
  \item{Ka}{Absorption rate constant}
  \item{DH}{Expanded dosing history table}
}
\details{
First compartment is the gut compartment for oral dosing. IV bolus and infusion dosing should be done at the second compartment.
}
\value{
This returns a table with the gut and the central compartment columns
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
DAT
DAT2 = ExpandDH(DAT)
X1 = Comp1(Ke=0.1, Ka=1, DAT2)
X1
matplot(DAT2[, "TIME"], X1, type="l")
}
\keyword{One compartment}
