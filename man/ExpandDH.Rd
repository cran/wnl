\name{ExpandDH}
\alias{ExpandDH}
\title{Expand Dosing History Table}
\description{
It expands dosing history table.
}
\usage{
ExpandDH(DH, Fo = 1)
}
\arguments{
  \item{DH}{Dosing history table of NONMEM type}
  \item{Fo}{Bioavailability of the first (gut) compartment}
}
\details{
It expands dosing history table of conventional NONMEM data format. It calculate bioavailable amount, then add time points of non-differentiable, e.g. stopping points of infusion. 
}
\value{
Returns expanded dosing history table.
}
\author{Kyun-Seop Bae <k@acr.kr>}

\keyword{Dosing history}
