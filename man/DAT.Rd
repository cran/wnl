\name{DAT}
\alias{DAT}
\docType{data}
\title{
An Example of Dosing History Table
}
\description{
This is a conventional NONMEM input data format. 
}
\usage{RptCfg}
\format{
  This data frame has 5 columns with 18 time-points for the simulation.
  \describe{
    \item{\code{TIME}}{Time}
    \item{\code{AMT}}{Amount given for the compartment of \code{CMT} column}
    \item{\code{RATE}}{Infusion rate}
    \item{\code{CMT}}{Compartment number, 1=gut, 2=central, 3=peripheral, etc.}
    \item{\code{DV}}{Currently blank and not used.}
  }
}
\details{
To be used at \code{Comp1} or \code{nComp}, expand dosing history with \code{ExpandDH} function. 
}

\keyword{datasets}
