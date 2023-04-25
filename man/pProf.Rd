\name{pProf}
\alias{pProf}
\title{Plot Likelihood or Objective Fnnction Value Profile}
\description{
  It plots estimated likelihood profile. This is not profile likelihood profile.
}
\usage{
  pProf(Bag = e, Title = "", ...)
}
\arguments{
  \item{Bag}{an environment or an object containing the objects of resultant environment e after nlr()}
  \item{Title}{title for the plot}
  \item{...}{arguments to pass to the plot function}
}
\details{
  This plots likelihood profile from the result of nlr() function. Bag should contain the results of nlr().
}
\value{
  No values but a plot.
}
\author{Kyun-Seop Bae <k@acr.kr>}
