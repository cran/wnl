\name{pProf}
\alias{pProf}
\title{Plot Likelihood or Objective Fnnction Value Profile}
\description{
  It plots likelihood profile.
}
\usage{
  pProf(Bag = e, Title = "", ...)
}
\arguments{
  \item{Bag}{an environment of a list}
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
