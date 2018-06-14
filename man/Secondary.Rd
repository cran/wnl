\name{Secondary}
\alias{Secondary}
\title{Get Secondary Parameter Estimates}
\description{
Get standard error and relative standard error (cv) of the secondary parameter estimate
}
\usage{
Secondary(Formula, PE, COV)
}
\arguments{
  \item{Formula}{Formula to calculate the secondary parameter estimate}
  \item{PE}{Point estimates of primary parameters with names}
  \item{COV}{Variance-covariance matrix of primary estimates}
}
\details{
Variables within \code{Formula} should exist in the names of \code{PE} vector.
}
\value{
This returns point estimate, standard error, relative standard error of the secondary parameter estimate.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
tData = Theoph
colnames(tData) = c("ID", "BWT", "DOSE", "TIME", "DV") # Table requires DV column

fPK = function(THETA) # Prediction function
{
  AMT  = 320000 # in microgram
  TIME = e$DATA[,"TIME"]
  V  = THETA[1]
  K  = THETA[2]
  Ka = THETA[3]
  Cp = AMT/V*Ka/(Ka - K)*(exp(-K*TIME) - exp(-Ka*TIME))
  return(Cp)
}
Data = tData[tData$ID == 1,]
Res = nlr(fPK, Data, pNames=c("V", "K", "Ka"), IE=c(30000, 0.1, 2))
Secondary(~V*K, Res$Est["PE",1:e$nPara], Res$Cov)
}
\keyword{Secondary}
