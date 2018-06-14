\name{wnl5}
\alias{wnl5}
\title{Old type WinNonlin - Least Square not MLE}
\description{
It performs old type Winnonlin regression.
}
\usage{
wnl5(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjLS)
}
\arguments{
  \item{Fx}{Function for structural model. It should return a vector of the same length to observations.}
  \item{Data}{Data table which will be used in Fx. Fx should access this with \code{e$DATA}.}
  \item{pNames}{Parameter names in the order of Fx arguments}
  \item{IE}{Initial estimates of parameters}
  \item{LB}{Lower bound for \code{optim} function. The default value is 0.}
  \item{UB}{Upper bound for \code{optim} function. The default value is 1e+06.}
  \item{Error}{Error model. One of \code{"POIS"} for Poisson error, \code{"PROP"} for proportional error, and others for additive error model.}
  \item{ObjFx}{Objective function to be minimized. The default is least square function.}
}
\details{
This uses scaled transformed parameters and environment \code{e} internally. Here we do not provide standard error. If you want standard error, use \code{nlr}.
}
\value{
  \item{PE}{Point estimates}
  \item{WRSS}{Weighted Residual Sum of Square}
  \item{run$m}{Count of positive residuals}
  \item{run$n}{Count of negative residuals}
  \item{run$run}{Count of runs of residuals}
  \item{run$p.value}{P value of run test with excluding zero points}
  \item{Objective Function Value}{Minimum value of the objective function}
  \item{AIC}{Akaike Information Criterion}
  \item{SBC}{Schwarz Bayesian Information Criterion}
  \item{Condition Number}{Condition number}
  \item{Message}{Message from \code{optim}.}
  \item{Prediction}{Fitted(predicted) values}
  \item{Residuals}{Residuals}
  \item{Elapsed Time}{Consumed time by minimization}
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
tData = Theoph
colnames(tData) = c("ID", "BWT", "DOSE", "TIME", "DV")

fPK = function(THETA) # Prediction function
{
  DOSE = 320000 # in microgram
  TIME = e$DATA[,"TIME"] # use data in e$DATA

  K  = THETA[1]
  Ka = THETA[2]
  V  = THETA[3]
  Cp  = DOSE/V*Ka/(Ka - K)*(exp(-K*TIME) - exp(-Ka*TIME))
  return(Cp)
}

IDs = unique(tData[,"ID"])
nID = length(IDs)
for (i in 1:nID) {
  Data = tData[tData$ID == IDs[i],]
  Res = wnl5(fPK, Data, pNames=c("k", "ka", "V"), IE=c(0.1, 3, 500))
  print(paste("## ID =", i, "##"))
  print(Res)
}

}
\keyword{Least Square Estimation (Old WinNonlin)}
