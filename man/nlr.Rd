\name{nlr}
\alias{nlr}
\title{Nonlinear Regression in R}
\description{
It performs nonlinear regression usually for pharmacokinetic and pharmacodynamic models.
}
\usage{
nlr(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjDef, SecNames, SecForms, 
    Method="L-BFGS-B", Sx, conf.level=0.95, k, fix=0)
}
\arguments{
  \item{Fx}{Function for structural model. It should return a vector of the same length to observations.}
  \item{Data}{Data table which will be used in Fx. Fx should access this with \code{e$DATA}.}
  \item{pNames}{Parameter names in the order of Fx arguments}
  \item{IE}{Initial estimates of parameters}
  \item{LB}{Lower bound for \code{optim} function. The default value is 0.}
  \item{UB}{Upper bound for \code{optim} function. The default value is 1e+06.}
  \item{Error}{Error model. One of \code{"A"} for additive error, \code{"POIS"} for Poisson error, \code{"P"} for proportional error, \code{"C"} for combined error model, \code{"S"} for general error model. With Error="S", \code{Sx} should be provieded.}
  \item{ObjFx}{Objective function to be minimized. The default is maximum likelihood estimation function(-2 log likelihood).}
  \item{SecNames}{Names of secondary parameter estimates}
  \item{SecForms}{Formula to calculate the secondary parameter estimates}
  \item{Method}{\code{"L-BFGS-B"} is default. See \code{optim} for more detail.}
  \item{Sx}{Scale function. This is usually the inverse of weight. It should return the same length(nrow) of Y. When Error="S", Scale function should be provided as \code{Sx}.}
  \item{conf.level}{Confidence level for confidence interval}
  \item{k}{1/k likelihood interval(LI) will be provided. Currently recommended value is exp(qf(1 - alpha, 1, nRec-nPara)/2) + 1.}
  \item{fix}{indices of parameters to fix}
}
\details{
This uses scaled transformed parameters and environment \code{e} internally.
}
\value{
  \item{Est}{Point estimate(PE) with standard error(SE) and relative standard error(RSE)}
  \item{LI}{1/k likelihood interval, at which likelihood drops to 1/k of maximum likelihood. This reflects asymmetry better than confidence interval. This is estimated likelihood interval, not profile likelihood interval.}
  \item{Skewness}{Hougaard's skewness measure. This is printed only with additive error model. See also \code{hSkew}}
  \item{Cov}{Variance-covariance matrix of the objective function at the value of point estimates}
  \item{run$m}{Count of positive residuals}
  \item{run$n}{Count of negative residuals}
  \item{run$run}{Count of runs of residuals}
  \item{run$p.value}{P value of run test with excluding zero points}
  \item{Objective Function Value}{Minimum value of the objective function}
  \item{-2LL}{-2 times log likelihood}
  \item{AIC}{Akaike Information Criterion}
  \item{AICc}{Corrected Akaike Information Criterion}
  \item{BIC}{Schwarz Bayesian Information Criterion}
  \item{Convergence}{Convergence code from \code{optim}}
  \item{Message}{Message from \code{optim}.}
  \item{Prediction}{Fitted(predicted) values}
  \item{Residuals}{Residuals}
  \item{Scale}{Scales with Error="S". Variances for each points are scale vector multiplied by \code{ScaleErrVar} in \code{Est}.}
  \item{Elapsed Time}{Consumed time by minimization}
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
  tData = Theoph
  colnames(tData) = c("ID", "BWT", "DOSE", "TIME", "DV")

  fPK = function(THETA) # Prediction function
  {
    DOSE = 320000 # in microgram
    TIME = e$DATA[, "TIME"] # use data in e$DATA

    K    = THETA[1]
    Ka   = THETA[2]
    V    = THETA[3]

    P  = DOSE/V*Ka/(Ka - K) * (exp(-K*TIME) - exp(-Ka*TIME))
    return(P)
  }

  IDs = unique(tData[,"ID"])
  nID = length(IDs)
  for (i in 1:nID) {
    Data = tData[tData$ID == IDs[i],]
    Res = nlr(fPK, Data, pNames=c("k", "ka", "V"), IE=c(0.1, 3, 500), 
              SecNames=c("CL", "Thalf", "MRT"), SecForms=c(~V*k, ~log(2)/k, ~1/k))
    print(paste("## ID =", i, "##"))
    print(Res)
  }

  # Another example from radioimmunoassay(RIA)
  d1 = data.frame(conc = c(200, 100, 50, 25, 12.5, 6.25, 3.125, 0),
                  DV = c(1.78, 1.5, 1.17, 0.74, 0.51, 0.31, 0.19, 0.04))

  PRED = function(TH) TH[1] + TH[2]*d1$conc^TH[4]/(TH[3]^TH[4] + d1$conc^TH[4])
  Scale = function(TH) 1/(PRED(TH) - (TH[1] + TH[2])/2)^2

  nlr(PRED, d1, pNames=c("R0", "Rmax", "RC50", "Hill"), IE=c(0.1, 3, 50, 1), 
      Error="S", Sx=Scale)
}
\keyword{Maximum Likelihood Estimation}
