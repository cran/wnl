wnl5 = function(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjLS)
{
#  e = new.env(parent=globalenv()) # environment should exist
  t1 = Sys.time()
  e$Fx = Fx # function of structural model. Fx should return a vector of the same length to Y
  e$DATA = Data # Fx use this
  e$Y = Data[,"DV"] # Observation values, Data should have "DV" column.
  e$nRec = length(e$Y)
  e$IE = IE # initial estimate of Fx arguments
  e$nPara = length(IE)
  e$Error = toupper(trimws(Error))
  e$Obj = ObjFx

  if (missing(LB)) {
    e$LB = rep(0, e$nPara) # lower bound
  } else {
    e$LB = LB
  }
  if (missing(UB)) {
    e$UB = rep(1e+6, e$nPara) # upper bound
  } else {
    e$UB = UB
  }

  e$alpha  = 0.1 - log((e$IE - e$LB)/(e$UB - e$LB)/(1 - (e$IE - e$LB)/(e$UB - e$LB)))
  e$r = optim(rep(0.1, e$nPara), ObjEst, method="L-BFGS-B")
  e$PE = exp(e$r$par - e$alpha)/(exp(e$r$par - e$alpha) + 1)*(e$UB - e$LB) + e$LB
  e$Hess = hessian(e$Obj, e$PE)
  e$Eigen = eigen(e$Hess)
  e$Cond = sqrt(max(e$Eigen$values)/min(e$Eigen$values))
  e$Pred = e$Fx(e$PE)
  e$Residual = e$Y - e$Pred
  e$run.test = run.test(e$Residual)
  e$WRSS = e$r$value
  e$AIC = e$nRec*log(e$WRSS) + 2*e$nPara
  e$SBC = e$nRec*log(e$WRSS) + e$nPara*log(e$nRec)
  e$Elapsed = difftime(Sys.time(), t1)

  names(e$PE) = pNames
  Result = list(e$PE, e$WRSS, e$run.test, e$AIC, e$SBC, e$Cond, e$r$covergence, e$r$message, e$Pred, e$Residual, e$Elapsed)
  names(Result) = c("PE", "WRSS", "run", "AIC", "SBC", "Condition Number", "Convergence", "Message", "Prediction", "Residual", "Elapsed Time")
  return(Result)
}

