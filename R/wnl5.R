wnl5 = function(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjLS)
{
#  e = new.env(parent=globalenv()) # environment should exist
  e$Fx = Fx # function of structural model. Fx should return a vector of the same length to Y
  e$DATA = Data # Fx use this
  e$Y  = Data[,"DV"] # Observation values, Data should have "DV" column.
  e$nRec = length(e$Y)
  e$IE = IE # initial estimate of Fx arguments
  e$nPara = length(IE)
  e$Error = UT(Error) # BasicUtil fx toupper and Trim
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
  e$Residual = e$Y - e$Fx(e$PE)
  e$run.test = run.test(e$Residual)
  e$WRSS = e$r$value
  e$AIC = e$nRec*log(e$WRSS) + 2*e$nPara
  e$SBC = e$nRec*log(e$WRSS) + e$nPara*log(e$nRec)

  names(e$PE) = pNames
  Result = list(e$PE, e$WRSS, e$run.test, e$AIC, e$SBC)
  names(Result) = c("PE", "WRSS", "run", "AIC", "SBC")
  return(Result)
}

