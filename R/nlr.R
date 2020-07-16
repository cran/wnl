nlr = function(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjDef, SecNames, SecForms, Method="L-BFGS-B", Sx)
{
#  e = new.env(parent=globalenv()) # environment should exist
  t1 = Sys.time()
  e$Fx = Fx # function of structural model. Fx should return a vector of the same length to Y
  if (toupper(Error)=="S") e$Sx = Sx # Scale (inverse weight) function. Sx should return a matrix of the same length rows to Y
  e$DATA = Data # Fx may use this
  if ("DV2" %in% colnames(Data)) {
    e$Y = Data[,c("DV1","DV2")]
    e$nRec = nrow(e$Y)
  } else { 
    if (!("DV" %in% colnames(Data))) stop("Data should have 'DV' column.")
    if (mean(abs(Data[,"DV"])) < 1e-6 | mean(abs(Data[,"DV"])) > 1e6) warning("DV is too large or too small. Rescale it.")  
    e$Y  = Data[,"DV"] # Observation values, Data should have "DV" column.
    e$nRec = length(e$Y)
    if (sum(is.na(Data[,"DV"])) > 0) stop("DV column should not have NAs.")
  }

  if (length(pNames) != length(IE)) stop("pNames and IE should match.")

  e$pNames = pNames # parameter names in the order of Fx arguments
  e$IE = IE # initial estimate of Fx arguments
  e$nTheta = length(IE)
  e$Error = UT(Error) # BasicUtil fx toupper and Trim
  e$AddErrVar = min(1e5, (min(e$Y[e$Y > 0])/4)^2) # initial estimate of addtive error varaince
  e$PoisErrVar = 1    # initial estimate of proportional error varaince
  e$PropErrVar = 0.1  # initial estimate of proportional error varaince
#  e$PowErrPow = 0.5   # initial estimate of power error power
#  e$PowErrVar = 0.1   # initial estimate of power error varaince
  e$ScaleErrVar = 1   # initial estimate of scale (inverse weight) vector of error varaince
  e$Obj = ObjFx
#  e$FASD = rep(FASD, e$nRec) # Fixed Additive Error Variance Vector

  if (e$Error == "A") {
    e$nEps = 1
    e$IE = c(e$IE, e$AddErrVar)
    e$pNames = c(e$pNames, "AddErrVar")
  } else if (e$Error == "POIS") {
    e$nEps = 1
    e$IE = c(e$IE, e$PoisErrVar)
    e$pNames = c(e$pNames, "PoisErrVar")
  } else if (e$Error == "P") {
    e$nEps = 1
    e$IE = c(e$IE, e$PropErrVar)
    e$pNames = c(e$pNames, "PropErrVar")
  } else if (e$Error == "C") {
    e$nEps = 2
    e$IE = c(e$IE, e$AddErrVar, e$PropErrVar)
    e$pNames = c(e$pNames, "AddErrVar", "PropErrVar")
#  } else if (e$Error == "CFA") { # Additive Fixed Proportional
#    e$nEps = 1
#    e$IE = c(e$IE, e$PropErrVar)
#    e$pNames = c(e$pNames, "PropErrVar")
  } else if (e$Error == "S") {
    e$nEps = 1
    e$IE = c(e$IE, e$ScaleErrVar)
    e$pNames = c(e$pNames, "ScaleErrVar")
  } else {
    e$nEps = 0
  }

  e$nPara = e$nTheta + e$nEps # number of parameters
  if (missing(LB)) {
    e$LB = rep(0, e$nPara) # lower bound
  } else {
    e$LB = c(LB, rep(0, e$nEps))
  }
  if (missing(UB)) {
    e$UB = rep(1e+6, e$nPara) # upper bound
  } else {
    e$UB = c(UB, rep(1e+6, e$nEps))
  }

#  if (e$Error == "POW") {
#    e$LB[e$nTheta+1] = 0
#    e$UB[e$nTheta+1] = 1
#  }

  e$alpha  = 0.1 - log((e$IE - e$LB)/(e$UB - e$LB)/(1 - (e$IE - e$LB)/(e$UB - e$LB))) # scaling constant
  if (e$nEps > 0) {
    e$SGindex = (e$nTheta + 1):(e$nTheta + e$nEps) # index of error variance(s)
  } else {
    e$SGindex = 0
#   e$SG1 = SG1
#   e$SG2 = SG2
#   e$Ci = cbind(rep(SG1*SG1, e$nRec), rep(SG2*SG2, e$nRec))  
#   e$SumLogCi = sum(log(e$Ci))
  }
  e$r = optim(rep(0.1, e$nPara), ObjEst, method=Method)
  e$PE = exp(e$r$par - e$alpha)/(exp(e$r$par - e$alpha) + 1)*(e$UB - e$LB) + e$LB

  e$InvCov = hessian(e$Obj, e$PE)/2  # FinalEst from EstStep()
  e$Cov = solve(e$InvCov)
  colnames(e$Cov) = e$pNames
  rownames(e$Cov) = e$pNames

  for (i in 1:e$nPara) {
    if (e$Cov[i,i] < 0) { 
      e$Cov[i,i] = -e$Cov[i,i]
      warning("Negative variance encountered. Variance estimation is unreliable!")
    }
  }
  e$SE = sqrt(diag(e$Cov))
  e$Correl = cov2cor(e$Cov)
  e$EigenVal = eigen(e$Correl)$values

  if (e$SGindex > 0) {  
    e$PE = c(e$PE, sqrt(e$PE[e$SGindex]))
    e$SE = c(e$SE, e$SE[e$SGindex]/2/sqrt(e$PE[e$SGindex])) # Delta method, See Wackerly p484, gr = (x^0.5)' = 0.5x^(-0.5), gr^2 = 1/(4*x)
  }
  e$RSE = e$SE / e$PE * 100
  e$Est = rbind(e$PE, e$SE, e$RSE)

#  colnames(e$Est) = c(e$pNames)
  if(e$Error == "A") {
    colnames(e$Est) = c(e$pNames, "AddErrSD")
  } else if (e$Error == "POIS") {
    colnames(e$Est) = c(e$pNames, "PoisErrSD")
  } else if (e$Error == "P") {
    colnames(e$Est) = c(e$pNames, "PropErrSD")
  } else if (e$Error == "C") {
    colnames(e$Est) = c(e$pNames, "AddErrSD", "PropErrSD")
#  } else if (e$Error == "POW") {
#    colnames(e$Est) = c(e$pNames, "PowErrSD")
#  } else if (e$Error == "CFA") {
#    colnames(e$Est) = c(e$pNames, "PropErrSD")
  } else if (e$Error == "S") {
    colnames(e$Est) = c(e$pNames, "ScaleErrSD")
  }
  rownames(e$Est) = c("PE", "SE", "RSE")

  if (!missing(SecNames)) {
    nSec = length(SecNames)
    tRes = matrix(nrow=3, ncol=nSec)
    colnames(tRes) = SecNames
    rownames(tRes) = c("PE", "SE", "RSE")
    for (i in 1:nSec) {
      tRes[,i] = t(Secondary(SecForms[[i]], e$Est["PE",1:e$nPara], e$Cov))
    }
    e$Est = cbind(e$Est, tRes)
  }
  e$Pred = e$Fx(e$PE[1:e$nTheta])
  e$Residual = e$Y - e$Pred
  if (Error == "NOSG") {
    e$Residual[e$Y == -1] = 0
    e$nRec = sum(e$Y != -1)
  } 

  e$run = run.test(e$Residual)
  e$'-2LL' = e$nRec*log(2*pi) + e$r$value
  e$AIC = e$'-2LL' + 2*e$nPara
  e$AICc = e$AIC + 2*e$nPara*(e$nPara + 1)/(e$nRec - e$nPara - 1)
  e$BIC = e$'-2LL' + e$nPara*log(e$nRec)
  e$Elapsed = difftime(Sys.time(), t1)
  if (toupper(Error) != "S") {
    Result = list(e$Est, e$Cov, e$run, e$r$value, e$'-2LL', e$AIC, e$AICc, e$BIC, e$r$covergence, e$r$message, e$Pred, e$Residual, e$Elapsed)
    names(Result) = c("Est", "Cov", "run", "Objective Function Value", "-2LL", "AIC", "AICc", "BIC", "Convergence", "Message", "Prediction", "Residual", "Elapsed Time")
  } else {
    Scale = Sx(e$Est["PE", 1:e$nTheta])
    Result = list(e$Est, e$Cov, e$run, e$r$value, e$'-2LL', e$AIC, e$AICc, e$BIC, e$r$covergence, e$r$message, e$Pred, e$Residual, Scale, e$Elapsed)
    names(Result) = c("Est", "Cov", "run", "Objective Function Value", "-2LL", "AIC", "AICc", "BIC", "Convergence", "Message", "Prediction", "Residual", "Scale", "Elapsed Time")    
  }
  return(Result)
}

