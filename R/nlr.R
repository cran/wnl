nlr = function(Fx, Data, pNames, IE, LB, UB, Error="A", ObjFx=ObjDef, SecNames, SecForms, Method="L-BFGS-B", Sx, k=20)
{
#  e = new.env(parent=globalenv()) # environment should exist
  t1 = Sys.time()
  e$Fx = Fx # function of structural model. Fx should return a vector of the same length to Y
  if (toupper(Error) == "S") e$Sx = Sx # Scale (inverse weight) function. Sx should return a matrix of the same length rows to Y
  e$DATA = Data # Fx may use this
  if ("DV2" %in% colnames(Data)) {
    e$Y = Data[, c("DV1", "DV2")]
    e$nRec = nrow(e$Y)
  } else {
    if (!("DV" %in% colnames(Data))) stop("Data should have 'DV' column.")
    if (mean(abs(Data[,"DV"])) < 1e-6 | mean(abs(Data[,"DV"])) > 1e6) warning("DV is too large or too small. Rescale it.")
    e$Y  = Data[, "DV"] # Observation values, Data should have "DV" column.
    e$nRec = length(e$Y)
    if (sum(is.na(Data[, "DV"])) > 0) stop("DV column should not have NAs.")
  }

  if (length(pNames) != length(IE)) stop("pNames and IE should match.")

  e$pNames = pNames # parameter names in the order of Fx arguments
  e$IE = IE # initial estimate of Fx arguments
  e$nTheta = length(IE)
  e$Error = toupper(trimws(Error))
  e$AddErrVar = min(1e5, (min(e$Y[e$Y > 0])/4)^2) # initial estimate of addtive error variance
  e$PoisErrVar = 1    # initial estimate of proportional error variance
  e$PropErrVar = 0.1  # initial estimate of proportional error variance
#  e$PowErrPow = 0.5   # initial estimate of power error power
#  e$PowErrVar = 0.1   # initial estimate of power error variance
  e$ScaleErrVar = 1   # initial estimate of scale (inverse weight) vector of error variance
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
  e$Cov = qr.solve(e$InvCov)
  colnames(e$Cov) = e$pNames
  rownames(e$Cov) = e$pNames

  for (i in 1:e$nPara) {
    if (e$Cov[i, i] < 0) {
      e$Cov[i, i] = -e$Cov[i, i]
#      e$Cov[i, i] = 0
      warning("Negative variance encountered. Variance estimation is unreliable!")
    }
  }
  e$SE = sqrt(diag(e$Cov))
  e$Correl = cov2cor(e$Cov)
  e$EigenVal = eigen(e$Correl)$values

  if (e$SGindex[1] > 0) {
    e$PE = c(e$PE, sqrt(e$PE[e$SGindex]))
    e$SE = c(e$SE, e$SE[e$SGindex]/2/sqrt(e$PE[e$SGindex])) # Delta method, See Wackerly p484, gr = (x^0.5)' = 0.5x^(-0.5), gr^2 = 1/(4*x)
  }
  e$RSE = e$SE/e$PE*100
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
      tRes[, i] = t(Secondary(SecForms[[i]], e$Est["PE", 1:e$nPara], e$Cov))
    }
    e$Est = cbind(e$Est, tRes)
  }
  e$Pred = e$Fx(e$PE[1:e$nTheta])
  e$Residual = e$Y - e$Pred
  if (Error == "NOSG") {
    e$Residual[e$Y == -1] = 0
    e$nRec = sum(e$Y != -1)
  }

## Hougaard Skewness
  if (e$Error == "A") {
    e$J = nGradient(e$Fx, e$PE[1:e$nTheta])
    e$H = nHessian(e$Fx, e$PE[1:e$nTheta])
    e$HouSkew = Hougaard(e$J, e$H, e$PE[e$SGindex])
  }

## Likelihood Profile
  cAdd = e$nRec*log(2*pi)
  nRes = 51
  Alpha = 0.05
  zCut = qnorm(1 - Alpha/2) # 1.959964
  e$zCI = rbind(e$PE - zCut*e$SE, e$PE + zCut*e$SE)[, 1:e$nPara]
  rownames(e$zCI) = c("Lower", "Upper")

  tCut = qt(1 - Alpha/2, e$nRec - e$nPara)
  e$tCI = rbind(e$PE - tCut*e$SE, e$PE + tCut*e$SE)[, 1:e$nPara]
  rownames(e$tCI) = c("Lower", "Upper")

  #e$fCut = qf(1 - Alpha, 1, e$nRec - e$nPara) ; e$fCut # PE change -> SE, sigma change
  #e$fCut = 2*exp(zCut^2/2) ; e$fCut
  #e$fCut = 2*exp(tCut^2/2) ; e$fCut
  e$fCut = 2*log(k) # 1/k likelihood interval

  e$mOFV = matrix(nrow=nRes, ncol=e$nPara)
  e$mPar = matrix(nrow=nRes, ncol=e$nPara)
  for (j in 1:e$nPara) {
    for (k1 in (1:9)/2) {
      tPar = e$PE
      tPar[j] = e$PE[j] - k1*zCut*e$SE[j]
      tOfv = e$Obj(tPar) - e$r$value - e$fCut
      if (is.finite(tOfv) & tOfv > 0) break
    }
    for (k2 in (1:50)/2) {
      tPar = e$PE
      tPar[j] = e$PE[j] + k2*zCut*e$SE[j]
      tOfv = e$Obj(tPar) - e$r$value - e$fCut
      if (is.finite(tOfv) & tOfv > 0) break
    }

    e$mPar[, j] = seq(e$PE[j] - k1*zCut*e$SE[j], e$PE[j] + k2*zCut*e$SE[j], length.out=nRes)
    for (i in 1:nRes) {
      tPar = e$PE
      tPar[j] = e$mPar[i, j]
      e$mOFV[i, j] = cAdd + e$Obj(tPar) # -2LL
    }
  }

## Lihood Interval (LI)
  mSign = sign(e$mOFV - cAdd -e$r$value - e$fCut)
  e$mParB = matrix(nrow=4, ncol=e$nPara)
  for (j in 1:e$nPara) {
    if (min(which(mSign[, j] < 0)) > 1) {
      e$mParB[1, j] = e$mPar[min(which(mSign[, j] < 0)) - 1, j]
      e$mParB[3, j] = 1
    } else {
      e$mParB[1, j] = e$mPar[1, j]
      e$mParB[3, j] = 0
    }
    if (max(which(mSign[, j] < 0)) < nRes) {
      e$mParB[2, j] = e$mPar[max(which(mSign[, j] < 0)) + 1, j]
      e$mParB[4, j] = 1
    } else {
      e$mParB[2, j] = e$mPar[nRes, j]
      e$mParB[4, j] = 0
    }
  }

  e$LI = matrix(nrow=2, ncol=e$nPara) # Likelihood interval
  colnames(e$LI) = e$pNames[1:e$nPara]
  rownames(e$LI) = c("Lower", "Upper")
  for (j in 1:e$nPara) {
    fx = function(x) {
      tPar = e$PE
      tPar[j] = x
      e$Obj(tPar) - e$r$value - e$fCut
    }
    if (e$mParB[3, j] == 1 & is.finite(e$mParB[1, j]) & is.finite(e$PE[j])) {
      rTemp = try(uniroot(fx, c(e$mParB[1, j], e$PE[j]))$root, silent=TRUE)
      if (!inherits(rTemp, "try-error")) { e$LI[1, j] = rTemp
      } else { e$LI[1, j] = NA_real_ }
    } else {
      e$LI[1, j] = NA_real_
    }
    if (e$mParB[4, j] == 1 & is.finite(e$mParB[2, j]) & is.finite(e$PE[j])) {
      rTemp = try(uniroot(fx, c(e$PE[j], e$mParB[2, j]))$root, silent=TRUE)
      if (!inherits(rTemp, "try-error")) { e$LI[2, j] = rTemp
      } else { e$LI[2, j] = NA_real_ }
    } else {
      e$LI[2, j] = NA_real_
    }
  }
##

  e$run = run.test(e$Residual)
  e$'-2LL' = e$nRec*log(2*pi) + e$r$value
  e$AIC = e$'-2LL' + 2*e$nPara
  e$AICc = e$AIC + 2*e$nPara*(e$nPara + 1)/(e$nRec - e$nPara - 1)
  e$BIC = e$'-2LL' + e$nPara*log(e$nRec)
  e$Elapsed = difftime(Sys.time(), t1)
  
  if (e$Error == "A") {
    Result = list(e$Est, e$LI, e$HouSkew, e$Cov, e$run, e$r$value, e$'-2LL', e$AIC, e$AICc, e$BIC, e$r$covergence, e$r$message, e$Pred, e$Residual)
    names(Result) = c("Est", "LI", "Skewness", "Cov", "run", "Objective Function Value", "-2LL", "AIC", "AICc", "BIC", "Convergence", "Message", "Prediction", "Residual")
  } else {
    Result = list(e$Est, e$LI, e$Cov, e$run, e$r$value, e$'-2LL', e$AIC, e$AICc, e$BIC, e$r$covergence, e$r$message, e$Pred, e$Residual)
    names(Result) = c("Est", "LI", "Cov", "run", "Objective Function Value", "-2LL", "AIC", "AICc", "BIC", "Convergence", "Message", "Prediction", "Residual")
  }
  Len0 = length(Result)
  Name0 = names(Result)

  if (toupper(Error) != "S") {
    Result[[Len0 + 1]] = e$Elapsed
    names(Result) = c(Name0, "Elapsed Time")
  } else {
    Scale = Sx(e$Est["PE", 1:e$nTheta])
    Result[[Len0 + 1]] = Scale
    Result[[Len0 + 2]] = e$Elapsed
    names(Result) = c(Name0, "Scale", "Elapsed Time")
  }
  return(Result)
}

