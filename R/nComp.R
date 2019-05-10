nComp = function(Sol, Ka, DH)
{
  L = Sol$L
  Co = Sol$Co
  NCOMP = length(L)
  if (NCOMP < 2) stop("Compartment count should be at least 2.")
  if (NCOMP != dim(Co)[1]) stop("Lengths of lambda and coefficients mismatch.")

  NTIME = nrow(DH)
  if (NTIME < 2) stop("Dosing history table should have at least two rows.")

  X = matrix(rep(0, (NCOMP + 1)*NTIME), ncol=(NCOMP + 1), nrow=NTIME)
  for (i in 2:NTIME) {
    pX = X[i - 1,]

    for (j in 1:(NCOMP + 1)) {
      if (DH[i - 1, "CMT"] == j & DH[i - 1, "BOLUS"] > 0) {
        pX[j] = pX[j] + DH[i - 1, "BOLUS"]
      }
    }

    dT = DH[i, "TIME"] - DH[i - 1, "TIME"]               # delta T
    cR = DH[i - 1, "RATE2"]                              # Infusion Rate
    E = exp(-L*dT)                                       # Exponentials

    Xo = rep(0, NCOMP)
    for (j in 1:NCOMP) Xo = Xo + pX[1+j] * Co[,,j] %*% E # Bolus

    if (cR > 0) Xo = Xo + ((cR*Co[,,1]) %*% ((1 - E)/L)) # Infusion

    Ea = exp(-Ka*dT)                                     # Oral
    if (pX[1] > 0) Xo = Xo + Ka*pX[1]*(Co[,,1] %*% ((E - Ea)/(Ka - L)))

    X[i,] = c(pX[1]*Ea, Xo)
  }
  return(X)
}
