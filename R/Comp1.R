Comp1 = function(Ke, Ka=0, DH) # DH: Dosing History
{
  NCOMP = 1

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

    dT = DH[i, "TIME"] - DH[i - 1, "TIME"]
    cR = DH[i - 1, "RATE2"]
    E  = exp(-Ke*dT)

    if (abs(Ka - Ke) > 1e-8) {
      Ea = exp(-Ka*dT)
      X[i, 1] = pX[1]*Ea
      X[i, 2] = pX[2]*E + cR*(1 - E)/Ke + pX[1]*Ka*(E - Ea)/(Ka - Ke)
    } else {
      X[i, 1] = pX[1]*E
      X[i, 2] = pX[2]*E + cR*(1 - E)/Ke + pX[1]*Ke*dT*E
    }
  }

  return(X)
}
