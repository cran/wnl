pProf = function(Bag = e, Title = "", ...)
{
  if (Bag$nPara <= 2) {
    mfRow = 1
    mfCol = 2
  } else {
    mfRow = ceiling(sqrt(Bag$nPara))
    mfCol = ceiling(Bag$nPara/mfRow)
  }
  oPar = par(mfrow=c(mfRow, mfCol))
  Args = list(...)
  if (is.null(Args$ylab)) Args$ylab = "-2LL"
  if (is.null(Args$type)) Args$type = "l"

  for (j in 1:Bag$nPara) {
    x = Bag$mPar[, j]
    y = Bag$mOFV[, j]
    if (is.finite(min(x)) & is.finite(max(x)) & is.finite(min(y, na.rm=T)) & is.finite(max(y, na.rm=T))) {
      Args$x = x
      Args$y = y
      RdUdL = format((Bag$LI[2, j] - Bag$PE[j])/(Bag$PE[j] - Bag$LI[1, j]), digits=3)
      Args$xlab = paste0(Bag$pNames[j], " = ", format(Bag$PE[j], digits=2), ", dU/dL = ", RdUdL)
      do.call(plot, Args)
      abline(h = Bag$'-2LL' + Bag$fCut, lty=2)
      abline(v = Bag$PE[j], lty=3)
      text(Bag$LI[, j], Bag$'-2LL', labels = format(Bag$LI[, j], digits=2))
    }
  }
  if (trimws(Title) != "") title(Title, outer=TRUE)
  par(oPar)
}

