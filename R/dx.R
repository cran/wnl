dx = function(r)
{
  nRec = length(r$Prediction)
  EstNames = colnames(r$Est)
  wt = rep(0, nRec)
  if ("AddErrVar" %in% EstNames) wt = wt + r$Est["PE","AddErrVar"]
  if ("PoisErrVar" %in% EstNames) wt = wt + r$Prediction*r$Est["PE","PoisErrVar"]
  if ("PropErrVar" %in% EstNames) wt = wt + r$Prediction^2*r$Est["PE","PropErrVar"]
  if ("ScaleErrVar" %in% EstNames) wt = r$Scale*r$Est["PE","ScaleErrVar"]

  if (attr(dev.cur(), "names") == "null device") {
    dev.new(width=12, height=6)
    DefPar = par(mfrow=c(1,2))
  }
  plot(r$Prediction, r$Prediction + r$Residual, xlab="Prediction", ylab="Observation", pch=16)
  abline(a=0, b=1, lty=3)

  if (min(wt) > 0) {
    plot(r$Prediction, r$Residual/sqrt(wt), xlab="Prediction", ylab="Normalized Residual", pch=16)
    abline(h=(-10:10), lty=3)
    abline(h=0)
  } else {
    plot(r$Prediction, r$Residual, xlab="Prediction", ylab="Residual", pch=16)
    abline(h=0)
  }
  if ("DefPar" %in% ls()) par(DefPar)
}
