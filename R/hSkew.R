hSkew = function(rx)
{ # rx: result of nls
  bx = coef(rx)
  pName = names(bx)
  z = length(bx)
  ssq = sigma(rx)^2
  m = length(resid(rx))
  d1 = eval(rx$data)
  fm1 = formula(rx)
  tm1 = intersect(attr(terms(fm1), "term.labels"), colnames(d1))
  if (length(tm1) > 0) {
    for (i in 1:length(tm1)) assign(tm1[i], d1[, tm1[i]], envir=as.environment("Autoloads"))
  }
  f1 = deriv(fm1, pName, function.arg=pName, func=TRUE, hessian=TRUE)
  rf1 = do.call(f1, as.list(bx))
  return(Hougaard(attr(rf1, "gradient"), attr(rf1, "hessian"), ssq))
}
