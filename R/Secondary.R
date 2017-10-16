Secondary = function(Formula, PE, COV)
{
  nArg = length(PE)
  gr0 = deriv(Formula, names(PE), function.arg=names(PE), func=TRUE)
  gr1 = do.call("gr0", as.list(PE))
  PE2 = gr1[1]
  gr2 = attr(gr1,"gradient") 
  SE2 = sqrt(gr2 %*% COV %*% t(gr2))
  CV2 = SE2 / PE2 *100  
  Result = c(PE2, SE2, CV2)
  names(Result) = c("PE", "SE", "RSE")
  return(Result)
}
