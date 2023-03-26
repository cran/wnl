EnvObj = function(envir = e)
{
  Name0 = ls(envir)
  nObj = length(Name0)
  Result = list()
  for (i in 1:nObj) Result[[i]] = get(Name0[i], envir=envir)
  names(Result) = Name0
  return(Result)
}
