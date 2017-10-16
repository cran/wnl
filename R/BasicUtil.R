#BasicUtil.R

#Note = function() file.show(system.file("NOTE.txt", package="NonCompart"))

UT  = function(x) toupper(gsub("^\\s+|\\s+$", "", x))

s2o = function(vPara) # scaled to original parameters
{
  b0 = exp(vPara - e$alpha)
  return(b0/(b0 + 1)*(e$UB - e$LB) + e$LB)
}

prun = function(m, n, r)
{
# INPUT
# m : count of fewer species (minimum value = 0)
# n : count of more frequent species (minimum value = 1)
# r : count of run (minimum value = 1)
# RETURNS probability of run count to be less than or equal to r with m and n
#         P(Run count <= r | m, n)
# Reference: 
#   Reviewed Work: Run Test for Randomness by Dorothy Kerr
#   Review by: J. W. W.
#   Mathematics of Computation
#   Vol. 22, No. 102 (Apr., 1968), p. 468
#   Published by: American Mathematical Society
#   DOI: 10.2307/2004702
#   Stable URL: http://www.jstor.org/stable/2004702

  if (m==0 & r==1) {
    return(1)
  } else if (m > n | m < 1 | n < 1 | r < 2 | (r > min(m + n, 2 * m + 1))) {
    return(0);
  }

  sumfu = 0
  for (u in 2:r) {
    if (u %% 2 == 0) {
      k = u/2
      fu = 2*choose(m-1, k-1)*choose(n-1, k-1)
    } else {
      k = (u + 1)/2
      fu = choose(m-1, k-1)*choose(n-1, k-2) + choose(m-1, k-2)*choose(n-1, k-1)
    }
    sumfu = sumfu + fu
  }
  return(sumfu / choose(m + n, m))
}

run.test = function(Residuals)
{
  Resid = Residuals[Residuals != 0] # Zeros are removed.
  nResid = length(Resid)
  r = Resid > 0
  m = sum(r)
  if (nResid > 1) {
    j = 2:nResid
    run = sum(abs(r[j] - r[j - 1])) + 1
  } else {
    run = 1
  }
  m = min(m, nResid - m)
  n = max(m, nResid - m)
  if (run > 1) {
    p = prun(m, n, run)
    if (p > 0.5) {
      p = 1 - prun(m, n, run - 1)
    }
  } else {
    p = 0.5^(nResid - 1)
  }
  
  Result = list(m, n, run, p)
  names(Result) = c("m", "n", "run", "p.value")
  return(Result)
}
