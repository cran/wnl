#BasicUtil.R

#Note = function() file.show(system.file("NOTE.txt", package="NonCompart"))

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

nGradient = function(func, x)
{
  nVar = length(x)
  nRec = length(func(x))
  x1 = vector(length = nVar)
  x2 = vector(length = nVar)
  mga = matrix(nrow=nRec, ncol = 4)
  mgr = matrix(nrow=nRec, ncol = nVar)

  for (i in 2:nVar) x1[i] = x2[i] = x[i]

  for (i in 1:nVar) {
    axi = abs(x[i])
    if (axi <= 1) { hi = 1e-04 
    } else { hi = 1e-04*axi }
    for (k in 1:4) {
      x1[i] = x[i] - hi
      x2[i] = x[i] + hi
      mga[, k] = (func(x2) - func(x1))/(2*hi)
      hi = hi/2
    }
    mga[, 1] = (mga[, 2]*4 - mga[, 1])/3
    mga[, 2] = (mga[, 3]*4 - mga[, 2])/3
    mga[, 3] = (mga[, 4]*4 - mga[, 3])/3
    mga[, 1] = (mga[, 2]*16 - mga[, 1])/15
    mga[, 2] = (mga[, 3]*16 - mga[, 2])/15
    mgr[, i] = (mga[, 2]*64 - mga[, 1])/63
    x1[i] = x2[i] = x[i]
  }

  return(mgr)
}

nHessian = function(fx, x)
{
  nVar  = length(x)
  h0 = vector(length=nVar)
  x1 = vector(length=nVar)
  x2 = vector(length=nVar)

  f0 = fx(x)
  nRec = length(f0)

  ha = matrix(nrow=nRec, ncol=4) # Hessian Approximation
  H = rep(0, nRec*nVar*nVar)     # Hessian Matrix
  dim(H) = c(nRec, nVar, nVar)

  for (i in 1:nVar) {
    x1[i] = x2[i] = x[i]
    axi   = abs(x[i])
    if (axi < 1) { h0[i] = 1e-4
    } else       { h0[i] = 1e-4*axi }
  }

  for (i in 1:nVar) {
    for (j in i:1) {
      hi = h0[i]
      if (i == j) {
        for (k in 1:4) {
          x1[i] = x[i] - hi
          x2[i] = x[i] + hi
          ha[, k] = (fx(x1) - 2*f0 + fx(x2))/(hi*hi)
          hi = hi/2
        }
      } else {
        hj = h0[j]
        for (k in 1:4) {
          x1[i] = x[i] - hi
          x1[j] = x[j] - hj
          x2[i] = x[i] + hi
          x2[j] = x[j] + hj
          ha[, k] = (fx(x1) - 2*f0 + fx(x2) - H[, i, i]*hi*hi - H[, j, j]*hj*hj)/(2*hi*hj)
          hi = hi / 2
          hj = hj / 2
        }
      }
      w = 4
      for (m in 1:2) {
        for (k in 1:(4 - m)) ha[, k] = (ha[, k + 1]*w - ha[, k])/(w - 1)
        w = w*4
      }
      H[, i, j] = (ha[, 2]*64 - ha[, 1]) / 63
      if (i != j) H[, j, i] = H[, i, j] 
      x1[j] = x2[j] = x[j]
    }
    x1[i] = x2[i] = x[i]
  }

  return(H)
}

g2inv = function(A, eps=1e-8)
{
  idx = abs(diag(A)) > eps
  p = sum(idx, na.rm=T)
  M = matrix(0, nrow=NCOL(A), ncol=NROW(A))
  if (p == 0) {attr(M, "rank") = 0 ; return(M) }
  B = A[idx, idx, drop=F]
  r = 0
  for (k in 1:p) {
    d = B[k, k]
    if (abs(d) < eps) { B[k, ] = 0 ; B[, k] = 0 ; next }
    B[k, ] = B[k, ]/d
    r = r + 1
    for (i in 1:p) {
      if (i != k) {
        c0 = B[i, k]
        B[i, ] = B[i, ] - c0*B[k, ]
        B[i, k] = -c0/d
      }
    }
    B[k, k] = 1/d
  }
  M[1:r, 1:r] = B[1:r, 1:r]
  attr(M, "rank") = r
  return(M)
}

Hougaard = function(J, H, ssq)
{# J : graident, H: hessian, ssq: sigma square  
  z = NCOL(J)
  m = NROW(J)
  if (z*m == 0) stop("No graident information!")

  L = g2inv(crossprod(J))
  if (attr(L, "rank") < ncol(L)) warning("Crossproduct of gradient is singular!")

  W = rep(0, z^3)
  dim(W) = rep(z, 3)
  for (k in 1:z) {
    for (p in 1:z) {
      for(j in 1:z) {
        for (i in 1:m) W[k, p, j] = W[k, p, j] + J[i, k]*H[i, p, j]
      }
    }
  }

  TM = rep(NA, z)
  for (i in 1:z) {
    tRes = 0
    for (j in 1:z) {
      for (k in 1:z) {
        for (p in 1:z) {
          tRes = tRes + L[i, j]*L[i, k]*L[i, p]*(W[j, k, p] + W[k, j, p] + W[p, k, j])
        }
      }
    }
    TM[i] = -ssq^2*tRes
  }
  
  SK = rep(NA, z)
  for (i in 1:z) SK[i] = TM[i]/(ssq*L[i, i])^1.5
  names(SK) = colnames(J)
  return(SK)
}
