SolComp2 = function(K10, K12, K21)
{
# Get Lambdas
  Sum = K10 + K12 + K21
  Disc = sqrt(Sum*Sum - 4*K10*K21)
  L1 = (Sum + Disc)/2
  L2 = (Sum - Disc)/2
  L = c(L1, L2)

# Get Coefficients
  Div = c(L2 - L1, L1 - L2) # Divisor, denominator
  if (prod(Div) == 0) stop("Roots should be distinct real values.")

  Co = array(dim=c(2, 2, 2)) # Dividend, numerator, NCOMP=2
  Co[1, , 1] = K21 - L
  Co[1, , 2] = K21
  Co[2, , 1] = K12
  Co[2, , 2] = K10 + K12 - L
  for (i in 1:2) Co[, i, ] = Co[, i, ]/Div[i]

  return(list(L=L, Co=Co))
}
