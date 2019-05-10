SolComp3 = function(K10, K12, K21, K13, K31)
{
# Get Lambdas
  A1 = K10 + K12 + K13 + K21 + K31
  A2 = K10*K21 + K10*K31 + K12*K31 + K13*K21 + K21*K31
  A3 = K21*K31*K10

  Q = (A1*A1 - 3*A2)/9
  RQ = 2*sqrt(Q)
  R = (2*A1*A1*A1 - 9*A1*A2 + 27*A3)/54
  M = Q*Q*Q - R*R
  if (M < 0) stop("Error: Not real roots.")

  Th = acos(8*R/(RQ*RQ*RQ))
  L1 = RQ*cos(Th/3) + A1/3
  L2 = RQ*cos((Th + 2*pi)/3) + A1/3
  L3 = RQ*cos((Th + 4*pi)/3) + A1/3
  L = c(L1, L2, L3)

# Get Coefficients
  D1 = (L[2] - L[1])*(L[3] - L[1])
  D2 = (L[1] - L[2])*(L[3] - L[2])
  D3 = (L[1] - L[3])*(L[2] - L[3])
  Div = c(D1, D2, D3)        # Divisor, denominator
  if (prod(Div) == 0) stop("Roots should be distinct real values.")

  Co = array(dim=c(3, 3, 3)) # Dividend, numerator, NCOMP=3
  Co[1,,1] = (K21 - L)*(K31 - L)
  Co[1,,2] = K21*(K31 - L)
  Co[1,,3] = K31*(K21 - L)
  Co[2,,1] = K12*(K31 - L)
  Co[2,,2] = (K10 + K12 + K13 - L)*(K31 - L) - K31*K13
  Co[2,,3] = K12*K31
  Co[3,,1] = K13*(K21 - L)
  Co[3,,2] = K21*K13
  Co[3,,3] = (K10 + K12 + K13 - L)*(K21 - L) - K21*K12
  for (i in 1:3) Co[,i,] = Co[,i,]/Div[i]

  return(list(L=L, Co=Co))
}
