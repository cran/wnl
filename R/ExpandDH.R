ExpandDH = function(DH, Fo = 1)
{
  DH[is.na(DH)] = 0
  DH[DH$AMT > 0 & DH$CMT == 1, "AMT"] =  Fo*DH[DH$AMT > 0 & DH$CMT == 1, "AMT"]

  DH[,"BOLUS"] = 0
  DH[DH$AMT > 0 & DH$RATE == 0, "BOLUS"] = DH[DH$AMT > 0 & DH$RATE == 0, "AMT"] 
  
  DH[,"RATE2"] = 0
  BlankRow = DH[1,]
  BlankRow[1,] = rep(0, ncol(DH))
  RateDat = DH[DH$RATE > 0,]
  nRate = nrow(RateDat)
  for (i in 1:nRate) {
    cRATE = RateDat[i, "RATE"]
    cCMT = RateDat[i, "CMT"]
    cTIME = RateDat[i, "TIME"]
    cDUR = RateDat[i, "AMT"] / RateDat[i,"RATE"]
    cEND = cTIME + cDUR
    if (!(cEND %in% DH[,"TIME"])) {

      BlankRow[1, "TIME"] = cEND
      DH = rbind(DH, BlankRow)
      DH = DH[order(DH$TIME),]
    }
    DH[DH$TIME >= cTIME & DH$TIME < (cTIME + cDUR), "RATE2"] = DH[DH$TIME >= cTIME & DH$TIME < (cTIME + cDUR), "RATE2"] + cRATE
    DH[DH$TIME >= cTIME & DH$TIME < (cTIME + cDUR), "CMT"] = cCMT
  }
  return(DH)
}
