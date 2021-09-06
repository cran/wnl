pComp = function(dComp, dRate, Shape="rect", Col=NA, Bx=0.3, By=0.2, Cex=1.0, Lwd=3, Radius=0.3, thIn=pi/2, thOut=pi/2, ...)
{
  if (!(is.data.frame(dComp) & is.data.frame(dRate))) stop("Two input data.frames are needed!")

  Shape = substr(toupper(trimws(Shape)), 1, 4)
  if (Shape == "CIRC") {
    Bx = Radius
    By = Radius
  }
  th0 = atan2(By, Bx)
  InOutL = 2*By

  d1 = dComp
  nComp = NROW(d1)
  if (nComp == 0) stop("There should be at least one compartment!")

  xmax = max(d1$xPos) + Bx
  xmin = min(d1$xPos) - Bx
  xRange = xmax - xmin
  tx = xRange/100 # tiny delta x

  d1$yPos = max(d1$Level) - d1$Level
  ymax = max(d1$yPos) + 2.5*By
  ymin = min(d1$yPos) - 2.5*By
  yRange = ymax - ymin
  ty = yRange/100 # tiny delta y

## Draw boxes & texts
  plot(0, 0, type="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", bty="n", axes=F, ...)
  if (Shape == "CIRC") {
    for (i in 1:nComp) {
      x0 = d1[i, "xPos"]
      y0 = d1[i, "yPos"]
      th = c(seq(0, 2*pi, length.out=200), 0)
      x = x0 + Radius*cos(th)
      y = y0 + Radius*sin(th)
      polygon(x, y, lwd=Lwd, col=Col)
    }
  } else {
    for (i in 1:nComp) {
      x0 = d1[i, "xPos"]
      y0 = d1[i, "yPos"]
      x = c(x0 - Bx, x0 + Bx, x0 + Bx, x0 - Bx, x0 - Bx)
      y = c(y0 + By, y0 + By, y0 - By, y0 - By, y0 + By)
      polygon(x, y, lwd=Lwd, col=Col)
    }
  }
  text(d1$xPos, d1$yPos, d1$Name, cex=Cex)

## Draw Rates
  d2 = dRate
  nRate = NROW(d2)
  if (nRate == 0) return()
  if (nrow(unique(dRate[,c("From", "To")])) != nrow(dRate)) stop("Rate data.frame should have unique rates!")

  d2$Both = FALSE
  for (i in 1:nrow(d2)) {
    cFrom = d2[i, "From"]
    cTo = d2[i, "To"]
    if (nrow(d2[d2$From == cTo & d2$To == cFrom,]) > 0) d2[i, "Both"] = TRUE
  }

  d2$s0 = NA #side
  d2$s1 = NA
  d2$q01 = NA # quadrant
  d2$q10 = NA
  d2$x0 = NA
  d2$y0 = NA
  d2$x1 = NA
  d2$y1 = NA

  for (i in 1:nRate) {
    n1 = d2[i, "From"] # From compartment number
    n2 = d2[i, "To"]   # To compartment number

  ## determine ceter of From point
    if (n1 == 0) { # input side
      cx0 = d1[d1$No == n2, "xPos"] - InOutL*cos(thIn)
      cy0 = d1[d1$No == n2, "yPos"] + InOutL*sin(thIn) + 1.5*By
    } else {      # input side
      cx0 = d1[d1$No == n1, "xPos"]
      cy0 = d1[d1$No == n1, "yPos"]
    }

  ## determine ceter of To point
    if (n2 > nComp) { # normal compartment
      cx1 = d1[d1$No == n1, "xPos"] + InOutL*cos(thOut)
      cy1 = d1[d1$No == n1, "yPos"] - InOutL*sin(thOut) - 1.5*By
    } else {           # output side
      cx1 = d1[d1$No == n2, "xPos"]
      cy1 = d1[d1$No == n2, "yPos"]
    }

    mxa = mean(c(cx0, cx1))
    mya = mean(c(cy0, cy1))

    th1 = atan2(cy1 - cy0, cx1 - cx0)
    th1 = ifelse(th1 >= 0, th1, th1 + 2*pi)
    q01 = floor(th1/ (pi/2)) + 1  # quadrant
    if (th1 >= th0 & th1 <= (pi - th0)) {
      s0 = 3
      x0 = cx0
      y0 = cy0 + By
    } else if (th1 > (pi - th0) & th1 < pi + th0) {
      s0 = 2
      x0 = cx0 - Bx
      y0 = cy0
    } else if (th1 >= pi + th0 & th1 <= 2*pi - th0) {
      s0 = 1
      x0 = cx0
      y0 = cy0 - By
    } else {
      s0 = 4
      x0 = cx0 + Bx
      y0 = cy0
    }

    th2 = atan2(cy0 - cy1, cx0 - cx1)
    th2 = ifelse(th2 >= 0, th2, th2 + 2*pi)
    q10 = floor(th2/ (pi/2)) + 1
    if (th2 >= th0 & th2 <= (pi - th0)) {
      s1 = 3
      x1 = cx1
      y1 = cy1 + By
    } else if (th2 > (pi - th0) & th2 < pi + th0) {
      s1 = 2
      x1 = cx1 - Bx
      y1 = cy1
    } else if (th2 >= pi + th0 & th2 <= 2*pi - th0) {
      s1 = 1
      x1 = cx1
      y1 = cy1 - By
    } else {
      s1 = 4
      x1 = cx1 + Bx
      y1 = cy1
    }

    if (d2[i, "Both"]) {
      if (y0 == y1) {
        y0 = y0 + ty*sign(x1 - x0)
        y1 = y1 + ty*sign(x1 - x0)
      } else if (x0 == x1) {
        x0 = x0 + 0.5*tx*sign(y1 - y0)
        x1 = x1 + 0.5*tx*sign(y1 - y0)
      } else if ((s0 == 3 & q01 == 1) | (s0 == 1 & q01 == 4)) {
        x0 = cx0 + 3*tx
        x1 = cx1 - tx
      } else if ((s0 == 4 & q01 == 1) | (s0 == 2 & q01 == 2) ) {
        y0 = cy0 + ty
        y1 = cy1 - 3*ty
      } else if ((s0 == 3 & q01 == 2) | (s0 == 1 & q01 == 3)) {
        x0 = cx0 - 3*tx
        x1 = cx1 + tx
      } else if ((s0 == 4 & q01 == 4) | (s0 == 2 & q01 == 3)) {
        y0 = cy0 - ty
        y1 = cy1 + 3*ty
      }
    }
    arrows(x0, y0, x1, y1, length=0.12, lwd=Lwd, angle=20)

#    d2[i, c("s0", "s1", "q01", "q10", "x0", "y0", "x1", "y1")] = c(s0, s1, q01, q10, x0, y0, x1, y1)

    mx = mean(c(x0, x1))
    my = mean(c(y0, y1))
    if (y0 == y1) { # same level
      text(mx, my + sign(x1 - x0)*2.5*ty, d2[i, "Name"], cex=Cex)
    } else if (y0 > y1) { # from high to low
      text(mx - tx, my, d2[i, "Name"], cex=Cex, adj=1)
    } else { # from low to high
      text(mx + tx, my, d2[i, "Name"], cex=Cex, adj=0)
    }
  }
}
