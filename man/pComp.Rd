\name{pComp}
\alias{pComp}
\title{Plot Compartment Model Diagram}
\description{
It plots the diagrom of a comparment model.
}
\usage{
pComp(dComp, dRate, Shape="rect", Bx=0.3, By=0.2, Cex=1.0, Lwd=3, Radius=0.3, 
      thIn=pi/2, thOut=pi/2, ...)
}
\arguments{
  \item{dComp}{data.frame for a compartment model. See the example.}
  \item{dRate}{data.frame for rate information. See the example.}
  \item{Shape}{Rectangle or cricle}
  \item{Bx}{half width of compartment box}
  \item{By}{half height of compartment box}
  \item{Cex}{character expansion}
  \item{Lwd}{line width}
  \item{Radius}{radius of compartment circle}
  \item{thIn}{Input angle in radian} 
  \item{thOut}{Output angle in radian}
   \item{...}{arguments to be passed to \code{plot} function}
}
\details{
  Flow direction is from the top to bottom.
}
\value{
  It plots.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
dA = data.frame(No = c(1, 2, 3, 4), Name=c("Gut Depot", "Skin Depot", "Central", "Peripheral"), 
                Level=c(1, 1, 2, 2), xPos=c(-0.5, 0.5, 0, 1))
dB = data.frame(From = c(1, 2, 3, 4, 3, 0, 0), To=c(3, 3, 4, 3, 5, 1, 2), 
                Name=c("KA", "KA2", "K12", "K21", "CL", "F1", "F2"))

#par(oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
pComp(dA, dB)
#par(oma=c(0, 0, 0, 0), mar=c(0, 0, 1, 0))
pComp(dA, dB, "circ", main="Compartmental Model Diagram")

pComp(dA, dB, "circ", main="Compartmental Model Diagram", thIn=pi/4, thOut=pi/4)
}
\keyword{Multi-compartment}
