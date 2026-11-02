## Plot+Colours----

myHeightEPS <- 6.5
myWidthEPS <- 10

myWidth <- 1070
myHeight <- 485

pdfWidth <- 14
pdfHeight <- 7

myCex <- 2
myLwd <- 5
myCexLab <- 1.5
myCexAxis <- 2.25
cexFactor <- 1.3

## Colours -----
colsBer <- c(committed    = "#BD8E17FF",
             universal    = "#E31A1CE6")
labsBer <- c(latex2exp::TeX(r"($b_{N,\alpha}=)"),
             latex2exp::TeX(r"($1/\alpha=20)"))

freqColours <- c("#E31A1CE6", "#FB9A9980")
bayesColours <- c("#B15928E6", "#FFFF9980")
eColours <- c("#1F78B4E6", "#A6CEE380")


committedColour <- "#BD8E17FF"
afterHorizonColour <- "darkorchid3"
universalColour <- "#E31A1CE6"

library(Matrix)
library(Rmpfr)

# Computes the exact crossing probability under an arbitrary H0 (p0)
computeVilleTerms <- function(nHorizon, bBoundary, pAlt, p0 = 0.5, precisionBits = 256) {
  uVal <- pAlt / p0
  dVal <- (1 - pAlt) / (1 - p0)

  totalNodes <- (nHorizon + 1) * (nHorizon + 2) / 2
  nodeVal    <- numeric(totalNodes)
  nodeN      <- integer(totalNodes)
  stateProb  <- numeric(totalNodes)

  getIdx <- function(n, k) {
    if (k < 0 || k > n) return(NULL)
    return((n * (n + 1) / 2) + k + 1)
  }

  for (n in 0:nHorizon) {
    for (k in 0:n) {
      idx <- getIdx(n, k)
      nodeN[idx]   <- n
      nodeVal[idx] <- (uVal^k) * (dVal^(n - k))
    }
  }

  stateProb[getIdx(0, 0)] <- 1.0

  for (n in 0:(nHorizon - 1)) {
    for (k in 0:n) {
      fromIdx <- getIdx(n, k)
      pCurrent <- stateProb[fromIdx]

      if (pCurrent == 0.0) next
      if (nodeVal[fromIdx] >= bBoundary) next

      idxUp   <- getIdx(n + 1, k + 1)
      idxDown <- getIdx(n + 1, k)

      stateProb[idxUp]   <- stateProb[idxUp] + (pCurrent * p0)
      stateProb[idxDown] <- stateProb[idxDown] + (pCurrent * (1 - p0))
    }
  }

  bMpfr    <- mpfr(bBoundary, precisionBits)
  oneMpfr  <- mpfr(1.0, precisionBits)
  vBn      <- mpfr(0.0, precisionBits)
  rBn      <- mpfr(0.0, precisionBits)

  # Clean, independent loop for overshoot tracking
  for (idx in 1:totalNodes) {
    if (nodeVal[idx] >= bBoundary) {
      if (stateProb[idx] > 0) {
        currentMass <- mpfr(stateProb[idx], precisionBits)
        currentNodeValue <- mpfr(nodeVal[idx], precisionBits)
        vBn <- vBn + (currentMass * (currentNodeValue - bMpfr))
      }
    }
  }

  # Clean, independent loop for residual tracking
  finalNodesIdx <- which(nodeN == nHorizon & nodeVal < bBoundary)
  for (idx in finalNodesIdx) {
    if (stateProb[idx] > 0) {
      currentMass <- mpfr(stateProb[idx], precisionBits)
      currentNodeValue <- mpfr(nodeVal[idx], precisionBits)
      rBn <- rBn + (currentMass * currentNodeValue)
    }
  }

  pCrossing <- (oneMpfr - rBn - vBn) / bMpfr

  result <- list("pCrossing"=pCrossing, "rBn"=rBn,
                 "vBn"=vBn, "b"=bMpfr, "power"=bMpfr*pCrossing+vBn)
  return(result)
}

# Root-finding algorithm utilizing the Bisection Method
findBoundaryForTargetProb <- function(nHorizon, targetProb, pAlt, p0 = 0.5, tol = 1e-6) {
  lowerB <- 1.0
  upperB <- 1.0 / targetProb

  # Clean structural safety expansion loop without any inline assignment leaks
  pAtUpper <- computeVilleTerms(nHorizon, upperB, pAlt, p0)$pCrossing
  while (pAtUpper > targetProb) {
    upperB <- upperB * 2.0
    pAtUpper <- computeVilleTerms(nHorizon, upperB, pAlt, p0)$pCrossing
  }

  while ((upperB - lowerB) > tol) {
    midB <- (lowerB + upperB) / 2.0
    currentP <- computeVilleTerms(nHorizon, midB, pAlt, p0)$pCrossing

    if (currentP > targetProb) {
      lowerB <- midB
    } else {
      upperB <- midB
    }
  }

  return((lowerB + upperB) / 2.0)
}

# 0.9 vs 0.5 bound ----
# Run a clean verification execution trace
exactBoundaryValue <- findBoundaryForTargetProb(
  nHorizon = 50, targetProb = 0.05, pAlt = 0.9, p0 = 0.5)


tempRes <- computeVilleTerms(50, 16.65838, pAlt=0.9)

round(as.numeric(tempRes$pCrossing), 4)
round(as.numeric(tempRes$b*tempRes$pCrossing), 4)
round(as.numeric(tempRes$vBn), 4)
round(as.numeric(tempRes$rBn), 4)
round(as.numeric(tempRes$power), 4)

tempRes <- computeVilleTerms(50, 20, pAlt=0.9)
round(as.numeric(tempRes$pCrossing), 4)
round(as.numeric(tempRes$b*tempRes$pCrossing), 4)
round(as.numeric(tempRes$vBn), 4)
round(as.numeric(tempRes$rBn), 4)
round(as.numeric(tempRes$power), 4)


# 0.6 vs 0.5 bound----
exactBoundaryValue06 <- findBoundaryForTargetProb(
  nHorizon = 150, targetProb = 0.05,
  pAlt = 0.6, p0 = 0.5)

tempRes <- computeVilleTerms(150, 13.15491, pAlt=0.6)
as.numeric(tempRes$pCrossing)

computePowerCommitted06 <- Vectorize(function(n){
  as.numeric(computeVilleTerms(n, 13.15491, 0.6, 0.5)$power)})

computePowerSavi06 <- Vectorize(function(n){
  as.numeric(computeVilleTerms(n, 20, 0.6, 0.5)$power)})

powerCommitted06 <- computePowerCommitted06(1:300)
powerSavi06 <- computePowerSavi06(1:300)


# myName <- "exampleBinomialPower06"
# pdf(paste0(myName, ".pdf"), width=pdfWidth, height=pdfHeight)
# graphics::par(cex.main=1.5, mar=c(5, 8.3, 2, 0)+0.1,
#               mgp=c(3.5, 1, 0), cex.lab=2.5,
#               font.lab=2, cex.axis=1.3, bty="n", las=1)
#
# plot(1:300, powerCommitted06, xlab="",
#      ylab="", col=NULL, xaxt="n", yaxt="n",
#      bty="n")
#
# lines(c(150, 150), c(0, 1),
#       lwd=myLwd, col="lightgrey", lty=2)
#
#
# lines(1:150, powerCommitted06[1:150], lwd=myLwd,
#       col=committedColour)
# lines(151:300, powerCommitted06[151:300], lwd=myLwd,
#       col=afterHorizonColour, lty=2)
# lines(1:300, powerSavi06, lwd=myLwd, col=universalColour)
# # abline(v=150)
#
#
#
# axis(1, cex.axis=cexFactor*myCex)
# axis(2, at=0:5*2/10, 0:5*20, cex.axis=cexFactor*myCex)
#
#
# mtext("Time horizon N", side = 1,
#       line = 2.5, las = 1, cex = cexFactor*myCex, padj=0.4)
#
# mtext(latex2exp::TeX(r"($Q(T_{b} \leq N)$ in %)"), side = 2, line = 4, las = 0,
#       cex = cexFactor*myCex, adj=0.5, padj=-0.4)
# dev.off()

N <- 150

# Power plot ----
pdf("exampleBinomialPower06.pdf", width = 14, height = 7)
par(mai = c(1.2, 1.5, 0.2, 1.7), cex.axis = 4/3, las = 1)

plot(1:300, powerCommitted06, xlab="",
     ylab="", col=NULL, xaxt="n", yaxt="n",
     bty="n", ylim=c(0, 1))

lines(c(N, N), c(0, 1),
      lwd=myLwd, col="lightgrey", lty=2)

lines(1:N, powerCommitted06[1:N], lwd=myLwd,
      col=committedColour)
lines((N+1):300, powerCommitted06[(N+1):300], lwd=myLwd,
      col=afterHorizonColour, lty=2)
lines(1:300, powerSavi06, lwd=myLwd, col=universalColour)
# abline(v=150)

axis(1, cex.axis=cexFactor*myCex)
axis(2, at=0:5*2/10, 0:5*20, cex.axis=cexFactor*myCex)


mtext("Time horizon N", side = 1,
      line = 2.5, las = 1, cex = cexFactor*myCex, padj=0.4)

mtext(latex2exp::TeX(r"($Q(T_{b} \leq N)$ in %)"), side = 2, line = 4, las = 0,
      cex = cexFactor*myCex, adj=0.5, padj=-0.4)

par(xpd = NA)


for (k in 1:2) {
  yCor <- if (k==1) powerCommitted06[300]+0.02 else if (k==2) powerSavi06[300]-0.02

  someLabel <- if (k==1) latex2exp::TeX(r"($b_{N,\alpha}=13.15)") else latex2exp::TeX(r"($1/\alpha=20)")
  someColour <- if (k==1) committedColour else universalColour
  text(300, yCor,
       someLabel, pos = 4, offset = 0.75, col = someColour,
       cex = 0.85*myCex)
}


text(N, par("usr")[4] + 0.0001, paste("N =", N),
     col = "gray40", cex = myCex)
# text(11, alpha + 0.0025, bquote(alpha == .(alpha)), adj = c(0, 0.5),
#      col = "gray40", cex = myCex)
dev.off()



computeFprCommitted06 <- Vectorize(function(n){
  as.numeric(computeVilleTerms(n, 13.15491, 0.6, 0.5)$pCrossing)})

computeFprSavi06 <- Vectorize(function(n){
  as.numeric(computeVilleTerms(n, 20, 0.6, 0.5)$pCrossing)})

fprCommitted06 <- computeFprCommitted06(1:300)
fprSavi06 <- computeFprSavi06(1:300)

# fpr plot ----
myName <- "exampleBinomialFpr06"
pdf(paste0(myName, ".pdf"), width = 14, height = 7)
par(mai = c(1.2, 1.5, 0.2, 1.7), cex.axis = 4/3, las = 1)

plot(1:300, fprCommitted06, xlab="",
     ylab="", col=NULL, xaxt="n", yaxt="n",
     bty="n")

lines(c(N, N), c(0, 1),
      lwd=myLwd, col="lightgrey", lty=2)

lines(c(1, 300), c(0.05, 0.05),
      lwd=myLwd, col="grey", lty=2)

lines(1:N, fprCommitted06[1:N], lwd=myLwd,
      col=committedColour)
lines((N+1):300, fprCommitted06[(N+1):300], lwd=myLwd,
      col=afterHorizonColour, lty=2)
lines(1:300, fprSavi06, lwd=myLwd, col=universalColour)
# abline(v=150)

axis(1, cex.axis=cexFactor*myCex)
axis(2, at=0:7/10, 0:7, cex.axis=cexFactor*myCex)


mtext("Time horizon N", side = 1,
      line = 2.5, las = 1, cex = cexFactor*myCex, padj=0.4)

mtext(latex2exp::TeX(r"($P(T_{b} \leq N)$ in %)"), side = 2, line = 4, las = 0,
      cex = cexFactor*myCex, adj=0.5, padj=-0.4)

par(xpd = NA)


for (k in 1:2) {
  yCor <- if (k==1) fprCommitted06[300] else if (k==2) fprSavi06[300]

  someLabel <- if (k==1) latex2exp::TeX(r"($b_{N,\alpha}=13.15)") else latex2exp::TeX(r"($1/\alpha=20)")
  someColour <- if (k==1) committedColour else universalColour
  text(300, yCor,
       someLabel, pos = 4, offset = 0.75, col = someColour,
       cex = 0.85*myCex)
}


text(N, par("usr")[4] + 0.0001, paste("N =", N),
     col = "gray40", cex = myCex)
text(1, alpha + 0.0025, bquote(alpha == .(alpha)), adj = c(0, 0.5),
     col = "gray40", cex = myCex)
dev.off()

# # fpr plot ----
# myName <- "exampleBinomialFpr06"
# pdf(paste0(myName, ".pdf"), width=pdfWidth, height=pdfHeight)
# graphics::par(cex.main=1.5, mar=c(5, 8.3, 2, 0)+0.1,
#               mgp=c(3.5, 1, 0), cex.lab=2.5,
#               font.lab=2, cex.axis=1.3, bty="n", las=1)
#
# plot(1:300, fprCommitted06, xlab="",
#      ylab="", col=NULL, xaxt="n", yaxt="n",
#      bty="n")
#
# lines(c(150, 150), c(0, 1),
#       lwd=myLwd, col="lightgrey", lty=2)
#
# lines(c(1, 300), c(0.05, 0.05),
#       lwd=myLwd, col="lightgrey", lty=2)
#
# lines(1:150, fprCommitted06[1:150], lwd=myLwd,
#       col=committedColour)
# lines(151:300, fprCommitted06[151:300], lwd=myLwd,
#       col=afterHorizonColour, lty=2)
# lines(1:300, fprSavi06, lwd=myLwd, col=universalColour)
# # abline(v=150)
#
#
#
# axis(1, cex.axis=cexFactor*myCex)
# axis(2, 0:7/100, 0:7, cex.axis=cexFactor*myCex)
# mtext("Time horizon N", side = 1,
#       line = 2.5, las = 1, cex = cexFactor*myCex, padj=0.4)
# mtext(latex2exp::TeX(r"($P(T_{b} \leq N)$ in %)"), side = 2, line = 4, las = 0,
#       cex = cexFactor*myCex, adj=0.5, padj=-0.4)
# dev.off()



tempRes <- computeVilleTerms(50, 16.658, pAlt=0.9)
round(as.numeric(tempRes$pCrossing), 4)
round(as.numeric(tempRes$b*tempRes$pCrossing), 4)
round(as.numeric(tempRes$vBn), 4)
round(as.numeric(tempRes$rBn), 4)
round(as.numeric(tempRes$power), 4)-0.6395

tempRes <- computeVilleTerms(150, 20, pAlt=0.6)
round(as.numeric(tempRes$pCrossing), 4)
round(as.numeric(tempRes$b*tempRes$pCrossing), 4)
round(as.numeric(tempRes$vBn), 4)
round(as.numeric(tempRes$rBn), 4)
round(as.numeric(tempRes$power), 4)




tempRes <- computeVilleTerms(150, 20, pAlt=0.6)
round(as.numeric(tempRes$pCrossing), 4)
round(as.numeric(tempRes$b*tempRes$pCrossing), 4)
round(as.numeric(tempRes$vBn), 4)
round(as.numeric(tempRes$rBn), 4)
# round(as.numeric(tempRes$power), 4)

termsPlot <- matrix(nrow=6, ncol=4)


