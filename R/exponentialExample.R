# Colours ----


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

naiveColour <- "#E31A1CE6"
universalColour <- "#1F78B4E6"
committedColour <- "#BD8E17FF"


processNullSamplePaths <- function(eValueMatrixNull, b, nEnd=NULL) {
  nSim <- dim(eValueMatrixNull)[1]

  if (is.null(nEnd))
    nEnd <- dim(eValueMatrixNull)[2]

  stopVector <- integer(nSim)
  whenStop <- integer(nSim)

  for (i in seq_along(stopVector)) {
    kip <- suppressWarnings(min(which(eValueMatrixNull[i, 1:nEnd] >= b)))

    if (is.finite(kip)) {
      stopVector[i] <- 1
      whenStop[i] <- kip
    }
  }

  stoppedIndeces <- which(stopVector==1)
  overShootObs <- numeric(length(stoppedIndeces))

  for (i in seq_along(stoppedIndeces)) {
    stoppedPath <- stoppedIndeces[i]
    stoppedN <- whenStop[stoppedPath]

    overShootObs[i] <- eValueMatrixNull[stoppedPath, stoppedN]-b
  }

  vB <- sum(overShootObs)/nSim

  notStoppedIndeces <- which(stopVector==0)

  rB <- sum(eValueMatrixNull[notStoppedIndeces, nEnd])/nSim

  pB <- (1-vB-rB)/b

  rawP <- mean(stopVector)
  return(list(pB=pB, rawP=rawP, rB=rB, vB=vB, b=b))
}

# Expo sim -------

lambda0 <- 1
lambda1 <- 0.8

nSim <- 1e4
n <- 500
K <- 100

bCandidates <- seq(15, 17, by=0.01)


estimatesPVille <- matrix(nrow=K, ncol=length(bCandidates))
estimatesPRaw <- matrix(nrow=K, ncol=length(bCandidates))
estimatesNoBound <- matrix(nrow=K, ncol=length(bCandidates))


pb <- txtProgressBar(min = 0, max = K, style = 3)

for (k in 1:K) {
  logEValues <- matrix(nrow=nSim, ncol=n)

  set.seed(k)
  for (i in 1:nSim) {
    x <- rexp(n, rate=lambda0)

    logEValues[i, ] <- cumsum(dexp(x, rate=lambda1, log=TRUE)-
                                dexp(x, rate=lambda0, log=TRUE))
  }

  eValues <- exp(logEValues)

  for (iter in seq_along(bCandidates)) {
    kip <- processNullSamplePaths(
      eValues, b=bCandidates[iter])

    estimatesPVille[k, iter] <- kip$pB
    estimatesPRaw[k, iter] <- kip$rawP
    estimatesNoBound[k, iter] <- kip$pB+kip$rB/kip$b
  }

  setTxtProgressBar(pb, k) # Update progress
}
close(pb)


# Results -----

pFree <- apply(estimatesNoBound, 2, mean)
pCommit <- apply(estimatesPVille, 2, mean)
pNaive <- apply(estimatesPRaw, 2, mean)



myName <- "exampleExponential"
pdf(paste0(myName, ".pdf"), width=pdfWidth, height=pdfHeight)
graphics::par(cex.main=1.5, mar=c(5, 8.3, 2, 0)+0.1,
              mgp=c(3.5, 1, 0), cex.lab=2.5,
              font.lab=2, cex.axis=1.3, bty="n", las=1)

plot(bCandidates, pFree,
     lwd=myLwd, xlab="",
     ylab="", col=NULL, xaxt="n", yaxt="n",
     bty="n")

lines(c(bCandidates[1], bCandidates[length(bCandidates)]), c(0.05, 0.05),
      lwd=myLwd, col="lightgrey", lty=2)

lines(c(16, 16), c(0, 1),
      lwd=myLwd, col="grey", lty=2)

lines(bCandidates, pFree,
      lwd=myLwd, col=eColours[1])

lines(bCandidates, pCommit,
      lwd=myLwd, col=committedColour)

lines(bCandidates, pNaive,
      col=naiveColour, lwd=myLwd)

axis(1)
axis(2)
mtext("Candidate threshold b", side = 1,
      line = 2.5, las = 1, cex = cexFactor*myCex, padj=0.4)
mtext(latex2exp::TeX(r"($P(T_{b} \leq 50)$)"), side = 2, line = 4, las = 0,
      cex = cexFactor*myCex, adj=0.5, padj=-0.4)
dev.off()


bCandidates[min(which(pNaive <=0.05))]
bCandidates[min(which(pComit <=0.05))]
bCandidates[min(which(pFree <=0.05))]
