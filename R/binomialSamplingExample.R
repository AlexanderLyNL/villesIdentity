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



naiveColour <- "#E31A1CE6"
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

# 0.9 vs 0.5  -------
theta0 <- 0.5
theta1 <- 0.9
# theta1 <- 0.52

nSim <- 1e4
n <- 50
K <- 200

bCandidates <- seq(13, 17, by=0.01)
# bCandidates <- seq(12, 19, by=0.1)


estimatesPVille <- matrix(nrow=K, ncol=length(bCandidates))
estimatesPRaw <- matrix(nrow=K, ncol=length(bCandidates))
estimatesNoBound <- matrix(nrow=K, ncol=length(bCandidates))


pb <- txtProgressBar(min = 0, max = K, style = 3) # style 3 shows percentage

for (k in 1:K) {
  logEValues <- matrix(nrow=nSim, ncol=n)

  set.seed(k)
  for (i in 1:nSim) {
    x <- rbinom(n, size=1, prob=theta0)

    logEValues[i, ] <- cumsum(dbinom(x, size=1, prob=theta1, log=TRUE)-
                                dbinom(x, size=1, prob=theta0, log=TRUE))
  }


  eValues <- exp(logEValues)

  for (iter in seq_along(bCandidates)) {
    kip <- processNullSamplePaths(
      eValues, b=bCandidates[iter])

    estimatesPVille[k, iter] <- kip$pB
    estimatesPRaw[k, iter] <- kip$rawP
    estimatesNoBound[k, iter] <- kip$pB+kip$rB/kip$b
  }

  setTxtProgressBar(pb, k)
}
close(pb)
# Result -----
# load("/Users/alexly/dropbox/notes/villeVictor/code/binomSim09.RData")
# list2env(result, envir = .GlobalEnv)


pFree <- apply(estimatesNoBound, 2, mean)
pCommit <- apply(estimatesPVille, 2, mean)
pNaive <- apply(estimatesPRaw, 2, mean)

min(c(pFree, pCommit, pNaive))
max(c(pFree, pCommit, pNaive))

myName <- "exampleBinomial090"
pdf(paste0(myName, ".pdf"), width=pdfWidth, height=pdfHeight)
graphics::par(cex.main=1.5, mar=c(5, 8.3, 2, 0)+0.1,
              mgp=c(3.5, 1, 0), cex.lab=2.5,
              font.lab=2, cex.axis=1.3, bty="n", las=1)

plot(bCandidates, pComit,
     lwd=myLwd, xlab="",
     ylab="", col=NULL, xaxt="n", yaxt="n",
     bty="n", ylim=c(0.048, 0.054))

lines(c(bCandidates[1], bCandidates[length(bCandidates)]), c(0.05, 0.05),
      lwd=myLwd, col="lightgrey", lty=2)

lines(c(16.658372, 16.658372), c(0, 1), col="grey", lwd=myLwd, lty=2)

# lines(c(16, 16), c(0, 1),
#       lwd=myLwd, col="grey", lty=2)

# lines(bCandidates, pFree,
#       lwd=1.5*myLwd, col=eColours[1])


lines(bCandidates, pCommit,
      lwd=myLwd, col=commitColour)


lines(bCandidates, pNaive,
      col=naiveColour, lwd=myLwd)


axis(1, at=2*(5:10))
axis(2)
mtext("Candidate threshold b", side = 1,
      line = 2.5, las = 1, cex = cexFactor*myCex, padj=0.4)


mtext(latex2exp::TeX(r"($P(T_{b} \leq 50)$)"), side = 2, line = 4, las = 0,
      cex = cexFactor*myCex, adj=0.5, padj=-0.4)
dev.off()



bCandidates[min(which(pNaive <=0.05))]
bCandidates[min(which(pComit <=0.05))]
bCandidates[min(which(pFree <=0.05))]



# Save --------

aap <- 9.98

kip <- processNullSamplePaths(
  eValues, b=15)

kip$rawP
kip$pB
kip$pB+kip$rB/kip$b

mean(estimatesNoBound[, which(bCandidates==aap)])
sd(estimatesNoBound[, which(bCandidates==aap)])

mean(estimatesPRaw[, which(bCandidates==aap)])
sd(estimatesPRaw[, which(bCandidates==aap)])


# result <- list(estimatesNoBound=estimatesNoBound, estimatesPRaw=estimatesPRaw,
#                estimatesPVille=estimatesPVille, theta0=theta0, theta1=theta1,
#                nSim=nSim, n=n, K=K, bCandidates=bCandidates)
#
# save(result, file="binomSim09.RData")



# load("/Users/alexly/Dropbox/Notes/villeVictor/code/binomSim09.RData")
# list2env(result, globalenv())

computeCrossingProbGeneralized(50, 16.658372, 0.9)

# Power ------

fpt <- matrix(nrow=K, ncol=nSim)

pb <- txtProgressBar(min = 0, max = K, style = 3) # style 3 shows percentage

fpt2 <- matrix(nrow=K, ncol=nSim)

b <- 20
b2 <- 16.658372

stoppedMatrix <- stoppedMatrix2 <- matrix(0, nrow=K, ncol=nSim)

for (k in 1:K) {
  logEValues <- matrix(nrow=nSim, ncol=n)

  set.seed(k)
  for (i in 1:nSim) {
    x <- rbinom(n, size=1, prob=theta1)

    logEValues[i, ] <- cumsum(dbinom(x, size=1, prob=theta1, log=TRUE)-
                                dbinom(x, size=1, prob=theta0, log=TRUE))

    iets <- min(which(logEValues[i, ] >= log(b)))

    if (is.finite(iets)) {
      fpt[k, i] <- iets
      stoppedMatrix[k, i] <- 1
    } else {
      fpt[k, i] <- n
      stoppedMatrix[k, i] <- 0
    }

    iets <- min(which(logEValues[i, ] >= log(b2)))

    if (is.finite(iets)) {
      fpt2[k, i] <- iets
      stoppedMatrix2[k, i] <- 1
    } else {
      fpt2[k, i] <- n
      stoppedMatrix2[k, i] <- 0
    }
  }

  setTxtProgressBar(pb, k)
}


mean(fpt)
mean(fpt2)

mean(stoppedMatrix)
mean(stoppedMatrix2)


# # 0.6 vs 0.5  -------
# theta0 <- 0.5
# theta1 <- 0.6
#
#
