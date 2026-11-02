nSim <- 1e6
n0 <- 1000
n1 <- 1000

eValueMatrixNull <- matrix(nrow=nSim, ncol=n0+n1)

progressBar <- txtProgressBar(min = 0, max = nSim, style = 3) # style 3 shows percentage
for (i in 1:nSim) {
  set.seed(i)
  # Generate stationary data (exchangeable)
  stationaryData <- rnorm(n0, mean = 0, sd = 1)

  # Introduce a mean location shift (non-exchangeable distribution shift)
  shiftedData <- rnorm(n1, mean = 0, sd = 1)
  fullDataSet <- c(stationaryData, shiftedData)

  # 1. Transform data stream into independent conformal p-values
  pVals <- computeRollingPValues(fullDataSet)

  # 2. Run the Simple Jumper Test
  jumperTest <- simpleJumper(J = 0.01)


  for (step in 1:length(pVals)) {
    p <- pVals[step]
    capital <- jumperTest$update(p)
    eValueMatrixNull[i, step] <- capital
  }
  setTxtProgressBar(progressBar, i)
}
close(progressBar)

# save(eValueMatrixNull, file="simpleJumperEValueMatrixNull1e6n2000.RData")


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


processTimeSamplePaths <- function(eValueMatrixNull, b, nEnd=NULL){
  nSim <- dim(eValueMatrixNull)[1]

  if (is.null(nEnd))
    nEnd <- dim(eValueMatrixNull)[2]

  rVec <- vVec <- numeric(n)

  trackIndicator <- integer(nSim)

  for (i in 1:n) {
    crossIndex <- which(eValueMatrixNull[, i] >= b)

    if (length(crossIndex) > 0) {
      someIndeces <- intersect(crossIndex, which(trackIndicator==0))

      trackIndicator[someIndeces] <- 1

      if (length(someIndeces) > 0){
        vVec[i] <- sum(eValueMatrixNull[someIndeces, i]-b)
      }
    }

    rVec[i] <- sum(eValueMatrixNull[trackIndicator==0, i])
  }

  vB <- 1/nSim*cumsum(vVec)
  rB <- 1/nSim*(rVec)

  pB <- (1-vB-rB)/b

  result <- list(vB=vB, rB=rB, pB=pB)
  return(result)
}
