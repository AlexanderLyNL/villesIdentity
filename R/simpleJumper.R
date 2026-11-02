sequentialRanks <- Rcpp::cppFunction('
NumericVector conformal_p_cpp(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double greater = 0.0, equal = 0.0;
    for (int j = 0; j <= i; ++j) {        // include i itself
      if (x[j] < x[i]) greater += 1.0;
      else if (x[j] == x[i]) equal += 1.0;
    }
    double theta = R::runif(0.0, 1.0);    // one theta per step, not per tie
    out[i] = (greater + theta * equal) / (i + 1);
  }
  return out;
}
')

sequentialRanksTop <- Rcpp::cppFunction('
NumericVector conformal_p_cpp(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double greater = 0.0, equal = 0.0;
    for (int j = 0; j <= i; ++j) {        // include i itself
      if (x[j] > x[i]) greater += 1.0;
      else if (x[j] == x[i]) equal += 1.0;
    }
    double theta = R::runif(0.0, 1.0);    // one theta per step, not per tie
    out[i] = (greater + theta * equal) / (i + 1);
  }
  return out;
}
')

computeRollingPValues <- function(data, top=TRUE)
  if (top) sequentialRanksTop(data) else sequentialRanks(data)

# Class definition for Simple Jumper Martingale Test
simpleJumper <- function(
    # epsilons = c(-1.0, -0.5, 0.0, 0.5, 1.0),
    epsilons = c(-1, 0, 1),
    J = 0.01) {
  k <- length(epsilons)

  # Distribute initial wealth (S_0 = 1) equally across components
  wealthDistribution <- rep(1 / k, k)
  capitalHistory <- c(1.0)

  # Return an environment acting as an object closure to manage mutable state
  self <- new.env()
  self$epsilons <- epsilons
  self$k <- k
  self$J <- J
  self$wealthDistribution <- wealthDistribution
  self$capitalHistory <- capitalHistory

  self$update <- function(pValue) {
    # 1. Apply the betting functions to the current p-value
    bettingMultipliers <- 1 + self$epsilons * (pValue - 0.5)

    # 2. Update wealth allocation per component
    updatedWealth <- self$wealthDistribution * bettingMultipliers

    # 3. Calculate total capital after the bet
    totalCapital <- sum(updatedWealth)

    # 4. Redistribute capital to allow components to 'jump' (mix step)
    if (totalCapital > 0) {
      # Each component keeps (1-J) of its value and gains J of the average value
      self$wealthDistribution <- (1 - self$J) * updatedWealth + (self$J * totalCapital / self$k)
    } else {
      self$wealthDistribution <- rep(0, self$k)
    }

    self$capitalHistory <- c(self$capitalHistory, totalCapital)
    return(totalCapital)
  }

  return(self)
}

# Helper to compute rolling conformal p-values

# --- Execution Example ---
set.seed(1)

# Generate stationary data (exchangeable)
stationaryData <- rnorm(200, mean = 0, sd = 1)

# Introduce a mean location shift (non-exchangeable distribution shift)
shiftedData <- rnorm(200, mean = 6, sd = 1)
fullDataSet <- c(stationaryData, shiftedData)

# 1. Transform data stream into independent conformal p-values
pVals <- computeRollingPValues(fullDataSet)

# 2. Run the Simple Jumper Test
jumperTest <- simpleJumper(J = 0.01)

alpha <- 0.05
threshold <- 1.0 / alpha  # Rejection threshold by Ville's inequality (20.0)
rejected <- FALSE

cat(sprintf("%-6s\t%-12s\t%-10s\t%-18s\n", "Step", "Data Value", "p-value", "Martingale Capital"))
cat(paste0(rep("-", 60), collapse = ""), "\n")

for (step in 1:length(pVals)) {
  p <- pVals[step]
  capital <- jumperTest$update(p)

  # Print every 40 steps, or immediately when the threshold is crossed
  if ((step - 1) %% 40 == 0 || (capital > threshold && !rejected)) {
    cat(sprintf("%-6d\t%-12.3f\t%-10.4f\t%-18.4f\n",
                step - 1,
                fullDataSet[step], p, capital))
  }

  if (capital > threshold && !rejected) {
    cat(sprintf("\n[!] EXCHANGEABILITY REJECTED at step %d! Capital (%.2f) exceeded threshold (%.2f).\n\n",
                step - 1, capital, threshold))
    rejected <- TRUE
  }
}
