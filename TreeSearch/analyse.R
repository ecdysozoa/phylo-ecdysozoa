source("common.R")

kValues <- c(10, 40, 6) # Concavity constants for implied weighting


# Load data from locally downloaded matrix
latest <- LatestMatrix()
message("* Reading ", latest)
dat <- ReadAsPhyDat(latest)

resultsFile <- ResultsFile(latest, "ew")

poolMaxSize <- 5000

best <- MaximizeParsimony(dataset = dat, effort = -1)
write.nexus(best, file = resultsFile)

best <- MaximizeParsimony(dataset = dat, tree = best, effort = 1,
                          poolMaxSize = poolMaxSize)
write.nexus(best, file = resultsFile)


for (k in kValues) {
  resultsFile <- ResultsFile(latest, "iw", k)
  best <- MaximizeParsimony(
    dataset = dat,
    concavity = k,
    effort = 1,
    poolMaxSize = poolMaxSize
  )
  write.nexus(best, file = resultsFile)
}
