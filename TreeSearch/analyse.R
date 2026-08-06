source("common.R")

kValues <- c(10, 40, 6) # Concavity constants for implied weighting


# Load data from locally downloaded matrix
latest <- LatestMatrix()
message("* Reading ", latest)
dat <- ReadAsPhyDat(latest)

resultsFile <- ResultsFile(latest, "ew")

best <- MaximizeParsimony(dataset = dat, effort = -1)
write.nexus(best, file = resultsFile)

best <- MaximizeParsimony(dataset = dat, tree = best, effort = 1,
                          targetHits = NTip(dat) / 2)
write.nexus(best, file = resultsFile)


for (k in kValues) {
  resultsFile <- ResultsFile(latest, "iw", k)
  best <- MaximizeParsimony(
    dataset = dat,
    concavity = k,
    targetHits = NTip(dat) / 2,
    effort = 1
  )
  write.nexus(best, file = resultsFile)
}
