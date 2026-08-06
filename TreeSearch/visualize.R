# Print current working directory, which should contain the scripts,
# matrix and trees.
getwd()

# If getwd() does not contain the relevant files, set wd to working directory
wd <- "./"
outgroup <- c("") # Specify taxa on which to root tree

# Number of trees to map in the tree-space panel.  Searches now run with
# poolMaxSize = 5000, and that panel is the only quadratic step: the distance
# matrix alone is ~21 min and 200 MB at 5000 trees, before cmdscale() attempts a
# 5000 x 5000 eigendecomposition.  The consensus is linear, so it still uses
# every tree.
nMap <- 100

source(paste0(wd, "/common.R"))
source(paste0(wd, "/plot.R"))

latest <- LatestMatrix(wd)
dat <- ReadAsPhyDat(latest)
if (outgroup == "") {
  outgroup <- names(dat)[1]
}
treeFiles <- list.files(
  path = wd,
  pattern = paste0(".+_", sub("^.*/", "", latest), ".trees"),
  full.names = TRUE
)

for (treeFile in treeFiles) {
  trees <- read.nexus(treeFile)
  
  # Ignore outgroup taxa that aren't in tree
  og <- intersect(outgroup, TipLabels(trees)[[1]])
  if (length(og)) {
    # Root trees on outgroup
    trees <- RootTree(trees, og)
  }
  rogues <- Rogue::QuickRogue(trees, p = 1)
  cons <- SortTree(ConsensusWithout(trees, rogues[-1, "taxon"]))
  
  pdf(gsub(".trees", ".pdf", treeFile, fixed = TRUE), 
      width = 8, height = 10)
  
  ColPlot(cons, ec = "black")
  if (nrow(rogues) > 1) {
    legend("topleft", rogues[-1, "taxon"], bty = "n", lty = 2)
  }
  k <- KValue(treeFile)
  # Resolve before scoring: with collapse = TRUE (the default) a search returns
  # a mixed-resolution set, and TreeLength() requires a binary tree.  Resolving
  # a soft polytomy does not change the parsimony score, so this reports the
  # same number whichever tree happens to come first.
  legend(
    "topright",
    c(
      sub("^(?:.*/)*([^/_]+)_.+", "\\1", treeFile, perl = TRUE),
      paste("Score:", signif(TreeLength(
        Preorder(MakeTreeBinary(trees[[1]])), dat, concavity = k
      )))
    ),
    bty = "n" # No bounding box
  )


  # Map a topologically diverse subset, not a random one: the number of MPTs in
  # a region of tree space reflects the density of the parsimony landscape, not
  # support, so a random draw over-represents broad plateaux and can miss an
  # isolated island entirely.  effort = 1 (FarFirst) is matrix-free, so it never
  # builds the N x N matrix the default tier would -- ~5 min at 5000 trees, and
  # it achieves a wider minimum separation than the best of 200 random draws.
  mapTrees <- if (length(trees) > nMap) {
    WideSample(trees, nMap, effort = 1)
  } else {
    trees
  }

  distances <- TreeDist::ClusteringInfoDistance(mapTrees)
  whenHit <- gsub("_\\d+$", "", names(mapTrees), perl = TRUE)
  firstHit <- table(whenHit)
  searchStages <- length(firstHit)
  map <- cmdscale(distances, k = 3)
  cols <- hcl.colors(searchStages, alpha = 0.8)
  presOrder <- c("seed", "start", paste0("ratch", 1:10000), "final")
  presOrder <- c(presOrder, setdiff(names(firstHit), presOrder))
  treeCols <- cols[match(whenHit, intersect(presOrder, whenHit))]

  # Prepare plotting area
  par(mar = rep(0, 4))
  plot(map, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
  
  # Add minimum spanning tree
  TreeTools::MSTEdges(distances, plot = TRUE, map[, 1], map[, 2],
                      col = '#00000030', lty = 2)
  
  # Connect trees by order found
  lines(map[, 1], map[, 2], col = "#ffccaa", lty = 1)
  
  # Add points
  TreeDist::Plot3(map,
                  col = treeCols,
                  pch = 16, cex = 2,
                  add = TRUE)
  
  # Add legends
  legend("topright", 
         intersect(presOrder, names(firstHit)),
         col = cols, pch = 16, bty = "n")
  legend("topleft",
         c("Minimum spanning tree (mapping distortion)",
           "Order in which trees found"),
         lty = c(2, 1),
         col = c("#00000030", "#ffccaaaa"),
         bty = "n")
  if (length(mapTrees) < length(trees)) {
    legend("bottomleft",
           paste0("WideSample: ", length(mapTrees), " of ", length(trees),
                  " trees (Max-Min diverse)"),
           bty = "n")
  }

  dev.off()
}
