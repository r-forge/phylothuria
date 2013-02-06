##'
##' ##' Graft unsampled taxa onto a tree by simulating random trees.
##'
##' Based on a data frame containing the number of species at the tip of a phylogeny,
##' this function will generate random phylogenies for each group of missing species
##' and graft them onto the tree. They are not grafted directly at the tip of the existing
##' tree. Instead a random fraction of the existing branch length is removed and the
##' simulated tree, representing the missing taxa, is scaled to have the same depth (length)
##' as the fraction that is removed. In the end, an ultrametric tree with the same depth (length)
##' as the original object is generated.
##' @title Graft unsample taxa onto a tree.
##' @param tr a phylogenetic tree stored as a \code{phylo} object.
##' @param listSp a single-column data frame which contains the number of species unsampled
##' at each tip. The labels of the tree \code{tr} and the row names of \code{listSp} need to
##' match.
##' @return a phylogenetic tree stored as a \code{phylo} object. The tips for which species have
##' been added too are appended a number (starting at 1) for each clade.
##' @author François Michonneau
##' @examples
##' tr <- rcoal(10)
##' listSp <- data.frame(sp=sample(1:10), row.names=tr$tip.label)
##' trWithGraft <- graftMissingTaxa(tr, listSp)
##' par(mfrow=c(1,2))
##' plot(tr)
##' plot(trWithGraft)
graftMissingTaxa <- function(tr, listSp) {

  if (!inherits(tr, "phylo"))
    stop("\'tr\' must be an object of class \'phylo\'")
  if (!is.ultrametric(tr))
    warning("This function is intended to be used on ultrametric trees.")
  
  spToAdd <- subset(listSp, listSp[,1] > 1)

  for (i in 1:nrow(spToAdd)) {
    trP4 <- as(tr, "phylo4")
    tipToAdd <- spToAdd[i,1]
    tipName <- rownames(spToAdd)[i]
    rTr <- as(rcoal(tipToAdd), "phylo4")   # draw random tree (coal to be ultrametric)
    tipLabels(rTr) <- paste(tipName, 1:tipToAdd, sep="_")
    distToTip <- edgeLength(trP4, getNode(trP4, tipName))
    newDist <- runif(1) * distToTip
    totalDist <- distToTip - newDist
    edgeLength(rTr) <- edgeLength(rTr) * totalDist / sumEdgeLength(rTr, ancestors(rTr, 1, "ALL"))
    edgeLength(trP4)[getEdge(trP4, tipName)] <- unname(newDist)
    tr <- bind.tree(as(trP4, "phylo"), as(rTr, "phylo"), match(tipName, tr$tip.label))
  }
  stopifnot(is.ultrametric(tr))
  tr
}
