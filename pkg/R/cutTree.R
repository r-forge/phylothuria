##' cutTree cuts a phylogenetic tree to a given age.
##'
##' This function cuts a tree to a given age by dropping all tips that
##' originated after that age and adjusts the terminal edge lengths. This function
##' only accept ultrametric trees.
##' @title cutTree
##' @param tr a phylogenetic tree (as a \code{phylo} object)
##' @param ageCut the age at which to cut the tree
##' @param classResults class of the returned object, can be either \code{phylo} or \code{phylo4}
##' @return a phylogenetic tree as a \code{phylo4} (default) or a \code{phylo} object.
##' @example
##'   data(skinktree)
##'   skinktree$node.label <- NULL # to avoid warning
##'   cutTree(skinktree, 10, classResults="phylo")
##' @author François Michonneau
cutTree <- function(tr, ageCut, classResults=c("phylo4", "phylo")) {

  classResults <- match.arg(classResults)
  stopifnot(is.ultrametric(tr))
  trP4 <- as(tr, "phylo4")
  newTrP4 <- trP4
  if (! all(is.na(nodeLabels(newTrP4)))) warning("all node labels are going to be erased")
  tr$node.label <- NULL
  
  bT <- branching.times(tr)
  ageRoot <- max(bT) - ageCut
  nodesToKeep <- names(bT)[bT > ageCut]
  nodesToKeep <- as.integer(nodesToKeep)

  desc <- descendants(trP4, nodesToKeep, "children")
  desc <- unique(unlist(desc))
 
  nodeLabels(newTrP4) <- as.character(nodeId(newTrP4, "internal"))
  for (i in 1:length(desc)) {
    if (desc[i] <= nTips(trP4)) next   # skip node if it's a tip
    ch <- descendants(trP4, desc[i], "children")
    if (all(sapply(ch, function(x) x %in% desc))) next # skip node if it's not "terminal"  
    toDrop <- descendants(trP4, desc[i])[-1]  
    newTrP4 <- subset(newTrP4, tips.exclude=names(toDrop))
  }

  for (i in 1:nTips(newTrP4)) {
    anc <- ancestor(newTrP4, i)
    dRoot <- sumEdgeLength(newTrP4, ancestors(newTrP4, anc, "ALL"))
    nextEdgeLength <- ageRoot - dRoot
    if (nextEdgeLength < 0) browser()
    edgeLength(newTrP4)[getEdge(newTrP4, i)] <- nextEdgeLength
  }
  
  if (classResults == "phylo4")
    newTrP4
  else as(newTrP4, "phylo")
}
