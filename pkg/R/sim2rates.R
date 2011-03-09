##' \code{sim2rate} simulates a tree under two regimes of diversification
##' starting from a backbone tree. 
##'
##' \code{sim2rate} can be used to simulate phylogenies created under
##' 2 rates of diversification such as the ones that are estimated by
##' \link{\code[laser]{fitNDR_2rate}} from the laser
##' package. Providing a backbone tree (\code{baseTree}) allows the
##' user to control the number of taxa present in the phylogeny at the
##' time of the shift in diversification rate and/or the time elapsed
##' before the shift in diversification occurred. The tips belonging
##' to the clades generated under the \code{newRate} of
##' diversification will be identified by the prefix \sQuote{new} and
##' will be numbered from 1. The tips belonging to the clades
##' generated under the \code{baseRate} will be identified by a letter
##' (from \sQuote{a}) and a number. The clades that share the same
##' letter as a prefix share the same common ancestor which is one of
##' the tip from the backbone tree provided (\code{baseTree}). This
##' means that if the backbone tree has more than 26 tips, the labels
##' originating from these tips will have the prefix \sQuote{NA},
##' which will lead to duplicated labels. If that doesn't fit your
##' needs, let me know. In the event that no speciation event is
##' generated before \code{timeTotal} is reached for a given branch,
##' the tip will keep the original label provided in \code{baseTree}. 
##' @title Simulating phylogenies under 2 rates of diversification.
##' @param baseTree a backbone tree stored as a \code{phylo} object.
##' @param baseRate the base (=background) net rate of diversification for the tree. 
##' @param newRate the new (=elevated) net rate of diversification.
##' @param eps the relative extinction rate (see \link{\code[laser]{fitNDR_2rate}}).
##' @param timeTotal the total depth of the tree (the amount of time between the root
##'   and the tip of tree).
##' @return an ultrametric phylogeny stored as a \code{phylo} object.
##' @author François Michonneau
##' @note This function is fairly slow, as it hasn't been optimized.
##' @examples
##' ## baseTr <- birthdeath.tree(b=.25, d=.1, time.stop=20)
##' ## baseTr <- as(baseTr, "phylo4")
##' ## baseTr <- dropExtinct(baseTr, classResult="phylo")
##' ## simTr  <- sim2rate(baseTr, baseRate=.13, newRate=.44, eps=.46, timeTotal=25)
##'   

sim2rate <- function(baseTree, baseRate, newRate, eps, timeTotal) {
  
  stopifnot(class(baseTree) == "phylo")
  stopifnot(eps < 1 && eps >= 0)

  dBaseRate <- eps * baseRate / (1 - eps)  
  dNewRate  <- eps * newRate / (1 - eps)
  bBaseRate <- baseRate / (1 - eps)
  bNewRate <- newRate / (1 - eps)

  tipGraft <- sample(1:Ntip(baseTree), 1)

  timeShift <- timeTotal - max(branching.times(baseTree))
  stopifnot(timeShift > 0)
  
  depthNewRateTree <- 0
  while (abs(depthNewRateTree - timeShift) > .Machine$double.eps^.5) {
    newRateTree <- birthdeath.tree(b=bNewRate, d=dNewRate, time.stop=timeShift,
                                   return.all.extinct=FALSE)
    newRateTree <- dropExtinct(as(newRateTree, "phylo4"), "phylo")
    depthNewRateTree <- depthTree(newRateTree)
  }

  newRateTree$tip.label <- paste("new", newRateTree$tip.label, sep="")
  newTree <- bind.tree(baseTree, newRateTree, tipGraft)
  newTreeP4 <- as(newTree, "phylo4")
  baseTreeP4 <- as(baseTree, "phylo4")

  ## graft trees at all other tips build with base rate
  for (i in 1:nTips(baseTreeP4)) {
    cat(i, "\n")
    tipGraft <- match(baseTree$tip.label[i], newTree$tip.label)
   
    if (is.na(tipGraft)) next
   
    baseEdLength <- rexp(1, rate=bBaseRate+dBaseRate)
    currentEdLength <- edgeLength(baseTreeP4)[getEdge(baseTreeP4, i)]
    nodeNewTree <- getNode(newTreeP4, tipLabels(baseTreeP4)[i])

    if (baseEdLength < timeShift) {
      ## ageTmpTree: total distance from node where tmpTree is grafted to tip
      ageTmpTree <- timeShift - baseEdLength
      tmpTree <- birthdeath.tree(b=bBaseRate, d=dBaseRate, time.stop=ageTmpTree,
                                 return.all.extinct=FALSE)
      tmpTree <- dropExtinct(as(tmpTree, "phylo4"), "phylo")
      tmpTree$tip.label <- paste(letters[i], tmpTree$tip.label, sep="")
      depthTmpTree <- depthTree(tmpTree)
      ## newEdgeLength
      newEdgeLength <- currentEdLength + baseEdLength + (ageTmpTree - depthTmpTree) 
      edgeLength(newTreeP4)[getEdge(newTreeP4, nodeNewTree)] <- newEdgeLength
      
      newTree <- as(newTreeP4, "phylo")
      newTree <- bind.tree(newTree, tmpTree, nodeNewTree)
    }
    else {  
      edgeLength(newTreeP4)[getEdge(newTreeP4, nodeNewTree)] <- currentEdLength + timeShift
      newTree <- as(newTreeP4, "phylo")
    }
    newTreeP4 <- as(newTree, "phylo4")
  }
  stopifnot(is.ultrametric(newTree))
  stopifnot(is.binary.tree(newTree))
  stopifnot(abs(sumEdgeLength(newTreeP4, ancestors(newTreeP4, 1, "ALL")) - timeTotal) < .Machine$double.eps^.5)
  newTree
}
