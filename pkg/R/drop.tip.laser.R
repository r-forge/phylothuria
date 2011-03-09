
##' Drop tips and prepare associated species richness table for analysis in laser
##'
##' This function returns a list that contains a tree (as a \code{phylo} object) and
##' a data frame with the species richness associated with the tips of the tree, after
##' pruning some tips in the original tree object.
##' Note that each element in the list of \code{toDrop} will be passed to the function
##'   \code{\link[ape]{drop.tip}} from the ape package. Also note that contrary to
##' \code{\link[laser]{getTipdata}}, \code{drop.tip.laser} only accept a data frame for
##' the species richness data.
##' @title drop.tip.laser: prepare pruned tree for analysis with the laser package.
##' @param tr a phylogenetic tree stored as a phylo object
##' @param toDrop a list where each element contains at least one tip label that need to
##'  be dropped.
##' @param new.name a vector of tip labels which will replace the dropped tips. The length
##'   of this vector must be the same as the length of \sQuote{toDrop}.
##' @param listSp a data frame with a single column, where the row names are the names of
##'   the tip labels of the original object \sQuote{tr}.
##' @return a list with two elements:
##'   \item{phy}{a phylogenetic tree where the specified tips have
##'   been dropped}
##'   \item{tipdata}{a data frame with a single column, where the number of species have
##'    been modified to account for dropped tips}
##' @examples
##'   data(skinktree)
##'   data(skinkdiversity)
##'   tipsToDrop <- list(c("Hemiergis", "Glaphyromorphus"), c("Lerista", "Ctenotus"))
##'   dropSkink <- drop.tip.laser(tr=skinktree, toDrop=tipsToDrop, new.name=c("H+G", "L+C"),
##'                               listSp=skinkdiversity)
##'   dropSkinkLaser <- getTipdata(dropSkink$tipdata, dropSkink$skinkdiversity)
##' @seealso \code{\link[laser]{getTipdata}}, \code{\link[ape]{drop.tip}}
##' @author François Michonneau
drop.tip.laser <- function(tr, toDrop, new.name="newTipName", listSp) {

  stopifnot(class(tr) == "phylo")
  stopifnot(length(toDrop) > 1)
  stopifnot(length(toDrop) == length(new.name))

  tmpTr <- tr
  tmpList <- listSp
  
  for (i in 1:length(toDrop)) {
    tmpTr <- drop.tip(tmpTr, tip=toDrop[[i]][-1])
    tmpTr$tip.label[match(toDrop[[i]][1], tmpTr$tip.label)] <- new.name[i]
    tmpList <- tmpList[- match(toDrop[[i]][-1], rownames(tmpList)),, drop=FALSE]
    tmpList[match(toDrop[[i]][1], rownames(tmpList)), "sp"] <-
      sum(listSp[match(toDrop[[i]], rownames(listSp)), "sp"])
    rownames(tmpList)[match(toDrop[[i]][1], rownames(tmpList))] <- new.name[i]
  }
  
  list(phy=tmpTr, tipdata=tmpList)
}
