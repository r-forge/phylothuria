##' \code{dropExtinct} prunes from a phylogenetic tree all the tips that are extinct,
##' that is the taxa for which the distance from the root is less than the maximal distance
##' from the root.
##'
##' This function is intended to be used after simulating trees from a birth-death process
##' such as \code{\link[geiger]{birthdeath.tree}} to recover only the taxa that are still
##' extant. This function differs from \code{\link[geiger]{prune.extinct.taxa}} is that
##' \code{dropExtinct} takes a \code{phylo4} object.
##' @title prunes extinct taxa
##' @param tr a phylogenetic tree stored as a \code{phylo4} or \code{phylo4d} object.
##' @param classResult the class of the object to be return (\code{phylo4} for \code{phylo4}
##'   or \code{phylo4d} objects, or \code{phylo})
##' @param tol tolerance in rounding error to avoid taxa to be pruned because of rounding error.
##' @return a phylogenetic tree as a \code{phylo4(d)} or a \code{phylo} object.
##' @export
##' @examples
##'   rTr <- birthdeath.tree(b=.4, d=.2, time=15)
##'   rTr <- as(rTr, "phylo4")
##'   pTr <- dropExtinct(rTr)
##'   nTips(rTr)
##'   nTips(pTr)
##' @seealso \code{\link[geiger]{prune.extinct.taxa}}, \code[geiger]{birthdeath.tree}
##' @author François Michonneau
dropExtinct <- function(tr, classResult=c("phylo4", "phylo"),
                        tol=.Machine$double.eps^0.5) {

  classResult <- match.arg(classResult)
  stopifnot(class(tr) == "phylo4" || class(tr) == "phylo4d")
  
  ageTip <- depthTips(tr)
  trDepth <- max(ageTip)
  
  toDrop <- (1:nTips(tr))[(max(ageTip) - ageTip) > tol]
  trRes <- subset(tr, tips.exclude=toDrop)

  if (nTips(trRes) < 2) warning("Less than 2 tips in your tree!")  
  if (classResult == "phylo")
    as(trRes, "phylo")
  else {
    trRes@annote$treeDepth <- max(ageTip)
    trRes
  }
}
