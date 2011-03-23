##' \code{ltt.2rates} extracts branching times and number of species
##' from simulated phylogenies (with the function \code{sim2rate}) in
##' preparation for lineage through time plots.
##'
##' This function extracts the number of species and branching times
##' from a phylogeny simulated under a 2-rate diversification
##' model. The clade simulated under the new diversification rate will
##' be identified based on a character string provided by the user
##' that is passed to the function \code{\link[phylobase]{MRCA}} and evaluated.
##' @title Extract number of species and branching times from
##' simulated trees (through the function \code{\link{replicate}}) in
##' preparation for lineage-through time plot.
##' @param trReplicates the object resulting from a call to the
##' function \code{replicate} where each element contains a
##' \code{phylo} object.
##' @param mrca.newRate a character that will be evaluated to identify
##' the MRCA node for the clade simulated under one of the rate of
##' diversification (by convention the \sQuote{new} rate).
##' @return a list containing the following elements for which the
##' length will be equal to the number of replicates in
##' \code{trReplicates}.
##' \item{timeNewRate}{a list containing the branching times for the
##' clade simulated under the \sQuote{new} rate in decreasing order.}
##' \item{timeBaseRate}{a list containning the branching times for the
##' clade simulated under the \sQuote{basal} rate in decreasing order.}
##' \item{offsetNbSp}{a list containing the number of species in the
##' \sQuote{basal} part of the tree at the time of the shift in
##' diversification rate.}
##' \item{nSpNewRate}{a list containing the number of species in the
##' clade simulated under the \sQuote{new} rate}
##' \item{nSpBaseRate}{a list containing the number of species in the
##' clade simulated under the \sQuote{basal} rate}
##' @export
##' @author François Michonneau
##' @note This function is fairly slow. See the vignette for some
##' examples on how to use this function.
ltt.2rates <- function(trReplicates, mrca.newRate) {

  nbRep <- ncol(trReplicates)
  timeNewRate <- timeBaseRate <- nSpNewRate <-
    nSpBaseRate <- offsetNbSp <- vector("list", nbRep)
  pb <- txtProgressBar(min=1, max=nbRep, style=3)
  for (i in 1:nbRep) {
    setTxtProgressBar(pb, i)
    trRep <- trReplicates[, i]
    class(trRep) <- "phylo"
    trRepP4 <- as(trRep, "phylo4")
    mRca <- MRCA(trRepP4, eval(parse(text=mrca.newRate)))
    trRepNewRate <- subset(trRepP4, mrca=mRca)
    trRepBaseRate <- prune(trRepP4, tips.exclude=tipLabels(trRepNewRate))
    trRepNewRate <- as(trRepNewRate, "phylo")
    trRepBaseRate <- as(trRepBaseRate, "phylo")
    ageMrca <- branching.times(trRep)[as.character(mRca)]
    timeNewRate[[i]] <- sort(branching.times(trRepNewRate), decreasing=TRUE)
    timeBaseRate[[i]] <- sort(branching.times(trRepBaseRate), decreasing=TRUE)
    offsetNbSp[[i]] <- sum(timeBaseRate[[i]] > timeNewRate[[i]][1])
    nSpNewRate[[i]] <- (1:length(trRepNewRate$tip.label))
    nSpBaseRate[[i]] <- 1:length(trRepBaseRate$tip.label)
  }
  close(pb)
  
  list(nSpNewRate=nSpNewRate, nSpBaseRate=nSpBaseRate,
       timeNewRate=timeNewRate, timeBaseRate=timeBaseRate,
       offsetNbSp=offsetNbSp)
}

##' Create lineage through time (LTT) plots for from the output of \code{ltt.2rates}
##'
##' \code{plot.2rates} takes the output from \code{ltt.2rates} and generates lineage
##' through time plots that illustrates the differences in diversification rates. The
##' first method creates separate curves for each diversification rate. The second method
##' puts the two curves on the same line.
##' @title Lineage through time plots for two diversification rates.
##' @param listRates the output from the \code{ltt.2rates} function.
##' @param new.plot Should a new plot be created (TRUE, default) or should the LTT be superimposed on
##'   an existing plot (FALSE)?
##' @param method Should the lineages appear on separate lines (1, default, recommended) or on the same
##'   line (2).
##' @param colBase color for the lineages simulated under the \sQuote{base} rate.
##' @param colNew color for the lineages simulated under the \sQuote{new} rate.
##' @param ... additional graphical parameters.
##' @return Nothing. Function to be used for its side effect.
##' @export
##' @author François Michonneau
##' @note see vignettes for examples.
plot.2rates <- function(listRates, new.plot=TRUE, method=c(1,2), colBase="blue", colNew="red", ...) {

  nbRep <- length(listRates[[1]])
  method <- match.arg(method)
  
  if (new.plot) {
    maxTime <- max(c(unlist(listRates$timeNewRate), unlist(listRates$timeBaseRate)))
    maxSp <- max(c(unlist(listRates$nSpNewRate), unlist(listRates$nSpBaseRate)))
    plot(-c(maxTime, 0), c(1, maxSp), type="n", ...)
  }

  if (method == 1) {
    for (i in 1:nbRep) {
      lines(-c(listRates$timeBaseRate[[i]], 0), listRates$nSpBaseRate[[i]], type="S", col=colBase)
      lines(-c(listRates$timeNewRate[[i]], 0), listRates$nSpNewRate[[i]], type="S", col=colNew)
    }    
  }
  if (method == 2) {
    for (i in 1:nbRep) {
      lines(-c(listRates$timeNewRate[[i]], 0), listRates$nSpNewRate[[i]]+listRates$offsetNbSp[[i]], type="S", col=colNew)
      lines(-c(listRates$timeBaseRate[[i]], 0), listRates$nSpBaseRate[[i]], type="S", col=colBase)
    }
  }
}
