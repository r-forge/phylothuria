### function barMonophyletic
barMonophyletic <- function(groupLabel, tree, cex.plot, cex.text=.8, include.tip.label=FALSE,
                            extra.space=0, coef.space=1, draw=TRUE, text.offset=1.02, font=1,
                            font.col=1, seg.col=1, srt=0, bar.at.tips=FALSE) {

  font.col <- rep(font.col, length(groupLabel))
  font <- rep(font, length(groupLabel))
  seg.col <- rep(seg.col, length(groupLabel))
  
  getTipOrderPlot <- function(tr) {
    tr$tip.label[tr$edge[tr$edge[, 2] <= length(tr$tip.label), 2]]
  }

  findMaxWidth <- function(tr, grpLbl, cex.plot, include.tip.label) {
    offLabel <- ifelse(include.tip.label, strwidth(tr$tip.label, cex=cex.plot), 0)
    max(branching.times(tr)) + extra.space + offLabel 
  }
  
  grpID <- vector("list", length(groupLabel))
  names(grpID) <- groupLabel
  
  for (i in 1:length(groupLabel)) {
    xx <- grep(groupLabel[i], tree$tip.label, value=TRUE)
    grpID[[i]] <- match(xx, getTipOrderPlot(tree))
  }

  maxTip <- length(tree$tip.label)
 
  barPos <- findMaxWidth(tr=tree, groupLabel, cex.plot=cex.plot, include.tip.label=include.tip.label) * coef.space
 
  if (draw) {
    for (i in 1:length(groupLabel)) {
      if (length(grpID[[i]]) == 1) {
        rg <- c(grpID[[i]] - 0.2, grpID[[i]] + 0.2)      
      }
      else {
        rg <- range(grpID[i])
      }
      if (bar.at.tips) {
        tips <- range(grpID[i])
        offLabel <- ifelse(include.tip.label, max(strwidth(tree$tip.label[tips], cex=cex.text)), 0)
        lastPP <- get("last_plot.phylo", envir=.PlotPhyloEnv)
        barPos <- max(lastPP$xx[tips]) + offLabel + extra.space
      }
      text(barPos*text.offset, mean(rg), groupLabel[i], adj=0, font=font[i], col=font.col[i], srt=srt)
      segments(barPos, rg[1], barPos, rg[2], lwd=3, col=seg.col[i])
    }
  }
  
  toRet <- barPos + max(strwidth(groupLabel, cex=cex.text))
  names(toRet) <- "TotalWidth"
  invisible(toRet)
}
