
recEnv <- new.env(parent=globalenv())
with(recEnv, i <- 1)
with(recEnv, listNodes <- NA)

record <- function(node) {
  assign("node", node, envir=recEnv)
  with(recEnv, listNodes[i] <- node)
  with(recEnv, i <- i+1)
}

postorder <- function(tr, node) {
  if (nodeType(tr)[node] != "tip") {
    desc <- descendants(tr, node, "children")
    postorder(tr, node=desc[1])
    postorder(tr, node=desc[2])
  }
  record(node)
}

preorder <- function(tr, node) {
  if (nodeType(tr)[node] != "root") {
    anc <- ancestors(tr, node, "parent")
    if (length(anc) == 1) {
      preorder(tr, node=anc)
    }
    else {
      preorder(tr, node=anc[1])
      preorder(tr, node=anc[2])
    }
  }
  record(node)
}



distFromTip <- function(tr, node) {
  lDesc <- descendants(tr, node)
  max(sapply(lDesc, function(x) {
    pth <- .Call("interCpp", descendants(tr, node, type="all"),
                 ancestors(tr, x, type="ALL"))$inter
    sumEdgeLength(tr, pth)
    }))
}

distFromTipFast <- function(tr) {
  storeList <- vector("list", nNodes(tr))
  names(storeList) <- nodeId(tr, "internal")
  for (i in 1:nTips(tr)) {
    with(recEnv, { i <- 1; listNodes <- NA})
    preorder(tr, i)
    lAnc <- with(recEnv, listNodes)
    for (j in 1:length(lAnc)) {
      ancTmp <- lAnc[j:length(lAnc)]
      storeList[[as.character(lAnc[j])]] <- c(storeList[[as.character(lAnc[j])]], sumEdgeLength(tr, lAnc[j:length(lAnc)]))
    }
  }
  sapply(storeList, max)
}

distFromTipCpp <- '
  Rcpp::NumericVector lNod(lNodR);
  Rcpp::IntegerVector lAnc(lAncR);
  Rcpp::NumericVector eLgt(eLgtR);

  std::vector<double> storeRes(lNod.size());
  
  for (int j=0; j < lAnc.size(); j++) {
    std::vector<int*> tmpAnc;
    int beg = lAnc[j];
    tmpAnc.insert(tmpAnc.begin(), beg, lAnc.end());
    for (int k=0; k < tmpAnc.size(); k++) {
      double sumEdgeLength;
      int l = k;
      while (l < tmpAnc.size()) {
         sumEdgeLength + eLgt[l];
         l++;
      }
      storeRes[k] = sumEdgeLength;
    }
  }

  return Rcpp::List::create(storeRes);
'

dCpp <- cfunction(signature(lNodR="numeric", lAncR="integer", eLgtR="numeric"),
                  distFromTipCpp, 
                  Rcpp=TRUE)

findGroups <- function(tr, threshold=.025) {

  stopifnot(class(tr) == "phylo4") # || class(tr) == "phylo4d")
  
  grp <- vector("list", nTips(tr))

  ## find the distance to the tips for each internal node and select the nodes
  ## below the threshold
  ## intNodes <- nodeId(tr, "internal")
  ## lGrp <- sapply(intNodes,
  ##                function(x) distFromTip(tr, x))
  ## names(lGrp) <- intNodes
  ## lGrp <- lGrp[lGrp <= threshold]

  ### New faster? way
  lGrp <- distFromTipFast(tr)
  lGrp <- lGrp[lGrp <= threshold]

  ## find all the descendants for the nodes below the threshold
  descGrp <- sapply(as.numeric(names(lGrp)), function(x) descendants(tr, x))
  
  ## graft all tips to be sure that even singletons are included in results
  descGrp <- c(lapply(tipLabels(tr), function(x) getNode(tr, x)), descGrp)

  ## merge all the identical results
  for (i in 1:length(descGrp)) {
    tmpGrp <- descGrp[sapply(descGrp, function(x) any(descGrp[[i]] %in% x))]
    grp[[i]] <- sort(unique(unlist(tmpGrp)))
  }
  grp <- unique(grp)
  grp <- grp[!sapply(grp, is.null)]
  grp <- sapply(grp, function(x) tipLabels(tr)[x]) # return the tip labels

  ## build a phylo4d object for the results
  dTip <- data.frame(Groups=rep(1:length(grp), sapply(grp, length)))
  rownames(dTip) <- unlist(grp)
  phylo4d(tr, tip.data=dTip)
}

### Example code
## spGrp <- findGroups(trP4, threshold=.025)
## spGrpCopy <- spGrp
## tipLabels(spGrp) <- paste(tipData(spGrp)$Group, tipLabels(spGrp), sep="_")

## grpLbl <- paste("^", 1:max(tipData(spGrp)), "_", sep="")

## pdf(file="treeWithBars-025.pdf", height=100, width=10)
## par(mai=c(0.5,0,2,0), xpd=T)
## plot(as(spGrp, "phylo"), cex=.5, show.tip.label=T, no.margin=F, label.offset=0)
## barMonophyletic(grpLbl, as(spGrp, "phylo"), extra.space=0.01, cex.plot=.5, cex.text=.5,
##                 bar.at.tips=TRUE, include.tip.label=TRUE)
## add.scale.bar()
## dev.off()
