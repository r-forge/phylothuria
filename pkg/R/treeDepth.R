treeDepth <- function(tr) {
  if (is.null(tr$root.edge))
    max(branching.times(tr))
  else
    max(branching.times(tr) + tr$root.edge)
}
