\name{distFromTip}
\alias{distFromTip}
\title{Distance from Tip}
\description{
  Returns the maximum edge-length distance from any node to the tips.
}
\usage{
  distFromTip(tr, node)
}
\arguments{
  \item{tr}{a \code{phylo4} or \code{phylo4d} object}
  \item{node}{a node label or node number.}
}
\value{
  Returns the edge-length between the node and the tip.
}
\details{
  Designed mostly for ultrametric trees.
  TODO: test more thouroughly and see if it can be optimized for
  ultrametric trees (i.e., we probably don't need to use max() in that case).
}
\seealso{
  getNode and other treewalk functions
}
\examples{
  data(geospiza)
  distFromTip(geospiza, 20)
}
\author{Francois Michonneau}
\keyword{distance}
