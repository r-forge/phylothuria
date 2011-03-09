\name{findGroups}
\alias{findGroups}
\title{Find evolutionary significant groups on a tree based on edge-length
  distances}
\description{
  Find evolutionary significant groups on a tree based on edge-length
  distances. All terminal tips that are less than a given threshold will
  be grouped together.
}
\usage{
  findGroups(tr, threshold=0.025)
}
\arguments{
  \item{tr}{a \code{phylo4}  object}
  \item{threshold}{the maximum distance between the terminal tip and the
  most recent common ancestor to the group. Twice the threshold gives
  the maximum distance between two tips in a group.}
}
\value{
  Returns a \code{phylo4d} object with a variable \sQuote{Groups} that
  contains the group number to which the tip belong. The group numbers
  are numbered from 1 to the total number of groups.
}
\details{
  TODO: needs to be modified to accept phylo4d objects
  TODO: needs to fix bug when only 2 sequences in a group
}
\seealso{
  distFromTip
}
\examples{
  ## Not run
  # findGroups(tr, threshold=0.024)
  ## End not run
}
\author{Francois Michonneau}
\keyword{cluster}
