\name{op.char.gs.fut}
\alias{op.char.gs.fut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the operating characteristics on the group sequential test with futility stopping statistics in paths
}
\description{
Compute the operating characteristics on the group sequential test with futility stopping statistics in paths: The expected stopping stage no., plus the probability of rejecting the null in favor of the 2-sided alternative. This is for a general GS test which uses the boundaries in bndry.vec.
}
\usage{
op.char.gs.fut(b.vec, a.vec, paths)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b.vec}{
"null-rejection" boundaries,  should be >= 0, and a.vec[n.stg] = b.vec[n.stg]. 
}
  \item{a.vec}{
futility boundaries, should be >= 0, and a.vec[n.stg] = b.vec[n.stg];a.vec[j]=0 means no futility stopping at stage j.
}
  \item{paths}{
matrix of test statistic sample paths, each row being a sample path, no. of columns is max number 
}
}

\value{
Returns a list:
          \item{EM}{ expected stopping stage number}
       \item{se.M}{standard error of stopping time}
        \item{prej}{prob. of rejecting the null (at any stage)}

 }

\keyword{internal}