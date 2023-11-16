\name{bdr.gs.mc.gen}
\alias{bdr.gs.mc.gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the boundaries for the group sequential tests
}
\description{
Returns the boundaries for any group sequential test of the null vs. 2-sided alternative whose boundaries take the form of a single constant times a known weight vector, which is w.vec.  These include Pocock (w.vec=(1,1,..)), O'Brien-Fleming (w.vec=(sqrt(n.stg/1), sqrt(n.stg/2), ..., 1)), etc.  It does this by returning quantiles of the sample paths of the (null) test statistic paths in mc.paths.   

}
\usage{
bdr.gs.mc.gen(alpha = 0.05, mc.paths, w.vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
desired rejection probability of the test; default is 0.05
}
  \item{mc.paths}{
matrix of sample paths, each row being a sample path, no. of columns is number of stages
}
  \item{w.vec}{
weight vector corresponding to desired test
}
}

\value{
Returns a list:
  \item{cons}{the constant in the boundary vector cons*w.vec}
 \item{bndry.vec}{the boundary vector cons*w.vec}

}

\author{
Jay Bartroff
}
\keyword{internal}
