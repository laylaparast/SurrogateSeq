\name{bdr.gs.mc.fut}
\alias{bdr.gs.mc.fut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the boundaries for the group sequential tests with futility stopping
}
\description{
Returns the boundaries for any group sequential test of the null vs. 2-sided alternative whose boundaries take the form of a single constant times a known weight vector, which is w.vec; allows for futility stopping.  These include Pocock (w.vec=(1,1,..)), O'Brien-Fleming (w.vec=(sqrt(n.stg/1), sqrt(n.stg/2), ..., 1)), etc.  It does this by returning quantiles of the sample paths of the (null) test statistic paths in mc.paths.   
}
\usage{
bdr.gs.mc.fut(c1 = NULL, c2 = NULL, pp = 0.4, n.stg, j.star = 1, alpha = 0.05, alpha0 = (j.star/n.stg) * alpha, mc.paths, inf.fraction = (1:n.stg)/n.stg, N.iter.max = 100, alpha.tol = 0.02 * alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c1}{
c1 and c2 are the constants determining the outer boundary b[j] = c1 * (j/J)^{pp-1/2} and futility boundaries a[j] = (c1+c2) * (j/J)^{1/2} - c2 * (j/J)^{pp-1/2} for j >= j.star, where J is the max no of stages (AKA n.stg). If c1 is null, it is found as the upper alpha0 quantile of the max over the first j.star stages.
}
  \item{c2}{
see description in c1
}
  \item{pp}{
power parameter for Wang-Tsiatis boundaries; default is 0.4
}
  \item{n.stg}{
maximum number of analyses
}
  \item{j.star}{
earliest stage at which futility stopping is allowed. Should be <= n.stg-1 (there is already "futility stopping" at the n.stg-th stage anyway). Default is 1.

}
  \item{alpha}{
desired rejection probability of the test; default is 0.05
}
  \item{alpha0}{
the part of alpha that c1 is chosen to spend in first j.star stages; default is (j.star/n.stg) * alpha
}
  \item{mc.paths}{
matrix of sample paths, each row being a sample path, no. of columns is number of stages
}
  \item{inf.fraction}{
information fraction vector of the same length as n.stg which reflects the fraction of information at each analysis, should be positive, non-decreasing, and the last entry should be 1; default is (1:n.stg)/n.stg, user may want to specify a different vector for unequal time points
}
  \item{N.iter.max}{
max no. of iterations for finding c2
}
  \item{alpha.tol}{
the tolerance for stopping search for c2
}
}

\value{
Returns a list:
  \item{a}{the futility boundary vectors}
    \item{b}{the null-rejection boundary vectors}
     \item{prej}{prob. of rejecting the null (at any stage)}
      \item{EM}{ expected stopping stage number}
       \item{se.M}{standard error of stopping time}
        \item{c1}{constants used in boundaries a, b}
        \item{c2}{constants used in boundaries a, b}
 }
\author{
Jay Bartroff
}
\keyword{internal}
