\name{gs.boundaries}
\alias{gs.boundaries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes group sequential boundaries }
\description{
Computes group sequential (and naive) boundaries for the nonparametric test for a treatment effect on the primary outcome using surrogate marker information. The boundaries and test statistic borrow information from a prior study (Study A) about the relationship between the surrogate and the primary outcome to test for a treatment effect in the current study (Study B).}
\usage{
gs.boundaries(szerop, sonep, yzerop, n.stg, B.norm = 1e+06, alpha = 0.05, pp = 0.4,
 inf.fraction = (1:n.stg)/n.stg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{szerop}{
surrogate marker in the control group in Study A
}
  \item{sonep}{
surrogate marker in the treated group in Study A
}
  \item{yzerop}{
primary outcome in the control group in Study A
}
  \item{n.stg}{
maximum number of analyses
}
  \item{B.norm}{
number of multivariate normal vectors to use in simulation for boundaries; default is 1e+06
}
  \item{alpha}{
desired rejection probability of the test; default is 0.05
}
\item{pp}{
	power parameter for Wang-Tsiatis boundaries; default is 0.4
}
\item{inf.fraction}{
	information fraction vector of the same length as n.stg which reflects the fraction of information at each analysis, should be positive, non-decreasing, and the last entry should be 1; default is (1:n.stg)/n.stg, user may want to specify a different vector for unequal time points
}

}

\value{
Returns a list of boundaries:
 \item{Naive}{Naive boundaries }
 \item{Bonf}{Bonferroni boundaries}
 \item{Pocock}{Pocock boundaries}
 \item{OBrien_Fleming}{O'Brien-Fleming boundaries}
 \item{Wang_Tsiatis}{Wang-Tsiatis boundaries}

}
\references{
Parast and Bartroff (2023+). Group Sequential Testing of a Treatment Effect Using a Surrogate Marker. Under Review.
}
\author{
Layla Parast and Jay Bartroff
}
\examples{
data(example.data)
data(StudyA.aids)
data(StudyB.aids)
s0.studya = StudyA.aids$s0
s1.studya = StudyA.aids$s1

bound = gs.boundaries(szerop = s0.studya, sonep = s1.studya, 
yzerop=StudyA.aids$y0, n.stg=4, B.norm=1e6, alpha=0.05)

bound
}
\keyword{nonparametric}
\keyword{test}
\keyword{boundaries}
