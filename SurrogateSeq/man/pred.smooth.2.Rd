\name{pred.smooth.2}
\alias{pred.smooth.2}
\title{
Calculates the conditional mean function
}
\description{
Helper function; calculates the condition mean of Y given S, based on Study A data
}
\usage{
pred.smooth.2(kernel.use,kernel.apply, bw,outcome)
}
\arguments{
  \item{kernel.use}{
surrogate values in the control group in Study A
}
  \item{kernel.apply}{
%%     ~~Describe \code{zi.one} here~~
surrogate values in Study B
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
bandwidth
}
  \item{outcome}{
%%     ~~Describe \code{y1} here~~
outcome in the control group in Study A
}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
expected outcome for each surrogate value 
}
\author{
%%  ~~who you are~~
Layla Parast
}

\keyword{internal}
