\name{qsdConverge}
\alias{qsdConverge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qsdConverge(matU, conv = 0.05, startLife = 1, nSteps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matU}{
%%     ~~Describe \code{matU} here~~
}
  \item{conv}{
%%     ~~Describe \code{conv} here~~
}
  \item{startLife}{
%%     ~~Describe \code{startLife} here~~
}
  \item{nSteps}{
%%     ~~Describe \code{nSteps} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (matU, conv = 0.05, startLife = 1, nSteps = 1000) 
{
    uDim = dim(matU)
    eig = eigen.analysis(matU)
    qsd = eig$stable.stage
    qsd = as.numeric(t(matrix(qsd/sum(qsd))))
    nzero = rep(0, uDim[1])
    nzero[startLife] = 1
    n = nzero
    dist = p = NULL
    survMatrix1 <- matU
    for (j in 1:nSteps) {
        p = n/sum(n)
        dist[j] = 0.5 * (sum(abs(p - qsd)))
        n = survMatrix1 \%*\% n
    }
    convage = min(which(dist < conv))
    return(convage)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
