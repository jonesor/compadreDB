\name{lifeTimeRepEvents}
\alias{lifeTimeRepEvents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lifeTimeRepEvents(matU, matF, startLife = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matU}{
%%     ~~Describe \code{matU} here~~
}
  \item{matF}{
%%     ~~Describe \code{matF} here~~
}
  \item{startLife}{
%%     ~~Describe \code{startLife} here~~
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
function (matU, matF, startLife = 1) 
{
    uDim = dim(matU)[1]
    surv = colSums(matU)
    repLifeStages = colSums(matF)
    repLifeStages[which(repLifeStages > 0)] = 1
    if (missing(matF) | missing(matU)) {
        stop("matU or matF missing")
    }
    if (sum(matF, na.rm = T) == 0) {
        stop("matF contains only 0 values")
    }
    Uprime = matU
    Uprime[, which(repLifeStages == 1)] = 0
    Mprime = matrix(0, 2, uDim)
    for (p in 1:uDim[1]) {
        if (repLifeStages[p] == 1) 
            Mprime[2, p] = 1
        else Mprime[1, p] = 1 - surv[p]
    }
    Bprime = Mprime \%*\% (ginv(diag(uDim) - Uprime))
    pRep = Bprime[2, startLife]
    out = data.frame(pRep = pRep)
    D = diag(c(Bprime[2, ]))
    Uprimecond = D \%*\% Uprime \%*\% ginv(D)
    expTimeReprod = colSums(ginv(diag(uDim) - Uprimecond))
    La = expTimeReprod[startLife]
    out$La = La
    firstRepLifeStage = min(which(repLifeStages == 1))
    N = solve(diag(uDim[1]) - matU)
    meanRepLifeExpectancy = colSums(N)[firstRepLifeStage]
    out$meanRepLifeExpectancy = meanRepLifeExpectancy
    remainingMatureLifeExpectancy = colSums(N)[startLife] - La
    out$remainingMatureLifeExpectancy = remainingMatureLifeExpectancy
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
