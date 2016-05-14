\name{vitalRatePerturbation}
\alias{vitalRatePerturbation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vitalRatePerturbation(matU, matF, matC = NULL, pert = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matU}{
%%     ~~Describe \code{matU} here~~
}
  \item{matF}{
%%     ~~Describe \code{matF} here~~
}
  \item{matC}{
%%     ~~Describe \code{matC} here~~
}
  \item{pert}{
%%     ~~Describe \code{pert} here~~
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
function (matU, matF, matC = NULL, pert = 0.001) 
{
    matA = matU + matF + matC
    aDim = dim(matA)[1]
    fakeA = matA
    sensA = elasA = matrix(NA, aDim, aDim)
    lambda = Re(eigen(matA)$values[1])
    propU = matU/matA
    propU[is.nan(propU)] = 0
    propProg = propRetrog = propU
    propProg[upper.tri(propU, diag = T)] = 0
    propRetrog[lower.tri(propU, diag = T)] = 0
    propStasis = matrix(Diagonal(aDim) * diag(propU), aDim, aDim)
    propF = matF/matA
    propF[is.nan(propF)] = 0
    propC = matC/matA
    propC[is.nan(propC)] = 0
    sensA = eigen.analysis(matA, zero = F)$sensitivities
    uIndep = matrix(NA, aDim, aDim)
    u = colSums(matU)
    for (j in which(u > 0)) uIndep[, j] = matA[, j]/u[j]
    sensSigmaA = uIndep * sensA
    uPrime = u
    elasSigmaA = t(t(sensSigmaA) * uPrime)/lambda
    elasA = sensA * matA/lambda
    uDistrib = matrix(0, ncol = aDim, nrow = aDim)
    for (j in which(u > 0)) uDistrib[, j] = matU[, j]/u[j]
    f = colSums(matF)
    fDistrib = matrix(0, ncol = aDim, nrow = aDim)
    for (j in which(f > 0)) fDistrib[, j] = matF[, j]/f[j]
    c = colSums(matC)
    cDistrib = matrix(0, ncol = aDim, nrow = aDim)
    for (j in which(c > 0)) cDistrib[, j] = matC[, j]/c[j]
    SuDistrib = sensA * uDistrib
    SfDistrib = sensA * fDistrib
    ScDistrib = sensA * cDistrib
    out = data.frame(SSurvival = NA, SGrowth = NA, SShrinkage = NA, 
        SReproduction = NA, SClonality = NA, ESurvival = NA, 
        EGrowth = NA, EShrinkage = NA, EReproduction = NA, EClonality = NA)
    out$SSurvival = sum(sensSigmaA, na.rm = T)
    out$SGrowth = sum(sensA * uDistrib * propProg, na.rm = T)
    out$SShrinkage = sum(sensA * uDistrib * propRetrog, na.rm = T)
    out$SReproduction = sum(sensA * fDistrib * propF, na.rm = T)
    out$SClonality = sum(sensA * cDistrib * propC, na.rm = T)
    EuDistrib = sensA * uDistrib * matrix(u, nrow = aDim, ncol = aDim, 
        byrow = T)/lambda
    EfDistrib = sensA * fDistrib * matrix(f, nrow = aDim, ncol = aDim, 
        byrow = T)/lambda
    EcDistrib = sensA * cDistrib * matrix(c, nrow = aDim, ncol = aDim, 
        byrow = T)/lambda
    out$ESurvival = sum(elasSigmaA, na.rm = T)
    out$EGrowth = sum(EuDistrib * propProg, na.rm = T)
    out$EShrinkage = sum(EuDistrib * propRetrog, na.rm = T)
    out$EReproduction = sum(EfDistrib * propF, na.rm = T)
    out$EClonality = sum(EcDistrib * propC, na.rm = T)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
