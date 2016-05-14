\name{makeLifeTable}
\alias{makeLifeTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeLifeTable(matU, matF = NULL, matC = NULL, startLife = 1, nSteps = 1000)
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
function (matU, matF = NULL, matC = NULL, startLife = 1, nSteps = 1000) 
{
    matDim = ncol(matU)
    matUtemp = matU
    survivorship = array(NA, dim = c(nSteps, matDim))
    for (o in 1:nSteps) {
        survivorship[o, ] = colSums(matUtemp \%*\% matU)
        matUtemp = matUtemp \%*\% matU
    }
    lx = survivorship[, startLife]
    lx = c(1, lx[1:(length(lx) - 1)])
    out = data.frame(x = 0:(length(lx) - 1), lx = lx)
    if (!missing(matF)) {
        if (sum(matF, na.rm = T) == 0) {
            warning("matF contains only 0 values")
        }
        ageFertility = array(0, dim = c(nSteps, matDim))
        fertMatrix = array(0, dim = c(nSteps, matDim))
        matUtemp2 = matU
        e = matrix(rep(1, matDim))
        for (q in 1:nSteps) {
            fertMatrix = matF \%*\% matUtemp2 * (as.numeric((ginv(diag(t(e) \%*\% 
                matUtemp2)))))
            ageFertility[q, ] = colSums(fertMatrix)
            matUtemp2 = matUtemp2 \%*\% matU
        }
        mx = ageFertility[, startLife]
        mx = c(0, mx[1:(length(mx) - 1)])
        out$mx = mx
    }
    if (!missing(matC)) {
        if (sum(matC, na.rm = T) == 0) {
            warning("matC contains only 0 values")
        }
        ageClonality = array(0, dim = c(nSteps, matDim))
        clonMatrix = array(0, dim = c(nSteps, matDim))
        matUtemp2 = matU
        e = matrix(rep(1, matDim))
        for (q in 1:nSteps) {
            clonMatrix = matC \%*\% matUtemp2 * (as.numeric((ginv(diag(t(e) \%*\% 
                matUtemp2)))))
            ageClonality[q, ] = colSums(clonMatrix)
            matUtemp2 = matUtemp2 \%*\% matU
        }
        cx = ageClonality[, startLife]
        cx = c(0, cx[1:(length(cx) - 1)])
        out$cx = cx
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
