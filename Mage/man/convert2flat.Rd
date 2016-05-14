\name{convert2flat}
\alias{convert2flat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
convert2flat(db = compadre, Aonly = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{
%%     ~~Describe \code{db} here~~
}
  \item{Aonly}{
%%     ~~Describe \code{Aonly} here~~
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
function (db = compadre, Aonly = TRUE) 
{
    db$metadata$Amatrix <- NULL
    for (i in 1:nrow(db$metadata)) {
        db$metadata$classnames[i] <- paste(db$matrixClass[[i]]$MatrixClassAuthor, 
            collapse = "| ")
        db$metadata$matrixA[i] <- paste("[", paste(t(db$mat[[i]]$matA), 
            collapse = " "), "]", sep = "")
        db$metadata$matrixU[i] <- paste("[", paste(t(db$mat[[i]]$matU), 
            collapse = " "), "]", sep = "")
        db$metadata$matrixF[i] <- paste("[", paste(t(db$mat[[i]]$matF), 
            collapse = " "), "]", sep = "")
        db$metadata$matrixC[i] <- paste("[", paste(t(db$mat[[i]]$matC), 
            collapse = " "), "]", sep = "")
    }
    return(db$metadata)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
