\name{subsetDB}
\alias{subsetDB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subsetDB(db = comadre, sub = 1:100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{
%%     ~~Describe \code{db} here~~
}
  \item{sub}{
%%     ~~Describe \code{sub} here~~
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
function (db = comadre, sub = 1:100) 
{
    subsetID <- sub
    ssdb <- db
    ssdb$metadata <- ssdb$metadata[subsetID, ]
    ssdb$mat <- ssdb$mat[subsetID]
    ssdb$matrixClass <- ssdb$matrixClass[subsetID]
    if ("version" \%in\% names(ssdb)) {
        ssdb$version$Version <- paste(ssdb$version$Version, " - subset created on ", 
            format(Sys.time(), "\%b_\%d_\%Y"), sep = "")
        ssdb$version$DateCreated <- paste(ssdb$version$DateCreated, 
            " - subset created on ", format(Sys.time(), "\%b_\%d_\%Y"), 
            sep = "")
        ssdb$version$NumberAcceptedSpecies <- length(unique(ssdb$metadata$SpeciesAccepted))
        ssdb$version$NumberStudies <- length(unique(ssdb$metadata$SpeciesAuthor))
        ssdb$version$NumberMatrices <- length(ssdb$mat)
    }
    return(ssdb)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
