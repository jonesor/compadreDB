\name{dbCompare}
\alias{dbCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dbCompare(db1, db2, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db1}{
%%     ~~Describe \code{db1} here~~
}
  \item{db2}{
%%     ~~Describe \code{db2} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (db1, db2, verbose = FALSE) 
{
    cat("Quick Summary\n")
    uniqueSource1 <- unique(paste(db1$metadata$Authors, " (", 
        db1$metadata$YearPublication, ") ", db1$metadata$Journal, 
        sep = ""))
    db1$metadata$binomial <- paste(db1$metadata$GenusAccepted, 
        db1$metadata$SpeciesEpithetAccepted, sep = " ")
    cat(paste("File-1 contains the demographic and associated data from ", 
        length(uniqueSource1), " source papers, corresponding to ", 
        length(unique(db1$metadata$binomial)), " accepted species, and ", 
        nrow(db1$metadata), " matrices.\n\n", sep = ""))
    uniqueSource2 <- unique(paste(db2$metadata$Authors, " (", 
        db2$metadata$YearPublication, ") ", db2$metadata$Journal, 
        sep = ""))
    db2$metadata$binomial <- paste(db2$metadata$GenusAccepted, 
        db2$metadata$SpeciesEpithetAccepted, sep = " ")
    cat(paste("File-2 contains the demographic and associated data for ", 
        length(uniqueSource2), " source papers, corresponding to ", 
        length(unique(db2$metadata$binomial)), " accepted species, and ", 
        nrow(db2$metadata), " matrices.\n\n", sep = ""))
    if (verbose == TRUE) {
        cat("Detailed summary\n")
        sp1 <- unique(db1$metadata$binomial)
        sp2 <- unique(db2$metadata$binomial)
        cat("Number of accepted species in File 1, based on latin binomial\n")
        print(length(sp1))
        cat("Number of accepted species in File 2, based on latin binomial\n")
        print(length(sp2))
        cat("Accepted species in File 1 that are not in File 2 (based on latin binomial)\n")
        print(sp1[which(!sp1 \%in\% sp2)])
        cat("Accepted species in File 2 that are not in File 1 (based on latin binomial)\n")
        print(sp2[which(!sp2 \%in\% sp1)])
        asp1 <- unique(db1$metadata$SpeciesAuthor)
        asp2 <- unique(db2$metadata$SpeciesAuthor)
        cat("Number of unique study-species combinations in File 1\n")
        print(length(asp1))
        cat("Number of unique study-species combinations in File 2\n")
        print(length(asp2))
        cat("\n\nSource papers in File 2 that are not in File 1\n")
        print(sort(uniqueSource2[which(!uniqueSource2 \%in\% uniqueSource1)]))
        cat("\n\nSource papers in File 1 that are not in File 2\n")
        print(sort(uniqueSource1[which(!uniqueSource1 \%in\% uniqueSource2)]))
        cat("See the User Guide for definitiions\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
