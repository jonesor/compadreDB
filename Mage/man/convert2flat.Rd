\name{convert2flat}
\alias{convert2flat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to convert from the list structured database object to a flat sheet.
}
\description{
The function converts from the list structured database object to a flat sheet by converting each matrix to a string.
}
\usage{
convert2flat(db, Aonly = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{
The COMPADRE or COMADRE matrix database object
}
  \item{Aonly}{
A logical value (TRUE/FALSE) indicatting whether ONLY the A matrix be included in the output.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A `data.frame` with the same columns as are present in the metadata part of the COMPADRE/COMADRE object, followed by the string-form matrix stage information and the matrices themselves.
}
\references{

}
\author{
Owen R. Jones <jones@biology.sdu.dk>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\dontrun{
newDB<-convert2flat(compadre,Aonly=FALSE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
