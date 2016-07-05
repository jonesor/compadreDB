\name{kentropy}
\alias{kentropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Keyfitz' life table entropy
}
\description{
This function calculates Keyfitz' life table entropy from an lx (surivorship) vector with even intervals.
}
\usage{
kentropy(lx, trapeze = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lx}{
A numerical vector of lx (survivorship). This is assumed to be with a constant interval (e.g. 1yr).
}
  \item{trapeze}{
A logical argument indicating whether the trapezoidal approximation should be used for approximating the definite integral.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an estimate of Keyfitz' life table entropy based on an lx (survivorship) vector.
}
\references{

}
\author{
Owen R. Jones <jones@biology.sdu.dk>

Roberto Salguero-Gomez <salguero@sheffield.ac.uk>

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#Survivorship (lx) with constant mortality (should have Keyfitz entropy of 1).

x <- 0:100
qx <- .4
px <- 1-qx
lx <- px^x

plot(x,lx,type="l",ylim=c(0,1),col="red")
kentropy(lx)
kentropy(lx,trapeze=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
