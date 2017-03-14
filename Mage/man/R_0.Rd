\name{R_0}
\alias{R_0}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate net reproductive value
}
\description{
This function calculates net reproductive value from a matU (survival-dependent processes) and a matF (sexual reproduction) and/or a matC (clonal reproduction.
}
\usage{
R_0(matU, matF, matC = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matU}{
A matrix containing only survival-dependent processes ( growth, stasis, shrinkage).
}
  \item{matF}{
A matrix containing only sexual reproduction, with zeros elsewhere.
}
  \item{matC}{
A matrix containing only clonal reproduction, with zeros elsewhere.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the net reproductive value of the matrix. When both 'matF' and 'matC' are provided, it outputs the net reprodutive value for sexual reproduction only, for clonal reproduction only, and for both types of reproduction together.
}
\references{
Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968
}
\author{
Roberto Salguero-Gómez <rob.salguero@zoo.ox.ac.uk>

Hal Caswell <h.caswell@uva.nl>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

matU <- matrix (c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
matF <- matrix (c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
matC <- matrix (c(0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)

R_0(matU, matF, matU)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
