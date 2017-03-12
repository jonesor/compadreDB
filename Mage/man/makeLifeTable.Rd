\name{makeLifeTable}
\alias{makeLifeTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to produce a life table from a matrix population model
}
\description{
This function uses age-from-stage decompositions to calculate a life table from a matrix population model
}
\usage{
makeLifeTable(matU, matF = NULL, matC = NULL, startLife = 1, nSteps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matU}{
The U matrix - survival-dependent transitions (e.g. changes in size/ageing/development)
}
  \item{matF}{
The F matrix - sexual reproduction
}
  \item{matC}{
The C matrix - clonal reproduction
}
  \item{startLife}{
The stage of of the life cycle that structures the matrix population model where life is considered to begin. This is usually 1 (i.e. the first stage), but in species with a permanent propagule bank (e.g. seedbank, or similar), users may wish to use the first non-permanent propagule bank stage along the life cycle.
}
  \item{nSteps}{
Number of time steps for which the life table be constructed. This is on the same units as the matrix population model - see MatrixPeriodicity in metadata of COMPADRE/COMADRE.
}
}
\details{
A detailed description of these methods can be found in section 5.3 of Caswell (2001) and the supplementary information of Jones et al. (2014).
}
\value{
A `data.frame` with between 2 and 4 columns.
If just `matU` is provided, the `data.frame` has columns `x` (age), `lx` (survivorship). 
If `matF` is provided (in addition to `matU`), `mx` (age-specific sexual reproduction) is included in the output `data.frame`. Likewise, if `matC` is provided, an addition column, `cx` (age-specific clonal reproduction) is included.
}
\references{

Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968

Caswell, H. (2006) Applications of Markov chains in demography. pp. 319-334 in A.N. Langville and W.J. Stewart (editors) MAM2006: Markov Anniversary Meeting. Boson Books, Raleigh, North Caroline, USA

Jones, O.R. et al. (2014) Diversity of ageing across the tree of life. Nature, 505(7482), 169–173

}
\author{
Roberto Salguero-Gómez <rob.salguero@sheffield.ac.uk>

Hal Caswell <h.caswell@uva.nl>

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

matU <- matrix(c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
matF <- matrix(c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.), nrow = 4, byrow = T)
matC <- matrix(c(0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)

makeLifeTable(matU, matF, matC, startLife = 1, nSteps = 100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
