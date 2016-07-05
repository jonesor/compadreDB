\name{makeLifeTable}
\alias{makeLifeTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to produce a life table from a matrix population model
}
\description{
This function uses age-from-stage approximation methods to calculate a life table.
}
\usage{
makeLifeTable(matU, matF = NULL, matC = NULL, startLife = 1, nSteps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matU}{
The U matrix - growth and survival transitions.
}
  \item{matF}{
The F matrix - sexual reproduction
}
  \item{matC}{
The C matrix - clonal reproduction
}
  \item{startLife}{
The stage of the matrix where life is considered to begin. This is usually 1 (the first stage), but in species with a seed stage (or other propagule or dormant stage), may be 2 or more.
}
  \item{nSteps}{
For how many time steps should the life table be constructed.
}
}
\details{
A detailed description of these methods can be found in the supplementary information of Jones et al. (2014) and in section 5.3 of Caswell (2001).
}
\value{
A `data.frame` with between 2 and 4 columns.
If just `matU` is provided, the `data.frame` has columns `x` (age), `lx` (survivorship). 
If `matF` is provided (in addition to `matU`), `mx` (age-specific sexual reproduction) is included in the output `data.frame`. Likewise, if `matC` is provided, an addition column, `cx` (age-specific clonal reproduction) is included.
}
\references{

Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968

Jones, O.R. et al. (2014) Diversity of ageing across the tree of life. Nature, 505(7482), 169–173

}
\author{
Owen R. Jones <jones@biology.sdu.dk>

Roberto Salguero-Gomez <r.salguero@sheffield.ac.uk>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(popbio)
data(whale)
matU <- splitA(whale)$T
matF <- splitA(whale)$F
makeLifeTable(matU, matF, startLife = 1, nSteps = 100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
