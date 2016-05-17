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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
