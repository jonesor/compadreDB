\name{longevity}
\alias{longevity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate measures of longevity.
}
\description{
A function to calculate the mean life expectancy and maximum longevity of individuals in a matrix population model
}
\usage{
longeivity(matU, startLife, initPop, run)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matU}{
A matrix containing only survival-dependent processes (growth, stasis, shrinkage).
}
  \item{startLife}{
The first stage at which the author consider the beginning of life.
}
  \item{initPop}{
Initial population size.
}
  \item{run}{
Number of iterations that the maximum longevity will be calculated for.
}
}
\details{
This function applies Markov chain approaches to obtain mean life expectancy, and a loop to calculate maximum longevity.
}
\value{
This function applies Markov chain approaches to obtain mean life expectancy, and a loop to calculate maximum longevity. Outputs are:

- 'eta': mean life expectancy conditional on entering the life cycle on life stage described by 'startLife'.

- 'Max': maximum longevity observed by iterating a population vector with 'initPop' individuals in stage 'startLife' up to 'run' times.
}
\references{
Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968

Doak & Morris (2002) BOOK INSERT CITATION
}
\author{
Roberto Salguero-Gomez <rob.salguero@zoo.ox.ac.uk>

Hal Caswell <hcaswell@whoi.edu>

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

longevity(matU, startLife = 1, initPop = 100, run = 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
