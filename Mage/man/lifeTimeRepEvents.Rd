\name{lifeTimeRepEvents}
\alias{lifeTimeRepEvents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate the timing of lifetime reproductive events.
}
\description{
A function to calculate the timing of lifetime reproductive events such as the probability of achieving maturity, age at first reproduction, mean life expectancy conditional on maturity, and life expectancy for mature individuals
}
\usage{
lifeTimeRepEvents(matU, matF, matC, startLife = 1)
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
  \item{startLife}{
The first stage at which the author consider the beginning of life.
}
}
\details{
This function applies Markov chain approaches to decompose various moments of along the age-based reproduction of individuals in a matrix population model.
}
\value{
This function applies Markov chain approaches to decompose various moments of along the age-based reproduction of individuals in a matrix population model. When both a 'matF' and a 'matC' are provided, the following outputs are calculated for both independently, and these are differentiated with the suffix "Fec" or "Clo", respectively:

- 'p': probability of achiving maturity, sexual or clonal.

- 'La': mean age at maturity (in the same units as the matrix population model).

- 'meanLifeExpectancy': mean life expectancy conditional on entering the life cycle in the first reproductive stage

- 'remainingMatureLifeExpectancy': Life expectancy from mean maturity. This is mean life expectancy - mean age at maturity ('La' above). This value can be negative because both mean life expectancy and mean age at maturity are means of their respective distributions.
}
\references{
Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968
}
\author{
Roberto Salguero-Gomez <rob.salguero@zoo.ox.ac.uk>

Hal Caswell <hcaswell@whoi.edu>

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
matU <- matrix (c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
matF <- matrix (c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
matC <- matrix (c(0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)

lifeTimeRepEvents(matU, matF, matU, startLife = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
