\name{vitalRates}
\alias{vitalRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to derive vital rates from the matrix population model.
}
\description{
A function to derive vital rates from the matrix population model for the separate demographic processes.
}
\usage{
vitalRates(matU, matF, matC, splitStages = FALSE, weighted = FALSE)
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
  \item{splitStages}{
Splits vital rates according to some pre-determined criteria (below).
}
  \item{weighted}{
Allows to weight mean vital rates according to various criteria (below).
}
}
\details{
This function decomposes vital rates of survival, progression, retrogression, sexual reproduction and clonal reproduction according to various ways of weighted means and organization of stages along the life cycle represented in the matrix population model.
}
\value{
- 'Weighted': This argument allows to weight mean values of vital rates (survival 'surv', progression 'prog', retrogression 'retr', sexual reproduction 'fec' and clonal reproduction 'clo') with an equal contribution for all stages (default), by the stable st/age distribution ('SSD'), or by a given population vector chosen by the user, so long as it is congruent with the dimensions of the chosen 'matU', 'matF', and 'matC'.

- 'splitStages': This argument allows to split the values of vital rates according to recognizable stages in the matrix. When 'all', all vital rates are averaged all existing stages, if 'ontogeny', they are averaged as juveniles ('Juv') and adults ('Adu'), and if by 'MatrixClassOrganized', it takes a vector with the pre-established stages of 'compadre$matrixClass[[i]]$MatrixClassOrganized' or 'compadre$matrixClass[[i]]$MatrixClassOrganized', where 'i' is the index of the chosen study in 'COMPADRE' or 'COMADRE'.

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

#Vital rate outputs without weights:
vitalRates(matU, matF, matC, splitStages = 'all', weighted = FALSE)
vitalRates(matU, matF, matC, splitStages = 'ontogeny', weighted = FALSE)
vitalRates(matU, matF, matC, splitStages = c('prop', 'active', 'active', 'active'), weighted = FALSE)



#Vital rate outputs weighted by the stable stage distribution of 'matA':
vitalRates(matU, matF, matC, splitStages = 'all', weighted = 'SSD')
vitalRates(matU, matF, matC, splitStages = 'ontogeny', weighted = 'SSD')
vitalRates(matU, matF, matC, splitStages = c('prop', 'active', 'active', 'active'), weighted = 'SSD')

#Vital rate outputs weighted by a chosen population vector of initial conditions:
initialConditions <- c(100, 10, 0, 1)

vitalRates(matU, matF, matC, splitStages = 'all', weighted = initialConditions)
vitalRates(matU, matF, matC, splitStages = 'ontogeny', weighted = initialConditions)
vitalRates(matU, matF, matC, splitStages = c('prop', 'active', 'active', 'active'), weighted = initialConditions)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
