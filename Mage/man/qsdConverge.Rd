\name{qsdConverge}
\alias{qsdConverge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculates the age at which a defined quasi-stationary stage distribution is reached.
}
\description{
Calculates the age at which a defined quasi-stationary stage distribution is reached.
}

\usage{
qsdConverge(matU, conv = 0.05, startLife = 1, nSteps = 1000)
}

\arguments{
  \item{matU}{
The survival/growth matrix (a subset of the A matrix), See Caswell 2001, or other references below.
}
  \item{conv}{
Departure from convergence. e.g. this value should be 0.05 if you want to get the time step when the stage distribution is within 5 percent of the QSD.
}
  \item{startLife}{
At which stage does the life start in the matrix. Some matrices contain seed banks, or similar.
}
  \item{nSteps}{
For how many time steps should the calculations run for. In most cases the default of 1000 is sufficient.
}
}
\details{
Matrix Population Models are typically parameterised with a stasis loop in the oldest/largest/most-developed stage (e.g. adult survival), which means that flat mortality and fertility plateaus may emerge as a mathematical artefact when examining age-specific patterns (see, Caswell 2001). The Quasi-Stationary stage Distribution (QSD) can be used to avoid this problem. The QSD can be defined as a stage distribution that is reached some time before the ultimate stable stage distribution (the normalised righted eigenvector of the MPM). Thus one can ask what is the age at which the cohort approximates, say,  95\% of its stable stage distribution, and only use age-based information from before this point. This is discussed and explained in detail in the supplementary information of Jones et al. (2014).
}

\value{
Returns a numeric value indicating the time step at which quasi-stationary stage distribution has been reached.
}

\references{
Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968

Jones, O.R. et al. (2014) Diversity of ageing across the tree of life. Nature, 505(7482), 169–173
}

\author{
Owen R. Jones <jones@biology.sdu.dk>

Roberto Salguero-Gomez <r.salguero@sheffield.ac.uk>

Hal Caswell

}

\note{}

\seealso{}

\examples{

U1 <- matrix(c(0.000, 0.691, 0.000, 0.000, 0.000, 0.880, 0.000, 0.000, 0.790),nrow=3)
qsdConverge(U1, 0.05)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

