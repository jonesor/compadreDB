\name{plotLifeCycle}
\alias{plotLifeCycle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot the life cycle diagram based on the matrix model
}
\description{
This function plots the life cycle diagram illustrated by a matrix model. It processes the matrix model and passes the information to the graphViz functionality of DiagrammeR. See http://rich-iannone.github.io/DiagrammeR/.
}
\usage{
plotLifeCycle(A, title = "my life cycle", shape = "egg", fontsize = 10, 
	nodefontsize = 12, edgecol = "grey")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
The A matrix of a matrix population model
}
  \item{title}{
A title for the plot
}
  \item{shape}{
The shape to be used for the stages of the diagram. Any node shape accepted by graphViz is acceptable.
}
  \item{fontsize}{
Size of the font used in the diagram.
}
  \item{nodefontsize}{
Size of the font used in the node part of the diagram.
}
  \item{edgecol}{
Colour of the arrows in the diagram.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Produces a plot consisting of a life cycle diagram, based on the A matrix.
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
M1 <- matrix(c(0.00, 0.69, 0.00, 0.10, 0.00, 0.88, 0.35, 0.00, 0.79),nrow=3)
colnames(M1) <- 1:3
plotLifeCycle(M1)
} 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
