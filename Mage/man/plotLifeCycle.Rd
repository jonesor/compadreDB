\name{plotLifeCycle}
\alias{plotLifeCycle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotLifeCycle(A, title = "my life cycle", shape = "egg", fontsize = 10, nodefontsize = 12, edgecol = "grey")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{shape}{
%%     ~~Describe \code{shape} here~~
}
  \item{fontsize}{
%%     ~~Describe \code{fontsize} here~~
}
  \item{nodefontsize}{
%%     ~~Describe \code{nodefontsize} here~~
}
  \item{edgecol}{
%%     ~~Describe \code{edgecol} here~~
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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

## The function is currently defined as
function (A, title = "my life cycle", shape = "egg", fontsize = 10, 
    nodefontsize = 12, edgecol = "grey") 
{
    require(DiagrammeR)
    rownames(A) <- colnames(A)
    Astages <- gsub("A", "", colnames(A))
    fromTo <- expand.grid(Astages, Astages)
    names(fromTo) <- c("From", "To")
    for (i in 1:nrow(fromTo)) {
        fromTo$quantity[i] <- A[fromTo$To[i], fromTo$From[i]]
    }
    temp <- subset(fromTo, quantity > 0)
    temp
    allNodes <- sort(unique(c(as.character(temp[, 1]), as.character(temp[, 
        2]))))
    allNodes <- paste(allNodes, collapse = "; ")
    temp$minLVal <- as.numeric(temp[, 2]) - as.numeric(temp[, 
        1])
    temp$minLVal <- temp$minLVal * 3
    temp
    allEdges <- paste(temp[, 1], "->", temp[, 2], "[minlen=", 
        temp[, "minLVal"], ",fontsize=", fontsize, ",color=", 
        edgecol, ",xlabel=", paste("\"", round(temp[, 3], 3)), 
        "\"]\n", collapse = "")
    grViz(paste("digraph {\n  {\n    graph[overlap=false];\n    rank=same;\nnode [shape=", 
        shape, ", fontsize=", nodefontsize, "];\n    \n", allNodes, 
        "\n  }\n  ordering=out\n  x [style=invis]\n  x -> {", 
        allNodes, "} [style=invis]", allEdges, "\n# title\n    labelloc=\"t\";\n    label=\"", 
        title, "\"\n}"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
