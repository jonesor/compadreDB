#' @export
#' @import DiagrammeR



#' A function to plot the life cycle diagram based on the matrix model
#' 
#' This function plots the life cycle diagram illustrated by a matrix model. It
#' processes the matrix model and passes the information to the graphViz
#' functionality of DiagrammeR. See http://rich-iannone.github.io/DiagrammeR/.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param A The A matrix of a matrix population model
#' @param title A title for the plot
#' @param shape The shape to be used for the stages of the diagram. Any node
#' shape accepted by graphViz is acceptable.
#' @param fontsize Size of the font used in the diagram.
#' @param nodefontsize Size of the font used in the node part of the diagram.
#' @param edgecol Colour of the arrows in the diagram.
#' @return Produces a plot consisting of a life cycle diagram, based on the A
#' matrix.
#' @note %% ~~further notes~~
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' \dontrun{
#' M1 <- matrix(c(0.00, 0.69, 0.00, 0.10, 0.00, 0.88, 0.35, 0.00, 0.79),nrow=3)
#' colnames(M1) <- 1:3
#' plotLifeCycle(M1)
#' } 
#' 
#' 
#' @export plotLifeCycle
plotLifeCycle <- function(A,title="my life cycle",shape="egg",fontsize=10,nodefontsize=12,edgecol="grey"){


#Identify the stages
rownames(A)<-colnames(A)

#Remove the unnecessary "A"
Astages<-gsub("A","",colnames(A))

#Construct a "from"->"to" dataset (edges)
fromTo<-expand.grid(Astages,Astages)
names(fromTo) <- c("From","To")

#Loop through the edges to get the quantities 
#(transition probabilities and fecundity)
for(i in 1:nrow(fromTo)){
fromTo$quantity[i] <- A[fromTo$To[i],fromTo$From[i]]  
}

#Subset to only include those where the quantity >0
temp<- subset(fromTo,fromTo$quantity > 0)
temp

#Create sorted vector of node names
allNodes <- sort(unique(c(as.character(temp[,1]),as.character(temp[,2]))))

#Add a semi-colon, for use by graphviz
allNodes <- paste(allNodes,collapse="; ")

#Manipulate minimim length of edge to make the plot pretty.
#Experimental!!
temp$minLVal <-as.numeric(temp[,2])-as.numeric(temp[,1])
temp$minLVal <- temp$minLVal*3
temp

#Create the edges argument for graphviz
#by pasting commands together.
#Note, one could modify this to alter the outputs away from my defaults.
allEdges <- paste(temp[,1],"->",temp[,2],"[minlen=",temp[,"minLVal"],",fontsize=",fontsize,",color=",edgecol,",xlabel=",
                  paste("\"",round(temp[,3],3)),"\"]\n",collapse="")

#The graphviz argument, pasted together
grViz(paste(
"digraph {
  {
    graph[overlap=false];
    rank=same;
node [shape=",shape,", fontsize=",nodefontsize,"];
    
",
allNodes
,"
  }
  ordering=out
  x [style=invis]
  x -> {",allNodes,"} [style=invis]",
  allEdges,
"
# title
    labelloc=\"t\";
    label=\"",title,"\"
}"
))
}
