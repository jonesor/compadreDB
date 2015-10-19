#Plot_LifeCycle_Diagram.R
#Author: Owen Jones (jones@biology.sdu.dk) & Rob Salguero-Gomez (r.salguero@uq.edu.au)
#
#About: This script is an attempt to allow automatic drawing of life cycle 
#diagrams based on matrix projection model data. It works pretty well for 
#simple life cycles, but struggles to produce pretty diagrams for complex ones.

#Load required library
library(DiagrammeR)

#Load data (you could also simply load an individual matrix of course).
load(file ="COMADRE_v.1.0.0.RData")

#Select a matrix by index
matNum <- 64

#Obtain the species name (removing underscore)
sp <- gsub("_"," ",comadre$metadata$SpeciesAccepted[matNum])

#Identify the stages
A<-comadre$mat[[matNum]]$matA
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
temp<- subset(fromTo,quantity > 0)
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
allEdges <- paste(temp[,1],"->",temp[,2],"[minlen=",temp[,"minLVal"],",fontsize=9,color=grey,xlabel=",
                  paste("\"",round(temp[,3],3)),"\"]\n",collapse="")

#The graphviz argument, pasted together
grViz(paste(
"digraph {
  {
    graph[overlap=false];
    rank=same;
node [shape=egg, fontsize=11];
    
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
    label=\"",sp,"\"
}"
))



