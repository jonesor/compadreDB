# Function to convert the matrix database into the flat-sheet format 
# used by Dave Hodgson's group in University of Exeter.
#
# use: x<-convert2flat(compadre)

convert2flat <- function(db=compadre, Aonly = TRUE){
 
  db$metadata$Amatrix <- NULL
  for (i in 1:nrow(db$metadata)){
    db$metadata$Amatrix[i] <- paste("[",paste(db$mat[[i]][[1]],collapse=" "),"]",sep="")
  }
 
  return(db$metadata)
}


