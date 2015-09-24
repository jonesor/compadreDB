# Function to convert the matrix database into the flat-sheet format 
# used by Dave Hodgson's group in University of Exeter.
#
# use: x<-convert2flat(compadre)
# write.csv(x,file = "compadreFlat.csv",row.names=FALSE)


convert2flat <- function(db=compadre, Aonly = TRUE){
 
  db$metadata$Amatrix <- NULL
  for (i in 1:nrow(db$metadata)){
    db$metadata$classnames[i] <- paste(db$matrixClass[[i]]$MatrixClassAuthor,collapse="| ")
    db$metadata$matrixA[i] <- paste("[",paste(db$mat[[i]]$matA,collapse=" "),"]",sep="")
    db$metadata$matrixU[i] <- paste("[",paste(db$mat[[i]]$matU,collapse=" "),"]",sep="")
    db$metadata$matrixF[i] <- paste("[",paste(db$mat[[i]]$matF,collapse=" "),"]",sep="")
    db$metadata$matrixC[i] <- paste("[",paste(db$mat[[i]]$matC,collapse=" "),"]",sep="")
  }
 
  return(db$metadata)
}
