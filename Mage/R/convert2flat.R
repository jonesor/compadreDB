#' @export

convert2flat <- function(db, Aonly = TRUE){
 
  db$metadata$Amatrix <- NULL
  for (i in 1:nrow(db$metadata)){
    db$metadata$classnames[i] <- paste(db$matrixClass[[i]]$MatrixClassAuthor,collapse=" | ")
    db$metadata$matrixA[i] <- paste("[",paste(t(db$mat[[i]]$matA),collapse=" "),"]",sep="")
    db$metadata$matrixU[i] <- paste("[",paste(t(db$mat[[i]]$matU),collapse=" "),"]",sep="")
    db$metadata$matrixF[i] <- paste("[",paste(t(db$mat[[i]]$matF),collapse=" "),"]",sep="")
    db$metadata$matrixC[i] <- paste("[",paste(t(db$mat[[i]]$matC),collapse=" "),"]",sep="")
  }
 
  return(db$metadata)
}
