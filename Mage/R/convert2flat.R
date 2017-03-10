#' @export



#' A function to convert from the list structured database object to a flat
#' sheet.
#' 
#' The function converts from the list structured database object to a flat
#' sheet by converting each matrix to a string.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param db The COMPADRE or COMADRE matrix database object
#' @param Aonly A logical value (TRUE/FALSE) indicatting whether ONLY the A
#' matrix be included in the output.
#' @return A `data.frame` with the same columns as are present in the metadata
#' part of the COMPADRE/COMADRE object, followed by the string-form matrix
#' stage information and the matrices themselves.
#' @note %% ~~further notes~~
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' \dontrun{
#' newDB<-convert2flat(compadre,Aonly=FALSE)
#' }
#' 
#' 
#' @export convert2flat
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
