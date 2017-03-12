#' @export

R0 <- function(matU, matF, matC=F){
  
  #Demetrius entropy (S):
  if(max(lx) > 1) stop("`lx` should be bounded between 0 and 1")
  if(sum(is.na(lx))>1) stop("There are missing values in `lx`")
  if(sum(!diff(lx) <= 0)!=0)stop("`lx` does not monotonically decline")
  if(sum(is.na(mx))>1) print("There are missing values in `mx`")
  if(sum(is.na(cx))>1) print("There are missing values in `cx`")
  
  R0=NULL
  matDim <- dim(matA)[1]
  N <- solve(diag(matDim)-matU)  #Fundamental matrix, which states the amount of time units spent in each stage on average
  
  if (sum(mx)>0){
    R0_matF <- matF%*%N
    R0$Fec <- R0_matF[1,1]
  }
  if (sum(cx)>0){
    R0_matC <- matC%*%N
    R0$Clo <- R0_matC[1,1]
  }
  if (sum(cx)>0){
    matFC <- matF + matC
    R0_matFC <- matFC%*%N
    R0$FecClo <- R0_matFC[1,1]
  }
  
return(R0)
}