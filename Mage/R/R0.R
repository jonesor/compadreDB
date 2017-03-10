#' @export



#' Calculate net reproductive value
#' 
#' This function calculates net reproductive value from a matU
#' (survival-dependent processes) and a matF (sexual reproduction) and/or a
#' matC (clonal reproduction.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param matU A matrix containing only survival-dependent processes ( growth,
#' stasis, shrinkage).
#' @param matF A matrix containing only sexual reproduction, with zeros
#' elsewhere.
#' @param matC A matrix containing only clonal reproduction, with zeros
#' elsewhere.
#' @return Returns the net reproductive value of the matrix. When both 'matF'
#' and 'matC' are provided, it outputs the net reprodutive value for sexual
#' reproduction only, for clonal reproduction only, and for both types of
#' reproduction together.
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' 
#' Hal Caswell <h.caswell@@uva.nl>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#' Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#' 978-0878930968
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' matU <- matrix (c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
#' matF <- matrix (c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
#' matC <- matrix (c(0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
#' 
#' R0(matU, matF, matU)
#' 
#' 
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
