#' Calculate Keyfitz' life table entropy
#' 
#' This function calculates Keyfitz' life table entropy from an lx
#' (surivorship) vector with even intervals.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param lx A numerical vector of lx (survivorship). This is assumed to be
#' with a constant interval (e.g. 1yr).
#' @param trapeze A logical argument indicating whether the trapezoidal
#' approximation should be used for approximating the definite integral.
#' @return Returns an estimate of Keyfitz' life table entropy based on an lx
#' (survivorship) vector.
#' @note %% ~~further notes~~
#' @author 
#' Owen R. Jones <jones@@biology.sdu.dk>
#' Roberto Salguero-Gomez <salguero@@sheffield.ac.uk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  %% ~~references~~
#' @examples
#'
#' #Survivorship (lx) with constant mortality (should have Keyfitz entropy of 1).
#' 
#' x <- 0:100
#' qx <- .4
#' px <- 1-qx
#' lx <- px^x
#' 
#' plot(x,lx,type="l",ylim=c(0,1),col="red")
#' kentropy(lx)
#' kentropy(lx,trapeze=FALSE)
#' 
#' 
#' @export kentropy
kentropy <- function(lx, trapeze = FALSE){
  
  if(max(lx) > 1) stop("`lx` should be bounded between 0 and 1")
  if(sum(is.na(lx))>1) stop("There are missing values in `lx`")
  #if(sum(!diff(lx) <= 0)) stop("`lx` does not monotonically decline")
  if(sum(!diff(lx) <= 0)!=0)stop("`lx` does not monotonically decline")
  
  
  if(trapeze == TRUE){
    ma <- function(x,n=2){filter(x,rep(1/n,n), sides=2)}
    lx2 <- na.omit(as.vector(ma(lx)))
    return(-sum(lx2*log(lx2))/sum(lx2))
    }else{
      return(-sum(lx*log(lx))/sum(lx))
    }
  }
