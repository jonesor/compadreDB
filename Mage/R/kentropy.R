#' @export

kentropy <- function(lx, trapeze = TRUE){
  
  if(max(lx) > 1) stop("`lx` should be bounded between 0 and 1")
  if(sum(is.na(lx))>1) stop("There are missing values in `lx`")
  if(sum(!diff(lx)) > 0) stop("`lx` does not monotonically decline")

  if(trapeze == TRUE){
    ma <- function(x,n=2){filter(x,rep(1/n,n), sides=2)}
    lx2 <- na.omit(as.vector(ma(lx)))
    return(-sum(lx2*log(lx2))/sum(lx2))
    }else{
      return(-sum(lx*log(lx))/sum(lx))
    }
  }