#' @export

longevity <- function(matU, startLife = 1, initPop = 100, run = 1000){
  #Function to calculate mean life expectancy and maximum longevity from
  # H. Caswell's matlab code, and Morris & Doak:

  if(missing(matU)){stop('matU missing')}
  
  out = NULL
  
  matDim <- dim(matU)[1]
  
  #Mean life expectancy
  N <- solve(diag(matDim) - matU)
  out$eta <- colSums(N)[startLife]
  
  #Maximum longevity up to 'run' interations.
  popVector <- c(rep(0,matDim))
  popVector[startLife] <- initPop
  lifespanLeftover=matrix(0,run,1)
  for (n in 1:1000)	{
    lifespanLeftover[n]=sum(popVector)
    popVector=matU%*%popVector
  }
  
  out$Max <- max(which(lifespanLeftover>1))
  if (out$Max==Inf) {out$Max=run}
  
	return(out)
 }
