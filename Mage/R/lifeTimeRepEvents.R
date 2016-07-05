#' @export
#' @import popbio
#' @import MASS


lifeTimeRepEvents <- function(matU, matF, startLife = 1){
  #Function to determine probability of reaching reproduction, age at 
  #maturity and reproductive lifespan (Code adapted from H. Caswell's 
  #matlab code):
  
  
  uDim <- dim(matU)[1]
  surv <- colSums(matU)
  repLifeStages <- colSums(matF)
  repLifeStages[which(repLifeStages>0)] <- 1
	
if(missing(matF) | missing(matU)){stop('matU or matF missing')}
if(sum(matF,na.rm=T)==0){stop('matF contains only 0 values')}

#Probability of survival to first reprod event
	Uprime <- matU
	Uprime[ , which(repLifeStages == 1)] <- 0
	Mprime <- matrix(0, 2, uDim)
	for (p in 1:uDim[1]) {
	  if (repLifeStages[p] == 1){ Mprime[2,p] <- 1 }else{
	    Mprime[1, p] = 1 - surv[p]}
	}
	Bprime <- Mprime %*% (MASS::ginv(diag(uDim) - Uprime))
	pRep <- Bprime[2, startLife]

	out = data.frame(pRep = pRep)
  
#Age at first reproduction (La; Caswell 2001, p 124)
	D <- diag(c(Bprime[2,]))
	Uprimecond <- D %*% Uprime %*% MASS::ginv(D)
	expTimeReprod <- colSums(MASS::ginv(diag(uDim) - Uprimecond))
	La <- expTimeReprod[startLife]
  out$La <- La

#Mean life expectancy conditional on entering the life cycle in the first reproductive stage
  firstRepLifeStage <- min(which(repLifeStages == 1))
	N <- solve(diag(uDim[1]) - matU)
	meanRepLifeExpectancy <- colSums(N)[firstRepLifeStage]
  out$meanRepLifeExpectancy <- meanRepLifeExpectancy

	
#Life expectancy from mean maturity
  remainingMatureLifeExpectancy <- colSums(N)[startLife] - La
  out$remainingMatureLifeExpectancy <- remainingMatureLifeExpectancy

	return(out)
 }
