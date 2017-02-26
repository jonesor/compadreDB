#' @export


vitalRates <- function(matU, matF, matC = FALSE, splitStages = FALSE, weighted = FALSE){
  #Function to quantify vital rates values
  
  if(missing(matU)){stop('matU missing')}
  if(missing(matF) & missing(matC)){warning('matF or matC missing. These have been coerced to matrices of zero')}
  if (sum(weighted)>0 & length(weighted) != dim(matU)[i]){stop('Population vector does not agree with matrix dimension')}
  
  matDim <- dim(matU)[1]
  matA <- matU + matF + matC
  
  surv <- colSums(matU)
  fec <- colSums(matF)
  clo <- colSums(matC)
  
  matUIndep <- matrix(NA, matDim, matDim)
  for (i in 1:matDim) {matUIndep[,i] = matU[,i]/surv[i]}
  prog <- retr <- matUIndep
  
  prog[which(upper.tri(matUIndep, diag = TRUE))]=0
  retr[which(lower.tri(matUIndep, diag = TRUE))]=0
  
  prog <- colSums(prog)
  retr <- colSums(retr)
  
  if (weighted[1] == FALSE) {
    weight <- rep(1,matDim)
  }
  
  if (weighted[1] == 'SSD') {
    weight <- Re(eigen(matA)$vectors[,which.max(Re(eigen(matA)$values))])
  }
  
  weight <- weight/sum(weight)
  
  surv1 <- surv * weight
  fec1  <- fec  * weight
  clo1  <- clo  * weight
  prog1 <- prog * weight
  retr1 <- retr * weight
  
  out = NULL
  
  if (splitStages[1] == 'all'){
    out$surv <- sum(surv1)
    out$retr <- sum(retr1)
    out$prog <- sum(prog1)
    out$fec  <- sum(fec1)
    out$clo  <- sum(clo1)
  }
  
  if (splitStages[1] == 'ontogeny'){
    adu <- which(colSums(matF)>0)  #This classification does not accoutn for non- and post-reproductive
    juv <- which(colSums(matF)==0)
    
    out$survJuv <- mean(surv1[juv], na.rm=T)
    out$retrJuv <- mean(retr1[juv], na.rm=T)
    out$progJuv <- mean(prog1[juv], na.rm=T)
    out$cloJuv  <- mean(clo1[juv], na.rm=T)
    
    out$survAdu <- mean(surv1[adu], na.rm=T)
    out$retrAdu <- mean(retr1[adu], na.rm=T)
    out$progAdu <- mean(prog1[adu], na.rm=T)
    out$fecAdu  <- mean(fec1[adu], na.rm=T)
    out$cloAdu  <- mean(clo1[adu], na.rm=T)
    }
  
  if (splitStages[1] %in% c('prop','active','dorm')){
    prop <- which(splitStages=="prop")
    active <- which(splitStages=="active")
    dorm <- which(splitStages=="dorm")
    
    out$survProp <- mean(surv1[prop], na.rm=T)
    out$progProp <- mean(prog1[prop], na.rm=T)
    
    out$survActive <- mean(surv1[active], na.rm=T)
    out$retrActive <- mean(retr1[active], na.rm=T)
    out$progActive <- mean(prog1[active], na.rm=T)
    out$fecActive  <- mean(fec1[active], na.rm=T)
    out$cloActive  <- mean(clo1[active], na.rm=T)
    
    out$survDorm <- mean(surv1[dorm], na.rm=T)
    out$retrDorm <- mean(retr1[dorm], na.rm=T)
    out$progDorm <- mean(prog1[dorm], na.rm=T)
  }
  
	return(out)
 }
