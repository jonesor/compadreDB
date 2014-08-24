makeLifeTable<-function(matU, matF = NULL, matC = NULL, startLife = 1, nSteps = 1000){
  
  matDim = ncol(matU)
  
  #Age-specific survivorship (lx) (See top function on page 120 in Caswell 2001):
  matUtemp = matU
  survivorship = array(NA, dim = c(nSteps, matDim))
  for (o in 1:nSteps){
    survivorship[o, ] = colSums(matUtemp %*% matU)
    matUtemp = matUtemp %*% matU
  }
  
  lx = survivorship[, startLife]
  lx = c(1, lx[1:(length(lx) - 1)])
  
  #Start to assemble output object
  out = data.frame(x = 0:(length(lx)-1),lx = lx)
  
  if(missing(matF) == FALSE){
  #Age-specific fertility (mx, Caswell 2001, p. 120)
  ageFertility = array(0, dim = c(nSteps, matDim))
  fertMatrix = array(0, dim = c(nSteps, matDim))
  matUtemp2 = matU
  e = matrix(rep(1, matDim))
  for (q in 1:nSteps) {
    fertMatrix = matF %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
    ageFertility[q, ] = colSums(fertMatrix)
    matUtemp2 = matUtemp2 %*% matU
  }  
  mx = ageFertility[, startLife]
  mx = c(0, mx[1:(length(mx) - 1)])
  out$mx = mx
  
  }
  
  return(out)
  }