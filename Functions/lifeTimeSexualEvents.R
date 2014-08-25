lifeTimeSexualEvents <- function(matU, matF, conv = 0.05, startLife = 1){
  #Function to determine probability of reaching reproduction, age at maturity and reproductive lifespan (Code adapted from H. Caswell's matlab code):
  
  uDim = dim(matU)[1]
  surv = colSums(matU)
  reprodLifeStages = colSums(matF)
  reprodLifeStages[which(reprodLifeStages>0)] = 1
	
if(!missing(matF)){
     if(sum(matF,na.rm=T)==0) stop('matF contains only 0 values')

#Probability of survival to first reprod event
	Uprime=matU
	Uprime[,which(reprodLifeStages==1)]=0
	Mprime=matrix(0,2,uDim)
	for (p in 1:uDim[1]) {
	  if (reprodLifeStages[p]==1) Mprime[2,p]=1 else
	    Mprime[1,p]=1-surv[p]
	}
	Bprime=Mprime%*%(ginv(diag(uDim)-Uprime))
	pFirstSexualReprod=Bprime[2,startLife]

	out = data.frame(prob1stReprod = prob1stReprod)
  
  #Age at first reproduction (La; Caswell 2001, p 124)
	D=diag(c(Bprime[2,]))
	Uprimecond=D%*%Uprime%*%ginv(D)
	expTimeReprod=colSums(ginv(diag(uDim)-Uprimecond))
	La=expTimeReprod[startLife]

	out = data.frame(La = La)

  #Reproductive lifespan
  
  
}