#' @export

dentropy <- function(matA, lx, fx, cx=F){
  
  #Demetrius entropy (S):
  if(max(lx) > 1) stop("`lx` should be bounded between 0 and 1")
  if(sum(is.na(lx))>1) stop("There are missing values in `lx`")
  if(sum(!diff(lx) <= 0)!=0)stop("`lx` does not monotonically decline")
  if(sum(is.na(mx))>1) print("There are missing values in `mx`")
  if(sum(is.na(cx))>1) print("There are missing values in `cx`")
  
  dentropy=NULL
  r <- log(max(Re(eigen(matA)$value)))  #population growth rate in log scale (rmax)
  
  if (sum(mx)>0){
    limiteFx <- min(length(mx[which(!is.na(mx))]), length(lx[which(!is.na(lx))]))   #Calculating the last time interval at which both mx and lx were able to be calculated
    lxmx <- lx[1:limiteFx]*mx[1:limiteFx]  #Step-by-step multiplication of lx and mx
    lxmx[which(lxmx==0)] <- 1 #Coertion of the first few values for which fx = 0 to take values lxmx = 1, so that it can be log-scaled, below
    loglxmx <- log(lxmx)
    loglxmx[which(lxmx==0)] <- NA
    
    dentropy$Fx <- abs(sum(lxmx*loglxmx)/sum(lxmx))
  }
  
  if (sum(cx)>0){
    limiteCx <- min(length(cx[which(!is.na(cx))]), length(lx[which(!is.na(lx))]))   #Calculating the last time interval at which both cx and lx were able to be calculated
    lxcx <- lx[1:limiteCx]*cx[1:limiteCx]  #Step-by-step multiplication of lx and cx
    lxcx[which(lxcx==0)] <- 1 #Coertion of the first few values for which cx = 0 to take values lxcx = 1, so that it can be log-scaled, below
    loglxcx <- log(lxcx)
    loglxcx[which(lxcx==0)] <- NA
    
    dentropy$Cx <- abs(sum(lxcx*loglxcx)/sum(lxcx))
  }
  
  if (sum(mx)>0 & sum(cx)>0){
    limiteMxCx <- min(length(mx[which(!is.na(mx))]), length(cx[which(!is.na(cx))]), length(lx[which(!is.na(lx))]))   #Calculating the last time interval at which both mx, cx and lx were able to be calculated
    mxcx=rowSums(cbind(mx,cx))
    lxmxcx <- lx[1:limiteMxCx]*mxcx[1:limiteMxCx]  #Step-by-step multiplication of lx and cx
    lxmxcx[which(lxmxcx==0)] <- 1 #Coertion of the first few values for which cx = 0 to take values lxcx = 1, so that it can be log-scaled, below
    loglxmxcx <- log(lxmxcx)
    loglxmxcx[which(lxmxcx==0)] <- NA
    
    dentropy$MxCx <- abs(sum(lxmxcx*loglxmxcx)/sum(lxmxcx))
  }
return(dentropy)
}