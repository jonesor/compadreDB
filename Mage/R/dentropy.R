#' @export



#' Calculate Demetrius' life table entropy
#' 
#' This function calculates Demetrius' life table entropy from an lx
#' (survivorship) vector and an mx (age-specific reproduction) or a cx
#' (age-specific clonality) vector with even intervals.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param lx A numerical vector of lx (survivorship). This is assumed to be
#' with a constant interval (e.g. 1yr).
#' @param fx A numerical vector of fx (age-specific sexual reproduction). This
#' is assumed to be with a constant interval (e.g. 1yr).
#' @param cx A numerical vector of cx (age-specific clonal reproduction). This
#' is assumed to be with a constant interval (e.g. 1yr).
#' @return Returns an estimate of Demetrius' life table entropy. When both 'fx'
#' and 'cx' are provided, it outputs Demetrius' entropy for sexual reproduction
#' only, for clonal reproduction only, and for both types of reproduction
#' together.
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references L. Demetrius. 1978. Adaptive value, entropy and survivorship
#' curves. Nature 275, 213 - 214. doi:10.1038/275213a0
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' matA <- matrix (c(0, 0, 5, 10, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
#' lx <-  c(1.0000, 0.1500, 0.0150, 0.0015, 0.0002, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)
#' mx <- c(0,  0,  5, 10, 10, 10, 10, 10, 10, 10)
#' cx <- c(0,  0,  2, 1, 1, 1, 1, 1, 1, 1)
#' 
#' dentropy(matA, lx, mx, cx)
#' 
#' 
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
