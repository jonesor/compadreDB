# ---------------------------------------------------------------------------- #
# - FILE NAME:   ternaryPlots.R        
# - VERSION:     01
# - DATE:        Sept 29 2015
# - DESCRIPTION: Code to produce ternary plots of elasticities
# - AUTHORS:     Rob Salguero-Gomez
# ---------------------------------------------------------------------------- #


# This code will produce a ternary plot a la Silvertown & Franco 1993 with various
# life history traits such as mean life expectancy, population growth rate or reactivity
# as the "fourth" dimension

# Set the working directory, then load the COMADRE data:
dir <- setwd(...)
load(paste(dir, "COMADRE_v.1.0.0.RData", sep=""))


# Subsetting to studies with a matrix dimension >= 3, matrices are mean of unmanipulated conditions,
# duration > 3 years, where sexual reproduction has been modeled
# explicitly, the matrices are split into U, F and C, and there are no issues with stage-specific survival >1

tempMetadata <- subset(comadre$metadata,
                        MatrixDimension >= 3 &
                        MatrixComposite == "Mean" &
                        MatrixTreatment == "Unmanipulated" &
                        StudyDuration > 3 &
                        MatrixFec == "Yes" &
                        MatrixSplit == "Divided" &
                        SurvivalIssue < 1)


# Use the row names from the subsetted dataframe to subset the matrices.

keep <- as.numeric(rownames(tempMetadata))

# The object tempMat is now a list object containing matrices in the same order that their metadata appears in tempMetadata.

tempMat <- comadre$mat[keep]

# These matrices can now be analyzed by applying functions in a loop, or by using lapply.

# To calculate population growth rate and damping ratio for the subset matrices, first crate a bummy variable to accommodate the output
output <- data.frame(species= rep(NA, length(tempMat)),
                      lambdas = rep(NA, length(tempMat)),
                      eta = rep(NA, length(tempMat)),
                      react = rep(NA, length(tempMat)),
                      EStasis = rep(NA, length(tempMat)),
                      EProgression = rep(NA, length(tempMat)),
                      ERetrogression = rep(NA, length(tempMat)),
                      EFecundity = rep(NA, length(tempMat)),
                      EClonality = rep(NA, length(tempMat)))

# Popbio is necessary to easily derive demographic output (but it can be done programmatically elsewhere)
require(popbio)
require(popdemo)
#Function to calculate element-level perturbations
matrixElementPerturbation <- function(matU, matF, matC=NULL,pert=0.001){
  #Function to calculate matrix element level sensitivities and elasticities
  
  matA=matU+matF+matC
  aDim=dim(matA)[1]
  fakeA=matA
  sensA=elasA=matrix(NA,aDim,aDim)
  lambda=Re(eigen(matA)$values[1])

  propU=matU/matA
    propU[is.nan(propU)]=NA
    propProg=propRetrog=propU
    propProg[upper.tri(propU,diag=T)]=NA
    propRetrog[lower.tri(propU,diag=T)]=NA
    propStasis=matrix(diag(aDim)*diag(propU),aDim,aDim)
  propF=matF/matA
    propF[is.nan(propF)]=NA
  propC=matC/matA
    propC[is.nan(propC)]=NA

  for (i in 1:aDim){
    for (j in 1:aDim){
       fakeA=matA
       fakeA[i,j]=fakeA[i,j]+pert
       lambdaPert=eigen(fakeA)$values[1]
       sensA[i,j]=(lambda-lambdaPert)/(matA[i,j]-fakeA[i,j])
    }
  }

  sensA=Re(sensA)
  elasA=sensA*matA/lambda
  
  out = data.frame("SStasis"=NA,"SProgression"=NA,"SRetrogression"=NA,"SFecundity"=NA,"SClonality"=NA,
                  "EStasis"=NA,"EProgression"=NA,"ERetrogression"=NA,"EFecundity"=NA,"EClonality"=NA)
    
    out$SStasis=sum(sensA*propStasis,na.rm=T)
    out$SRetrogression=sum(sensA*propRetrog,na.rm=T)
    out$SProgression=sum(sensA*propProg,na.rm=T)
    out$SFecundity=sum(sensA*propF,na.rm=T)
    out$SClonality=sum(sensA*propC,na.rm=T)
    out$EStasis=sum(elasA*propStasis,na.rm=T)
    out$EProgression=sum(elasA*propProg,na.rm=T)
    out$ERetrogression=sum(elasA*propRetrog,na.rm=T)
    out$EFecundity=sum(elasA*propF,na.rm=T)
    out$EClonality=sum(elasA*propC,na.rm=T)

  return(out) 
}


# Loop to examine each matrix
for (i in 1:length(tempMat)){
 tryCatch({
    matA=tempMat[[i]]$matA
    matU=tempMat[[i]]$matU
    matF=tempMat[[i]]$matF
    matC=tempMat[[i]]$matC
    output$species[i] <- tempMetadata$SpeciesAuthor[i]
    output$lambdas[i] <- max(Re(eigen(matA)$value))
      uDim=dim(matA)[1]
      N = solve(diag(uDim[1])-matU)
      output$eta[i] = eta = colSums(N)[1]
      output$react[i] <- reactivity(matA)
    output[i,c("EStasis","EProgression","ERetrogression","EFecundity","EClonality")]=matrixElementPerturbation(matU=matU,matF=matF,matC=matC)[6:10]
    print(paste("Species number ", i,": ", tempMetadata$SpeciesAuthor[i], sep = ""))
      }, error = function(e){})
}


#Plotting things on a ternary plot
output$S=output$EStasis+output$ERetrogression
output$G=output$EProgression
output$R=output$EFecundity+output$EClonality

#Scaling to 1 (rounding issues)
output$S=output$S/rowSums(output[,c("S","G","R")])
output$G=output$G/rowSums(output[,c("S","G","R")])
output$R=output$R/rowSums(output[,c("S","G","R")])


#Temporary space holder while Sterna_hirundo_2 is being fixed
#output=output[-which(is.na(output[,2])),]


#Color-coding according to range of lambda, eta and reactivity
lambdaData=output[which(log(output$lambdas)<=2),]
x_norm=log(lambdaData$lambda)
x_norm = (lambdaData$lambda - min(lambdaData$lambda)) / (max(lambdaData$lambda) - min(lambdaData$lambda))
col_fun <- colorRamp(c("white","yellow","orange","red","dark red"))
rgb_cols <- col_fun(x_norm)
colsLambda <- rgb(rgb_cols, maxColorValue = 256)

etaData=output
etaData$etalog=log(etaData$eta)
x_norm = (etaData$etalog - min(etaData$etalog)) / (max(etaData$etalog) - min(etaData$etalog))
col_fun <- colorRamp(c("white","yellow","orange","red","dark red"))
rgb_cols <- col_fun(x_norm)
colsEta <- rgb(rgb_cols, maxColorValue = 256)

reactData=output[which(log(output$react)<=5),]
reactData$reactlog=log(reactData$react)
x_norm = (output$react - min(output$react)) / (max(output$react) - min(output$react))
col_fun <- colorRamp(c("white","yellow","orange","red","dark red"))
rgb_cols <- col_fun(x_norm)
colsReact <- rgb(rgb_cols, maxColorValue = 256)


install.packages('vcd')
require(vcd)
install.packages('scales')
require(scales)


pdf("Ternary plot lambda.pdf")
  zr <- range(c(lambdaData$lambda,na.rm=T))
  colCode <- colorRampPalette(c("white","yellow","orange","red","dark red"))(n = 999)  
  image.plot(legend.only=TRUE, zlim= zr, col=colCode, smallplot=c(.75,.8, .5,.75),cex.axis=0.2) 
  ternaryplot(lambdaData[,c("R","S","G")],scale=1,col=alpha(colsLambda,0.7),bg="black", newpage=F, dimnames=c("Stasis","Growth","Reproduction"),dimnames_position="edge", main=expression(paste("Population growth rate - ", lambda,)))

dev.off()

pdf("Ternary plot mean life expectancy.pdf")
  zr <- range(c(etaData$eta,na.rm=T))
  colCode <- colorRampPalette(c("white","yellow","orange","red","dark red"))(n = 999)  
  image.plot(legend.only=TRUE, zlim= zr, col=colCode, smallplot=c(.75,.8, .5,.75),cex.axis=0.2) 
  ternaryplot(etaData[,c("R","S","G")],scale=1,col=alpha(colsEta,0.7),bg="black", newpage=F, dimnames=c("Stasis","Growth","Reproduction"),dimnames_position="edge", main=expression(paste("Mean life expectancy - ", eta["e"])))
dev.off()

pdf("Ternary plot reactivity.pdf")
  zr <- range(c(reactData$react,na.rm=T))
  colCode <- colorRampPalette(c("white","yellow","orange","red","dark red"))(n = 999)  
  image.plot(legend.only=TRUE, zlim= zr, col=colCode, smallplot=c(.75,.8, .5,.75),cex.axis=0.2) 
  ternaryplot(reactData[,c("R","S","G")],scale=1,col=alpha(colsReact,0.7),bg="black", newpage=F, dimnames=c("Stasis","Growth","Reproduction"),dimnames_position="edge",main=expression(paste("Reactivity - ||", hat(A),"||"[1])))
dev.off()

