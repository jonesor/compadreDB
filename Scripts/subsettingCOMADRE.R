# ---------------------------------------------------------------------------- #
# - FILE NAME:   subsettingCOMADRE.R    		
# - VERSION:     01
# - DATE:        May 29 2015
# - DESCRIPTION: Code with some examples to subset searches in COMADRE
# - AUTHORS:     Rob Salguero-Gomez
# ---------------------------------------------------------------------------- #


# Set the working directory, then load the COMADRE data:
dir <- setwd(...)
load(paste(dir, "COMADRE_v.1.0.0.RData", sep=""))


# Subsetting to data of interest:

# Example 1: Plot the population growth rate and damping ratio of a given subset (below) of species in COMADRE:
# Subset: mean matrrices for bony fish of studies 5+ years long, with a dimension of 4+.
# Use subset() to subset the metadata part of the compadre object to rows that match my logical criteria. 

tempMetadata <- subset(comadre$metadata,MatrixComposite == "Mean" & Class == "Actinopterygii" & StudyDuration >= 5 & MatrixDimension > 4)

# Use the row names from the subsetted dataframe to subset the matrices.

keep <- as.numeric(rownames(tempMetadata))

# The object tempMat is now a list object containing matrices in the same order that their metadata appears in tempMetadata.

tempMat <- comadre$mat[keep]

# These matrices can now be analyzed by applying functions in a loop, or by using lapply.

# To calculate population growth rate and damping ratio for the subset matrices, first crate a bummy variable to accommodate the output
output <- data.frame(lambdas = rep(NA, length(tempMat)), damps = rep(NA, length(tempMat)))

# Popbio is necessary to easily derive demographic output (but it can be done programmatically elsewhere)
require(popbio)

# Loop to examine each matrix
for (i in 1:length(tempMat)){
 tryCatch({
    output$lambdas[i] <- max(Re(eigen(tempMat[[i]]$matA)$value))
    output$damps[i] <- damping.ratio(tempMat[[i]]$matA)
    print(paste("Species number ", i,": ", tempMetadata$SpeciesAuthor[i], sep = ""))
      }, error = function(e){})
}

#Plotting population growth rates and damping ratio of the subset matrices. The vertical, dashed red line indicates population growth rate = 1.
par(mfrow = c(1,2))
hist(log(output$lambdas), xlab = "Log population growth rate", col = "gold", main = "")
    abline(v=0,col = "red", lwd = 4, lty = 3)
hist(output$damps, xlab = "Damping ratio", col = "brown", main = "")


# Example 2: to plot on a world map populations of a given subset (below) by viability (population growth rate >, =, < 1)
# Subset mean matrices for all Carnivora in the wild in the northern hemisphere, with no issues for survival >1, for which matrices have been split into A = U + F + C, and for which reproduction was modeled

tempMetadata <- subset(comadre$metadata, MatrixComposite == "Mean" & Order == "Carnivora" & MatrixCaptivity == "W" & LatNS == "N" & SurvivalIssue < 1 & MatrixSplit == "Divided" & MatrixFec == "Yes")

keep <- as.numeric(rownames(tempMetadata))

tempMat <- comadre$mat[keep]

output <- data.frame(lambdas = rep(NA, length(tempMat)))

#Create dummy variables to convert geographic information to be plotted in map (below)
tempMetadata$LAT=NA
tempMetadata$LON=NA

for(i in 1:dim(tempMetadata)[1]){
  #if(tempMetadata$LatNS[i] == "S") {GPS$LatDeg[i] <- -GPS$LatDeg[i]}
  if(tempMetadata$LonWE[i] == "W") {tempMetadata$LonDeg[i] <- -tempMetadata$LonDeg[i]}
  tempMetadata$LAT[i] <- tempMetadata$LatDeg[i] + tempMetadata$LatMin[i]/60 + tempMetadata$LatSec[i]/3600
  tempMetadata$LON[i] <- tempMetadata$LonDeg[i] + tempMetadata$LonMin[i]/60 + tempMetadata$LonSec[i]/3600
}

# Create dummy variable to accommodate output from lambda calculations
tempMetadata$lambdas <- NA

# Popbio is necessary to easily derive demographic output (but it can be done programmatically elsewhere)
require(popbio)

# Loop to examine each matrix
for (i in 1:length(tempMat)){
 tryCatch({
    tempMetadata$lambdas[i] <- max(Re(eigen(tempMat[[i]]$matA)$value))
    print(paste("Species number ", i,": ", tempMetadata$SpeciesAuthor[i], sep = ""))
      }, error = function(e){})
}

# Creating color-codes for blue (lambda > 1), white (lambda = 1), and red (lambda < 1)
tempMetadata$LambdaCol <- NA
tempMetadata$LambdaCol[which(tempMetadata$lambdas > 1)]="Blue"
tempMetadata$LambdaCol[which(tempMetadata$lambdas == 1)]="White"
tempMetadata$LambdaCol[which(tempMetadata$lambdas < 1)]="Red"

# Maps is necessary to plot the output
require(maps)
map("world",col="gray",fill=T,bg="light blue",xlim=c(-175,176),ylim=c(-60,85),border="white",mar=rep(0,4))
# GPS data are jittered a bit for spatial visualization
points(jitter(tempMetadata$LON,amount=.6),jitter(tempMetadata$LAT,amount=.6),col="black",bg=alpha(tempMetadata$LambdaCol,.5),ylim=c(-90,90),xlim=c(-180,180),type="p",pch=21,cex=0.7)
