# ---------------------------------------------------------------------------- #
# - FILE NAME:   subsettingCOMPADRE.R    		
# - VERSION:     01
# - DATE:        July 27 2014
# - DESCRIPTION: Code with some examples to subset searches in COMPADRE
# - AUTHORS:     Owen Jones & Rob Salguero-Gomez
# ---------------------------------------------------------------------------- #


# Set the working directory, then load the COMPADRE data:
load("COMPADRE Jul 21 2014.RData")

# Subsetting to data of interest:

# Example 1:
# I want to look at the mean matrrices for shrubs that are from studies that are 
# 5+ years long, with a dimension of 5+.
# I first use subset() to subset the metadata part of the compadre object to rows 
# that match my logical criteria. 

tempMetadata <- subset(compadre$metadata,MatrixComposite == "Mean" & 
    GrowthType =="Shrub" & StudyDuration >= 5 & MatrixDimension >= 4)

# Now I can use the row names from the subsetted dataframe to subset the matrices.

keep <- as.numeric(rownames(tempMetadata))
tempMat <- compadre$mat[keep]

# The object tempMat is now a list object containing matrices in the same order 
# that their metadata appears in tempMetadata.
# I could analyse these matrices by applying functions in a loop, or by using 
# lapply.

#The following calculates population growth rate and damping ratio for the subset matrices
output <- data.frame(lambdas=rep(NA,length(tempMat)),damps=rep(NA,length(tempMat)))

# We need this package to calculate damping ratios
require(popbio)
for (i in 1:length(tempMat)){
    output$lambdas[i] <- max(Re(eigen(tempMat[[i]]$matA)$value))
    output$damps[i] <- damping.ratio(tempMat[[i]]$matA)
    print(paste("Species number ", i,": ",tempMetadata$SpeciesAuthor[i],sep=""))
}

par(mfrow=c(1,2))
hist(log(output$lambdas), xlab = "Log population growth rate", col = "gold", main = "")
    abline(v=0,col = "red", lwd = 4, lty = 3)
hist(output$damps, xlab = "Damping ratio", col = "brown", main = "")

# Example 2:
# I want to look at the individual matrices for trees in the southern hemisphere
# where each study contained just one population

tempMetadata <- subset(compadre$metadata,MatrixComposite == "Individual" & 
    GrowthType =="Tree" & StudyDuration >= 5 & MatrixDimension >= 4
    & NumberPopulations == 1 & LatNS == "S")

keep <- as.numeric(rownames(tempMetadata))
tempMat <- compadre$mat[keep]


