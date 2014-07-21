# ---------------------------------------------------------------------------- #
# - DESCRIPTION: Code to demonstrate subsetting of the COMPADRE database to 
#                species of interest.
# - AUTHORS:     Owen Jones
# ---------------------------------------------------------------------------- #

# Set the working directory, then load the COMPADRE data:
load("COMPADRE Jul 21 2014.RData")

# Subsetting to data of interest:
# I want to look at the mean matrrices for shrubs that are from studies that are 
# 5+ years long, with a dimension of 5+.
# I first use subset() to subset the metadata part of the compadre object to rows 
# that match myn criteria. 

tempMetadata <- subset(compadre$metadata,MatrixComposite == "Mean" & 
    GrowthType =="Shrub" & StudyDuration >= 5 & MatrixDimension >= 4)

# Now I can use the row names from the subsetted dataframe to subset the matrices.

keep <- rownames(tempMetadata)
tempMat <- compadre$mat[keep]

# The object tempMat is now a list object containing matrices in the same order 
# that their metadata appears in tempMetadata.
# I could analyse these matrices by applying functions in a loop, or by using 
# lapply.

