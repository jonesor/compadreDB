#This script checks the taxonomy in COMADRE dataset against
#Catalogue of Life (www.catalogueoflife.org)

#Load the data
load("COMADRE_v.1.0.0.RData")

#Load the required library
library(taxize)

#Turn SpeciesAuthor into latin binomial (i.e. without subsp. sp. _2 etc.)
comadre$metadata$SpeciesBinomial <- paste(comadre$metadata$GenusAccepted,comadre$metadata$SpeciesEpithetAccepted)

#Some species did not have an epithet (e.g. Tribolium sp.), for these, the epithet is listed as NA.
#This removes the NA to allow the search of CoL.
comadre$metadata$SpeciesBinomial <- gsub(" NA","",comadre$metadata$SpeciesBinomial)

#Make a temporary dataframe containing the taxonomy from COMADRE
temp <- unique(comadre$metadata[,c("SpeciesBinomial","GenusAccepted","Family","Order","Class","Phylum","Kingdom")])

#Make a dataframe to hold the results of the CoL check
checkResults<-as.data.frame(matrix(NA,nrow=nrow(temp),ncol=7))
names(checkResults) <- c("SpeciesBinomial","Genus","Family","Order","Class","Phylum","Kingdom")

#Within a loop, look up the taxon on CoL and retrieve the results
for(i in 1:nrow(temp)){
  x<-classification(temp$SpeciesBinomial[i],db='col')
  
  #Print results to screen
  print(x)
  cat("\n")
  
  #Add the result to the appropriate row/column of the checkResults dataframe,
  #If there is no result found in CoL, insert NA.
  #Note, some user attention is required!
  checkResults[i,1]<-temp$SpeciesBinomial[i]
  if(!sum(is.na(x[[1]]))>0){
    checkResults[i,7] <- subset(x[[1]],rank=="Kingdom")$name
    checkResults[i,6] <- subset(x[[1]],rank=="Phylum")$name
    checkResults[i,5] <- subset(x[[1]],rank=="Class")$name
    checkResults[i,4] <- subset(x[[1]],rank=="Order")$name
    checkResults[i,3] <- subset(x[[1]],rank=="Family")$name
    checkResults[i,2] <- subset(x[[1]],rank=="Genus")$name
  }else{checkResults[i,2:7] <- rep(NA,6)}
}

#You can now compare the results from CoL with those already in COMADRE.
checkSummary <- NULL
for(i in 1: nrow(temp)){
  checkSummary[i]<- sum(checkResults[i,]==temp[i,],na.rm=TRUE)
}

#Which taxonomies are not the same as in COL
#Compare the following two subsets of data to see which species have 
#different taxonomies than in CoL, and why.
#In COMADRE:
temp[which(checkSummary != 7),]

#The results from CoL:
checkResults[which(checkSummary != 7),]

#You can check individual entries like this:
classification("Daphnia magna",db='col')
classification("Platynothrus peltifer",db='col')



