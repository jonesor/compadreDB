#Function to compare two COMPADRE data objects
#Author: Owen Jones
#
#Example use
#
#compadreCompare(file1="/Users/orj/Documents/Dropbox/ComPADRe intro ms/Data/COMPADRE_10_11_2014_version_3.0.RData"
#                , file2 = "/Users/orj/Documents/Dropbox/ComPADRe intro ms/Data/COMPADRE_Nov_12_2014_Version_3.1.RData")
#
#

compadreCompare <- function(file1,file2){ 

load(file1)
compadre1 <- compadre
rm(compadre)
load(file2)
compadre2 <- compadre
rm(compadre)

#Quick summary
cat("QUICK SUMMARY\n")
#File 1
uniqueSource <- unique(paste(compadre1$metadata$Authors,compadre1$metadata$Journal,compadre1$metadata$YearPublication))                      
cat(paste("File-1 contains the demographic and associated data for ", 
       length(uniqueSource), " studies, corresponding to ",
       length(unique(compadre1$metadata$SpeciesAccepted))," accepted species, (",length(unique(compadre1$metadata$SpeciesAuthor)) ," according to authors) and ",
       nrow(compadre1$metadata), " matrices.\n\n",sep=""))

#File 2
uniqueSource <- unique(paste(compadre2$metadata$Authors,compadre2$metadata$Journal,compadre2$metadata$YearPublication))                      
cat(paste("File-2 contains the demographic and associated data for ", 
      length(uniqueSource), " studies, corresponding to ",
      length(unique(compadre2$metadata$SpeciesAccepted))," accepted species, (",length(unique(compadre2$metadata$SpeciesAuthor)) ," according to authors) and ",
      nrow(compadre2$metadata), " matrices.\n\n",sep=""))

cat("DETAILED SUMMARY\n")
#Accepted species in File 1 that are not in File 2
sp1 <- unique(compadre1$metadata$SpeciesAccepted)
sp2 <- unique(compadre2$metadata$SpeciesAccepted)

cat("Number of accepted species in File 1\n")
print(length(sp1))

cat("Number of accepted species in File 2\n")
print(length(sp2))

cat("Accepted species in File 1 that are not in File 2\n")
print(sp1[which(!sp1%in%sp2)])

cat("Accepted species in File 2 that are not in File 1\n")
print(sp2[which(!sp2%in%sp1)])


asp1 <- unique(compadre1$metadata$SpeciesAuthor)
asp2 <- unique(compadre2$metadata$SpeciesAuthor)
cat("Number of author species in File 1\n")
print(length(asp1))

cat("Number of author species in File 2\n")
print(length(asp2))

cat("Author species in File 1 that are not in File 2\n")
print(asp1[which(!asp1%in%asp2)])

cat("Author species in File 2 that are not in File 1\n")
print(asp2[which(!asp2%in%asp1)])
}



