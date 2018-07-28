# Build a phylogeny for compadre.
#Script by Owen Jones, with species replacements based on work at sApropos
#workshop by Owen Jones, Jean Burns, Tiffany Knight and Judy Che-Castaldo.
#
#Description: This code builds a phylogeny based on Qian tree and the
#accompanying S.Phylomaker code that was included as supp info. Code is
#available at https://github.com/jinyizju/S.PhyloMaker, but included in this
#repository as a convenience. It would be a good idea to check the repository
#for updates before running.

#Load libraries
library(ape)
library(dplyr)
library(ggtree)
library(phytools)
library(taxize)

source("Code/S.PhyloMaker.R")

#Backbone tree and nodes information from PhytoPhylo repository
qian <- read.tree("BackboneTrees/PhytoPhylo.tre")
nodes<-read.csv("Data/nodes.csv",header=T,sep="\t") 

#Read in compadre data.
#load("/Users/jones/Dropbox/MatrixDatabases/COMPADRE/v.X.X.X/COMPADRE_v.X.X.X.RData")
load("/Users/jones/Dropbox/MatrixDatabases/COMPADRE/v.4.0.1/COMPADRE_v.4.0.1.RData")

#Obtain the species list from COMPADRE
compadre$metadata %>%
  dplyr::filter(Kingdom == "Plantae" & !is.na(SpeciesAccepted) & OrganismType!="Algae") %>%
  dplyr::select(SpeciesAccepted) %>%
  unique() -> spList 

#Corrections based on first pass TPL (these corrections should be made in
#COMPADRE for next version!)
spList %>%
  dplyr::filter(SpeciesAccepted != "Vulpicida pinastri") -> spList

spList$SpeciesAccepted <- gsub( "Mircothlaspi perfoliatum" ,"Microthlaspi perfoliatum",spList$SpeciesAccepted)
spList$SpeciesAccepted <- gsub( "Mammillaria napia" ,"Mammillaria napina",spList$SpeciesAccepted)
spList$SpeciesAccepted <- gsub( "Verticosa staminosa staminosa" ,"Verticordia staminosa staminosa",spList$SpeciesAccepted)
spList$SpeciesAccepted <- gsub( "Choerospodnias axillaris" ,"Choerospondias axillaris",spList$SpeciesAccepted)

spList <- spList$SpeciesAccepted

length(spList)

#Query The Plant List to check the taxonomy is correct.
tplResult <- Taxonstand::TPL(spList)

#Where/what are the problems?
filter(tplResult,Plant.Name.Index == FALSE)

#Check for non-matches, which indicate an error in COMPADRE
#These errors can be fixed in code here, but should be corrected in COMPADRE
#ASAP.
which(tplResult$Genus != tplResult$New.Genus)
which(tplResult$Species != tplResult$New.Species)

#process the TPL results to produce input file for the phylogeny code
tplResult %>%
  dplyr::mutate(species = paste(New.Genus,New.Species)) %>%
  dplyr::rename(genus = New.Genus, family = Family) %>%
  dplyr::select(species, genus, family) -> PhyloMakerSpeciesList

#Generate a new tree from our species list.
#This takes about 25mins on my laptop.
system.time(
  newTree<-S.PhyloMaker(qian,spList=PhyloMakerSpeciesList,nodes=nodes)
  )

#Summary of matching.
table(newTree$Species.list$status)
newTree$Species.list %>%
  dplyr::filter(status == "unmatch") -> unmatched
unmatched$species

#There are some unmatched species.
#The most parsimonious way of matching up the species with suitable candidates
#in the tree is to change the species name in the COMPADRE species list, 
#then change it back in the phylogeny produced. This is because the nodes file
#and the tree file must match. After doing this, the S.PhyloMaker code is run
#again, whereupon it should match perfectly.

speciesReplace <- data.frame(compadreSpecies = rep(NA,5), qianTreeSpecies = rep(NA,5))

#"Echinospartum ibericum"   
speciesReplace[1,] <- c("Echinospartum ibericum","Genista sp")
#"Tetraneuris herbacea"
speciesReplace[2,] <- c("Tetraneuris herbacea","Psilostrophe cooperi")
#"Periandra mediterranea"
speciesReplace[3,] <- c("Periandra mediterranea","Centrosema pubescens")
#"Argyroxiphium sandwicense"
speciesReplace[4,] <- c("Argyroxiphium sandwicense","Dubautia plantaginea")
#"Eremosparton songoricum"
speciesReplace[5,] <- c("Eremosparton songoricum","Colutea arborescens")
speciesReplace

for(i in 1:nrow(speciesReplace)){
  PhyloMakerSpeciesList$species[PhyloMakerSpeciesList$species==speciesReplace$compadreSpecies[i]] <- speciesReplace$qianTreeSpecies[i]
}

#spList genus needs to match the species, so resplit the vector.
PhyloMakerSpeciesList %>% 
  dplyr::select(-genus) %>% 
  tidyr::separate(.,species, c("genus", "species"),sep=" ",remove=F) %>%
  dplyr::mutate(species = paste(genus,species)) %>%
  dplyr::select(species,genus,family) -> PhyloMakerSpeciesList2

#Now re-run the S.PhyloMaker code based on the updated species list.
newTree<-S.PhyloMaker(qian,spList=PhyloMakerSpeciesList2,nodes=nodes)

#Now, for each tree, replace the temporary "tree" names with the original
#compadre names
for(i in nrow(speciesReplace)){
newTree$Scenario.1$tip.label[speciesReplace$qianTreeSpecies[i]] <-speciesReplace$compadreSpecies[i] 
newTree$Scenario.2$tip.label[speciesReplace$qianTreeSpecies[i]] <-speciesReplace$compadreSpecies[i]
newTree$Scenario.3$tip.label[speciesReplace$qianTreeSpecies[i]] <-speciesReplace$compadreSpecies[i]
}

#Add the species replacement information to the object that is created by
#S.PhyloMaker
newTree$speciesReplace <- speciesReplace

#Rename as compadrePhylogeny
compadrePhylogeny <- newTree

#Write out the new trees, in various formats.
save(compadrePhylogeny,file = paste("GeneratedTrees/",paste(gsub("-","",Sys.Date()),"compadrePhylogeny.RData",sep="_"),sep="")) 

write.tree(compadrePhylogeny$Scenario.1,file=paste("GeneratedTrees/",paste(gsub("-","",Sys.Date()),"compadrePhylogeny.tre",sep="_"),sep=""))
write.nexus(compadrePhylogeny$Scenario.1,file=paste("GeneratedTrees/",paste(gsub("-","",Sys.Date()),"compadrePhylogeny.nex",sep="_"),sep=""))
