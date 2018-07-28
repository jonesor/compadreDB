#This code is part of the supplementary information for Qian & Jin 2016.
#Qian, H. & Jin, Y. (2016) An updated megaphylogeny of plants, a tool for
#generating plant phylogenies and an analysis of phylogenetic community
#structure. Journal of Plant Ecology, 9, 233–239.
#
#
# S.PhyloMaker This repository includes four files: R_codes for S.PhyloMaker,
# PhytoPhylo, nodes, and example.splist, in addition to the Readme file. Below
# are descriptions for each of these files.
#
# 1. R_codes for S.PhyloMaker S.PhyloMaker is a function that generates
# phylogenies for seed plants. (Note: Before running this function, please
# install and load the "phytools" package.)
#
# 1.1 Description This R function (i.e., S.PhyloMaker) generates phylogenies for
# seed plants under three scenarios, which are described in the paper. This
# function should be used in conjunction with PhytoPhylo, a megaphylogeny for
# vascular plants.
#
# 1.2 Usage S.PhyloMaker(splist, tree, nodes, output.splist = T, scenarios =
# c("S1", "S2", "S3"))
#
# 1.3 Arguments splist: a user specified data frame which contains a list of
# seed plant species for which S.PhyloMaker generates phylogenies. This data
# frame should include three columns with field names of species, genus, and
# family. The three columns are, respectively, for species name (e.g. Acacia
# berlandieri), genus name (e.g. Acacia), and family name (e.g. Fabaceae).
#
# tree: a phylo object; the megaphylogeny of vascular plants (i.e. PhytoPhylo).
#
# nodes: a data.frame; it contains the basal node information of every family
# and genus on the megaphylogeny (i.e. PhytoPhylo).
#
# output.splist: logical; the output includes the columns of the "splist" data
# frame and an extra column showing which species in "splist" have matched with
# species in the megaphylogeny.
#
# scenarios: optional; defining which scenario(s) will be used to generate
# phylogenies. By default, three phylogenies corresponding to the three
# scenarios are generated.
#
# 1.4 Details Using the PhytoPhylo megaphylogeny as a backbone, this function
# takes a user-specified species list and matches species in the species list
# with species in the megaphylogeny. If a species is found in the megaphylogeny,
# the species will be selected for pruning; if the species is not found but its
# genus or family is found in the megaphylogeny, the species will be added to
# the backbone at a certain place within the genus or family depending on which
# scenarios are taken. Finally, the selected and added species are pruned from
# the megaphylogeny, which results in three phylogenies corresponding to the
# three scenarios (when default is chosen).
#
# 1.5 Value The output of a run on the function contains a set of three
# phylogenies corresponding to the three scenarios (when default is chosen) and
# a species list, which includes the original columns of the "splist" data frame
# and an additional column (called "status") showing which species in the user’s
# species list have matched with species in the megaphylogeny. Species matched
# between the species list and the megaphylogeny are indicated as
# "match(prune)"; species that are added to the megaphylogeny before the pruning
# are indicated as "match(add)"; species whose families are not found in the
# megaphylogeny are indicated as "unmatch".
#
# 1.6 Example library("phytools")                       # load the "phytools"
# package. example<-read.csv("example.splist.csv",header=T)       # read in the
# example species list. phylo<-read.tree("PhytoPhylo.tre")      # read in the
# megaphylogeny. nodes<-read.csv("nodes.csv",header=T)     # read in the nodes
# information of the megaphylogeny. result<-S.PhyloMaker(spList=example,
# tree=phylo, nodes=nodes)      # run the function S.PhyloMaker. str(result)
# # the structure of the ouput of S.PhyloMaker. par(mfrow=c(1,3),mar=c(0,0,1,0))
# # show the phylogenies of the three scenarios.
# plot(result$Scenario.1,cex=1.1,main="Scenarion One")
# plot(result$Scenario.2,cex=1.1,main="Scenarion Two")
# plot(result$Scenario.3,cex=1.1,main="Scenarion Three")
#
# 2. PhytoPhylo This is the megaphylogeny that is used as the backbone for
# S.PhyloMaker to generate phylogenies. The content of PhytoPhylo is also
# included in Appendix S3.
#
# 3. nodes This is the data frame that contains the information of the basal
# node information of every family and genus on the megaphylogeny (i.e.
# Phytophylo), as described in the “nodes” subsection of section 1.3.
#
# 4. example.splist This is an example of species list on which S.PhyloMaker can
# run. It includes three columns, which are species, genus and family names.
#
# Citation: Qian, H. and Y. Jin. (2016) An updated megaphylogeny of plants, a
# tool for generating plant phylogenies and an analysis of phylogenetic
# community structure. Journal of Plant Ecology 9(2): 233–239.



S.PhyloMaker<-function (tree, spList, nodes, output.spList = T, scenarios = c("S1", "S2", "S3")) 
{
    options(scipen=999)
    spList[sapply(spList, is.factor)] <- lapply(spList[sapply(spList, is.factor)], as.character)
    if (any(duplicated(spList$species))) 
      {
        warning("Duplicated species detected and removed.")
        print(spList$species[duplicated(spList$species)])
      }
    spList <- spList[!duplicated(spList$species), ]
    spList.original <- spList
    spList$species <- gsub(" ", "_", spList$species)
    spList$species <- gsub("(^[[:alpha:]])", "\\U\\1", spList$species, perl = TRUE)
    spList$genus <- gsub("(^[[:alpha:]])", "\\U\\1", spList$genus, perl = TRUE)
    spList$family <- gsub("(^[[:alpha:]])", "\\U\\1", spList$family, perl = TRUE)
    rnN <- data.frame(node.label = paste("N", 1:length(tree$node.label), sep = ""), oriN = tree$node.label, stringsAsFactors = FALSE)
    nodes[,c("level","family","genus","rn","bn","taxa")]<-lapply(nodes[,c("level","family","genus","rn","bn","taxa")], as.character)
    tree$node.label <- paste("N", 1:length(tree$node.label), sep = "")
    kk<-c()
    for (i in 1:length(tree$tip.label)) {
    kk<-c(kk,substring(tree$tip.label[i],1,gregexpr("_",tree$tip.label[i])[[1]][1]-1))
    }
    m<-data.frame(num=1:length(kk),genus=kk,species=tree$tip.label)
    m<-merge(m,nodes[,c("genus","family")])
    m<-m[order(m$num),]
    mm<-merge(data.frame(family=unique(spList$family)),m)
    mm1<-spList
    mm1$mark<-paste(mm1$genus,mm1$family,sep="_")
    g0<-intersect(mm1$genus,m$genus)
    g1<-intersect(mm1$genus,mm$genus)
    g<-setdiff(g0,g1) 
    if (length(g)>0)
      {
        gg0<-m[,c("genus","family")][match(g,m$genus),];dimnames(gg0)[[2]][2]<-"family_in_PhytoPhylo"
        gg1<-mm1[,c("genus","family")][match(g,mm1$genus),];dimnames(gg1)[[2]][2]<-"family_in_spList" 
        gg<-merge(gg0,gg1)
        warning("Taxnomic classification not consistent between species list and PhytoPhylo.")
        print(gg)  
      }
    mm<-tree$tip.label[setdiff(1:length(tree$tip.label),mm$num)]
    tree<-drop.tip(tree,mm)
    add.tip <- spList[which(is.na(match(spList$species, tree$tip.label))), ]
    status <- rep("match(prune)", dim(spList)[1])
    status[which(is.na(match(spList$species, tree$tip.label)))] <- "match(add)"
    if (dim(add.tip)[1] == 0 & length(na.omit(match(spList$species, tree$tip.label))) == 0)
        stop("Incorrect format of species list.")
    if (dim(add.tip)[1] == 0 & length(na.omit(match(spList$species, tree$tip.label))) > 0)
        stop("There is no species needs to be added, just prune the species from PhytoPhylo.")
    add.tip$sort <- ""
    add.tip$sort[which(!is.na(match(add.tip$genus, nodes[nodes$level == "G", ]$genus)))] <- "G1"
    add.tip$sort[which(is.na(match(add.tip$genus, nodes[nodes$level == "G", ]$genus)) & !is.na(match(add.tip$family, nodes[nodes$level == "F", ]$family)))] <- "F1"
    add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort == "F1", ]$genus)] <- "F2"
    a <- which(add.tip$sort == "")
    if (length(a) > 0)
      {  
        print(paste("Note:", length(a), "taxa unmatch:",sep=" "))
        print(add.tip$species[a])
        status[match(add.tip$species[a], spList$species)] <- "unmatch"         
      }  
    spList.original$status <- status
    if ("S1" %in% scenarios) {
        t1 <- tree
        rnN1<-rnN
        nG <- nodes[nodes$level == "G", ]
        nF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1", ]
       if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- match(data$family[i], nF$family)
                g <- nF$gen.n[n]
                s <- nF$sp.n[n]
                if (g == 1 & s == 1) {                                                                          
                  num <- grep(nF$taxa[n], t1$tip.label)
                  len <- t1$edge.length[match(num, t1$edge[, 2])]
                  t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num, position = len)
                  nF$gen.n[n] <- g + 1
                  nF$sp.n[n] <- s + 1
                  num <- grep(nF$taxa[n], t1$tip.label)
                  t1$node.label[match(t1$edge[match(num, t1$edge[, 2]), 1], unique(t1$edge[, 1]))] <- paste("NN", t1$Nnode + 1, sep = "")
                  rnN1$node.label[match(nF$bn[n], rnN1$node.label)]<- paste("NN", t1$Nnode + 1, sep = "")
                  nF$bn[n] <- paste("NN", t1$Nnode + 1, sep = "")
                  nF$bn.bl[n] <- len
                }             
                else {
                  num <- unique(t1$edge[, 1])[match(nF$bn[n], t1$node.label)]
                  len <- nF$bn.bl[n]
                  t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num)
                }
            }
        }
        data <- add.tip[add.tip$sort != "F1", ]
       if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- grep(paste(data$genus[i], "_", sep = ""), t1$tip.label)
              if (length(n) == 1) {
                  num <- n
                  len <- t1$edge.length[match(num, t1$edge[, 2])]
                  t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num, position = len)
                  }
              if (length(n) > 1) {
                  num <- fastMRCA(t1, t1$tip.label[min(n)], t1$tip.label[max(n)])
                  len <- fastDist(t1, t1$tip.label[min(n)], t1$tip.label[max(n)])/2
                  t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num)
                  }   
               }
            }
        toDrop <- setdiff(1:length(t1$tip.label), which(!is.na(match(t1$tip.label, spList$species))))
        t1 <- drop.tip(t1, tip = toDrop)
        Re <- which(!is.na(match(t1$node.label, rnN1$node.label)))
        noRe <- which(is.na(match(t1$node.label, rnN1$node.label)))
        t1$node.label[Re] <- rnN1$oriN[match(t1$node.label, rnN1$node.label)[Re]]
        t1$node.label[noRe] <- ""
    }
    else {
        t1 <- NULL
    }
    if ("S2" %in% scenarios) {
        t2 <- tree
        rnN2<-rnN
        nG <- nodes[nodes$level == "G", ]
        nF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- match(data$family[i], nF$family)
                g <- nF$gen.n[n]
                s <- nF$sp.n[n]
                if (g == 1 & s == 1) {
                  num <- grep(nF$taxa[n], t2$tip.label)
                  len <- t2$edge.length[match(num, t2$edge[,2])] * sample((1:99)/100,1)
                  t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)    
                  nF$gen.n[n] <- g + 1
                  nF$sp.n[n] <- s + 1
                  num <- grep(data$species[i], t2$tip.label)
                  t2$node.label[match(t2$edge[match(num, t2$edge[, 2]), 1], unique(t2$edge[,1]))] <- paste("NN", t2$Nnode + 1, sep = "")
                  rnN2$node.label[match(nF$bn[n], rnN2$node.label)]<- paste("NN", t2$Nnode + 1, sep = "")
                  nF$bn[n] <- paste("NN", t2$Nnode + 1, sep = "")
                  nF$bn.bl[n] <- len
                }
                else {
                  num <- unique(t2$edge[, 1])[match(nF$bn[n], t2$node.label)]
                  len <- t2$edge.length[match(num,t2$edge[,2])] * sample((1:99)/100,1)
                  t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
                  nF$gen.n[n] <- g + 1
                  nF$sp.n[n] <- s + 1
                  num <- grep(data$species[i], t2$tip.label)
                  t2$node.label[match(t2$edge[match(num, t2$edge[, 2]), 1], unique(t2$edge[, 1]))] <- paste("NN", t2$Nnode + 1, sep = "")
                  rnN2$node.label[match(nF$bn[n], rnN2$node.label)]<- paste("NN", t2$Nnode + 1, sep = "")
                  nF$bn[n] <- paste("NN", t2$Nnode + 1, sep = "")
                  nF$bn.bl[n] <- nF$bn.bl[n]+len
                }
            }
        }
       data <- add.tip[add.tip$sort != "F1", ]
       if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- grep(paste(data$genus[i], "_", sep = ""), t2$tip.label)
              if (length(n) == 1) {
                  num <- n
                  len <- t2$edge.length[match(num, t2$edge[, 2])] * sample((1:99)/100,1)
                  t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
                  }
              if (length(n) > 1) {
                  num <- sample(n,1)
                  len <- t2$edge.length[match(num, t2$edge[, 2])] * sample((1:99)/100,1)
                  t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
                  }   
               }
         }
        toDrop <- setdiff(1:length(t2$tip.label), which(!is.na(match(t2$tip.label, spList$species))))
        t2 <- drop.tip(t2, tip = toDrop)
        Re <- which(!is.na(match(t2$node.label, rnN2$node.label)))
        noRe <- which(is.na(match(t2$node.label, rnN2$node.label)))
        t2$node.label[Re] <- rnN2$oriN[match(t2$node.label, rnN2$node.label)[Re]]
        t2$node.label[noRe] <- ""
    }
    else {
        t2 <- NULL
    }
          if ("S3" %in% scenarios) {
        t3 <- tree
        rnN3<-rnN
        nG <- nodes[nodes$level == "G", ]
        nF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- match(data$family[i], nF$family)
                g <- nF$gen.n[n]
                s <- nF$sp.n[n]
                if (g == 1 & s == 1) {                                                                          
                  num <- grep(nF$taxa[n], t3$tip.label)
                  len <- t3$edge.length[match(num, t3$edge[, 2])] * (2/3)
                  t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num, position = len)
                  nF$gen.n[n] <- g + 1
                  nF$sp.n[n] <- s + 1
                  num <- grep(nF$taxa[n], t3$tip.label)
                  t3$node.label[match(t3$edge[match(num, t3$edge[, 2]), 1], unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1, sep = "")
                  rnN3$node.label[match(nF$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
                  nF$bn[n] <- paste("NN", t3$Nnode + 1, sep = "")
                  nF$bn.bl[n] <- len
                }
                if (g == 1 & s > 1) {
                  num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
                  if ((2/3)*nF$rn.bl[n] <= nF$bn.bl[n])  { len <-(nF$rn.bl[n]-nF$bn.bl[n])/2 }
                  if ((2/3)*nF$rn.bl[n] > nF$bn.bl[n])   { len <-nF$rn.bl[n]*2/3-nF$bn.bl[n] }
                  t3 <- bind.tip(t3, tip.label = data$species[i], where = num, position = len)
                  nF$gen.n[n] <- g + 1
                  nF$sp.n[n] <- s + 1
                  num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
                  t3$node.label[match(t3$edge[match(num, t3$edge[, 2]), 1], unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1, sep = "")
                  rnN3$node.label[match(nF$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
                  nF$bn[n] <- paste("NN", t3$Nnode + 1, sep = "")
                  nF$bn.bl[n] <- len+nF$bn.bl[n]
                }
                if (g > 1) {
                  num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
                  len <- nF$bn.bl[n]
                  t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num)
                }
            }
        }
        data <- add.tip[add.tip$sort != "F1", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- grep(paste(data$genus[i], "_", sep = ""), t3$tip.label)
                if (length(n) == 1) {
                  len <- t3$edge.length[match(n, t3$edge[, 2])]/2
                  t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = n, position = len)
                  nG$sp.n[match(data$genus[i], nG$genus)] <- length(n) + 1
                  nu <- grep(paste(data$genus[i], "_", sep = ""), t3$tip.label)
                  num <- fastMRCA(t3, t3$tip.label[nu[1]], t3$tip.label[nu[2]])
                  t3$node.label[match(num, unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1,  sep = "")
                  rnN3$node.label[match(nG$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
                  nG$bn[match(data$genus[i], nG$genus)] <- paste("NN", t3$Nnode + 1, sep = "")
                  nG$bn.bl[match(data$genus[i], nG$genus)] <- len
                }
                if (length(n) > 1) {
                  num <- fastMRCA(t3, t3$tip.label[min(n)], t3$tip.label[max(n)])
                  len <- as.numeric(branching.times(t3))[match(num, unique(t3$edge[, 1]))]
                  t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num)
                }
            }
        }
        toDrop <- setdiff(1:length(t3$tip.label), which(!is.na(match(t3$tip.label, spList$species))))
        t3 <- drop.tip(t3, tip = toDrop)
        Re <- which(!is.na(match(t3$node.label, rnN3$node.label)))
        noRe <- which(is.na(match(t3$node.label, rnN3$node.label)))
        t3$node.label[Re] <- rnN3$oriN[match(t3$node.label, rnN3$node.label)[Re]]
        t3$node.label[noRe] <- ""
    }
    else {
        t3 <- NULL
    }
    if (output.spList == FALSE) 
        spList <- NULL
    phylo <- list(Scenario.1 = t1, Scenario.2 = t2, Scenario.3 = t3, Species.list = spList.original)
    phylo[sapply(phylo, is.null)] <- NULL
    return(phylo)
}

