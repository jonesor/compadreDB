#Code to extract lx, mx and L_alpha from stage-based population projection matrices.

#Presented as part of the Supplementary Methods for the manuscript:
#Owen R. Jones, Alexander Scheuerlein, Roberto Salguero-Gómez, Carlo Giovanni Camarda, Ralf Schaible, Brenda B. Casper, Johan P. Dahlgren, 
#Johan Ehrlén, María B. García, Eric Menges, Pedro F. Quintana-Ascencio, Hal Caswell, Annette Baudisch, James W. Vaupel
#Varieties of Ageing Across the Tree of Life. Submitted to Nature.

#Code developed by R. Salguero-Gomez (University of Queensland & Max Planck Institute for Demographic Research), O. R. Jones (University of Southern Denmark & Max-Planck Odense Center on the Biodemography of Aging) and Hal Caswell (Woods Hole Oceanographic Institute & Max-Planck Institute for Demographic Research)
#Email: Roberto Salguero-Gomez <r.salguero@uq.edu.au>
#Using equations from H. Caswell (2001) Matrix Population Models. 2nd Edition. Sinauer, Sunderland, MA. Specific equations and pages within the references are cited in each of the functions below.

# Last modified: August 11th, 2013

# For this script, the matrix data must be arranged in a CSV as follows:
# The first column, "classOrganize" indicates the stage-type for each row/column of the matrix.
# This is a re-classification of the stages indicated by the authors (see variable "classAuthor" below) according to five possible general classes:
#                 - prop: propagule that has not yet been established (seeds in the seedbank of a plant, spores in some sessile animals)
#                 - pre_rep: pre-reproductive stages.
#                 - rep: sexually reproductive stages.
#                 - post_rep: post-reproductive stages, where the individual is not dormant/hybernating.
#                 - dorm: dormant/hybernating individuals. By default individuals in this class are non-reproductive.
# The second column, "classAuthor" gives the description of the stages in the population matrix model as defined by the author in the pertinent publication.
# The second column, "classNumber" gives the numerical ordination of the classes, from 1 to n, where n is the last class.
# There then follow n rows each for the U, F and C matrices. U gives the survival probabilities while F and C give the sexual reproduction and clonal reproduction respectively
#   - U1-Un: matrix of survival probabilities (u[i,j]), where n is the largest class number in the matrix model.
#   - F1-Fn: matrix of per-capita sexual contributions (f[i,j]), where n is the largest class number in the matrix model.
#   - C1-Cn: matrix of per-capita clonal contributions (c[i,j]), where n is the largest class number in the matrix model.
# 
# Note: there are several examples in the Supplementary Information (Source Data).

#Packages necessary:
require(MASS, popbio)

#Read in the csv file containing the matrix for the species of interest
Dataset = read.csv("PATH TO MATRIX CSV FILE")

#Define possiblie class categories
Classes = c("prop", "pre_rep", "rep", "post_rep", "dorm")

#Simple re-statement of each of the variables as "character" or as "numeric" for later calculations.
for(i in 1:dim(Dataset)[2]){
  if(i < which(colnames(Dataset) == "classNumber")){Dataset[ , i] = as.character(Dataset[ , i])
    }else{
      Dataset[,i]=as.numeric(Dataset[,i])
    }
}

#Matrix dimension:
matDim = dim(Dataset)[1]
	
#Read the different sub-matrices:
  #U matrix (transition probabilities):
	  Umat = as.matrix(Dataset[ , which(colnames(Dataset) == "U1"):(which(colnames(Dataset) == "U1") + matDim - 1)])
  #F matrix (per-capita sexual contributions):
    Fmat = as.matrix(Dataset[ , which(colnames(Dataset) == "F1"):(which(colnames(Dataset) == "F1") + matDim - 1)])
  #C matrix (per-capita asexual contributions):
    Cmat = as.matrix(Dataset[ , which(colnames(Dataset) == "C1"):(which(colnames(Dataset) == "C1") + matDim - 1)])
  #The full projection matrix is thus calculated as:
    Amat = Umat + Fmat + Cmat
  
#The life stages of this model are:
lifeStages = Dataset[ , "classOrganized"]
     
#The calculations here employed define the beginning of life when an individual become established. Thus, we do not consider transitions from the "prop" stages
notProp = min(which(lifeStages != "prop"))

#Mean life expectancy based on the fundamental matrix (See bottom equation of page 118 in Caswell 2001):
N = solve(diag(matDim) - Umat)

#The life expectancy conditional on entering the life-cycle of the species in the first non-propagule stage as per the fundamental matrix of A, N, is:
lifespanFundamental = colSums(N)[notProp]
  
#Age-specific survivorship (lx) (See top function on page 120 in Caswell 2001):
    	Umat2 = Umat
    	survivorship <- array(NA, dim = c(1000, matDim))
      for (o in 1:1000){
        survivorship[o, ] = colSums(Umat2 %*% Umat)
        Umat2 = Umat2 %*% Umat
      }
     
      lx = survivorship[, notProp]
      
      #The following line makes sure that survivorship at age 0 is 1:
      lx = c(1, lx[1:(length(lx) - 1)])
  
  #Probability of survival to first sexual reproductive event (See Eq 5.19 onwards on page 114 of Caswell 2001)
      u = colSums(Umat)
    	Uprime = Umat
    	Uprime[, (lifeStages == "rep")] = 0
    	Mprime = matrix(0, 2, matDim)
    	for (p in 1:matDim){
        if (lifeStages[p]=="flow") Mprime[2,p]=1
        }else{
        Mprime[1, p] = 1 - u[p]
    	}

    	Bprime = Mprime %*% (ginv(diag(matDim) - Uprime))
    	prob1stReprod = Bprime[1, notProp]
      
  #Mean age at sexual maturity  (See pages 124-125 in Caswell 2001)
    	D = diag(c(Bprime[2, ]))
    	Uprimecond = D%*% Uprime %*% ginv(D)
    	expTimeReprod = colSums(ginv(diag(matDim) - Uprimecond))
    	La=expTimeReprod[notProp]
      
	#Age-specific fertility (mx, Caswell 2001, p. 120)
    	ageFertility = array(0, dim = c(1000, matDim))
      fertMatrix = array(0, dim = c(1000, matDim))
     	Umat3 = Umat
    	e = matrix(rep(1, matDim))
    	for (q in 1:1000) {
    	  fertMatrix = Fmat %*% Umat3 * (as.numeric((ginv(diag(t(e) %*% Umat3)))))
    	  ageFertility[q, ] = colSums(fertMatrix)
    	  Umat3 = Umat3 %*% Umat
    	}  
      mx = ageFertility[, notProp]

  #The following line ensures that mx at age 0 is 0:
    mx = c(0, mx[1:(length(mx) - 1)])

#Function to determine the cutoff age at quasi-convergence for lx and mx (Code adapted from H. Caswell's matlab code):
    qsdConvergence <- function(survMatrix, beginLife){
      uDim = dim(survMatrix)
      eig = eigen.analysis(survMatrix)
      qsd = eig$stable.stage
      qsd = as.numeric(t(matrix(qsd / sum(qsd))))
      
      #Set up a cohort
      nzero = rep(0, uDim[1]) #Set a population vector of zeros
      nzero[beginLife] = 1 #Set the first stage to =1
      n = nzero #Rename for convenience
      
      #Iterate the cohort (n= cohort population vector, p = proportional structure)
      dist = p = NULL
      survMatrix1 <- survMatrix
      for (j in 1:1500){ #j represent years of iteration
        p = n / sum(n) #Get the proportional distribution
        dist[j] = 0.5 * (sum(abs(p - qsd)))
        n = survMatrix1 %*% n #Multiply the u and n matrices to iterate
      }
      #Find the ages for convergence to 0.1, 0.05, and 0.01
          pick1 = min(which(dist < 0.1))
          pick2 = min(which(dist < 0.05))
          pick3 = min(which(dist < 0.01))
          convage = c(pick1, pick2, pick3)
        return(convage) 
    }
    
  Convergence=qsdConvergence(Umat, notProp)
  
#Save the age-specific trajectories as a lifetable, and write out as a CSV file.
  lifetable = matrix(NA, nrow = 1000, ncol = 5)
  colnames(lifetable) = c("x", "lx", "mx", "La", "QSD")
  lifetable[ , "x"] = c(0:999)
  lifetable[ , "lx"] = lx
  lifetable[ , "mx"] = mx
  lifetable[1, "La"] = La
  lifetable[1:3, "QSD"] = Convergence
  
	write.csv(lifetable, "Lifetable.csv", row.names = FALSE)


