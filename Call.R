# R code to run R functions for predicting spatial effects and for mapping QTL
#-----------------------------------------------------------------------------

# For sourcing two mian functions
#--------------------------------
source("Spatial.R")
source("MappingQTL.R")
source("UtilFunctions.R")

# Main script for importing data, relationship matrix 
# (IBD and based on pedigree)
#------------------------------------------------
library("INLA")
library("coxme")
library("kinship2")
library("gdata")

# To load IBD matrices for all positions
#---------------------------------------
load("KIN_Eg9PP_10.Rdata")

# To load pedigree file
#----------------------
Ped <- read.csv("Eg9PP_Pedigree.csv", sep=";", na.string= 0)

# To load data with all palms 'data', and for fully genotyped palms 'data.mapping'
#---------------------------------------------------------------------------------
data <- read.csv("Eg9PP_Phenotypes.csv", sep=";")
n <- nrow(data)
m <- nrow(Ped)-n

data.mapping <- read.csv("Eg9PP_Phenotypes_Mapping.csv", sep=";")
n.mapping <- nrow(data.mapping)
m.mapping <- nrow(LG[[1]])-n.mapping
# For the coxme function, need to create new identifiers for individuals
#-----------------------------------------------------------------------
nom.ord <- paste(rep(letters, each=676),paste(rep(letters,each=26),rep(letters),sep=""), sep="")[1:nrow(data)]
nom.ord.mapping <- paste(rep(letters,each=26),rep(letters),sep="")[1:nrow(data.mapping)]

data$iden <- nom.ord
data.mapping$iden <- nom.ord.mapping

# For performing the model (2)
SPA <- SpatialFunction(data, Ped, nb=5.9 ,m )$T1S

# For performing scan genome with models (3) and (4)
#---------------------------------------------------
data.mapping$SPA <- SPA[data.mapping$PLOT]
res.qtl <- MappingQTL(data = data.mapping, LG,
                      type.resp = 1 ,spatial = TRUE, m.mapping)


