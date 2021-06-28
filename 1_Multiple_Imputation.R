###################################################################	
### <April 20 2021> Yifan Sha and Min Lu
###
###  
###
###################################################################	
### 
### R code for "LR hunting: a Random Forests based cell-cell interaction 
###             discovery method for single-cell gene expression data"
### ------------------------------------------------------------
###  Yifan Sha, PhD                  y.sha@umiami.edu 
###  Min Lu, PhD                     m.lu6@umiami.edu
###                                  phone: 305-243-5473
###  Research Assistant Professor, Division of Biostatistics           
###  Clinical Research Building
###  1120 NW 14th Street, Room 1059
###  University of Miami, Miami FL 33136
###
###  -------------------------------------------------------------
###  THIS PROGRAM SHOULD NOT BE COPIED, USED, MODIFIED, OR 
###  DISSEMINATED IN ANY WAY WITHOUT SPECIFIC WRITTEN PERMISSION 
###  FROM THE AUTHOR.
###
###################################################################
library(randomForestSRC) 

load("LogNormalized_9types.data")
nam <- rownames(LogNormalized.data)

Lig.cell <- "Myeloid" 
Recp.cell <-"Epithelial_Basal"
pair <- 1
patient <- 1
setwd(paste("/Volumes/Extreme SSD/all/folders/pair", pair, "patient", patient, sep = ""))

Lig.id <- which((info$patientID == "P1") & (info$celltype_final == Lig.cell) )
Recp.id <- which((info$patientID == "P1") & (info$celltype_final == Recp.cell) )
setwd("/Volumes/HDD/pair7patient1")

Lig.n <- length(Lig.id)
Recp.n <- length(Recp.id) 

lig <- unique(LRdb$ligand)
recp <- unique(LRdb$receptor)

varsel.lig <- which((rownames(LogNormalized.data) %in% lig)==TRUE)
varsel.recp <- which((rownames(LogNormalized.data) %in% recp)==TRUE)

lig.p <- length(varsel.lig)
recp.p <- length(varsel.recp)

dat <- cbind(t(LogNormalized.data[varsel.lig, Lig.id]), matrix(NA, Lig.n, recp.p))
colnames(dat) <- c(paste(nam[varsel.lig], ".Lig", sep = ""),
                   paste(nam[varsel.recp], ".Recp", sep = "") )
dat <- rbind(dat, cbind(matrix(NA,Recp.n, lig.p), t(LogNormalized.data[varsel.recp, Recp.id])) )
dat <- dat[,-which(colSums(dat,na.rm = TRUE) <= 0)]
dat <- as.data.frame(dat)

rm(LogNormalized.data)

for (kk in 1:20){
  obj <- rfsrc(Unsupervised() ~., data = dat, na.action = "na.impute")
  options(rf.cores = 1)
  options(mc.cores = 1)
  pmd <- find.interaction(obj, method="maxsubtree", sorted=FALSE, verbose=FALSE)
  
  lig.P <- length(which(grepl(".Lig",colnames(pmd))==TRUE))
  recp.P <- length(which(grepl(".Recp",colnames(pmd))==TRUE))
  
  
  save(pmd, obj, Lig.n, Recp.n, lig.P, recp.P, Lig.cell, Recp.cell, file = paste(c(c(paste("Ligand (",Lig.cell,")", sep =""),
                                                                                     paste("Receptor (", Recp.cell,")", sep =""),"PMDvalue"),kk, ".RData"), collapse = "_"))
}