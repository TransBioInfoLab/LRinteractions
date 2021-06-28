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
pvalues <- function(pmd, x){
  xx <- pmd
  diag(xx) <- NA
  xx <- sort(na.omit(c(xx)))
  n <- length(xx)
  
  unlist(lapply(1:length(x), function(i){
   max(which(xx <= x[i]))[1]/n
  }))
}
load("/Volumes/HDD/LRdb.RData")
Lig <- c("Myeloid",
  "Epithelial_Basal",
   "Epithelial_Basal",
  "Endothelial",
  
  "Epithelial_Basal",
  "Myoepithelial",
  "Epithelial_Basal",
  "Epithelial_Luminal_Mature",
  
  "Myeloid",
  "iCAFs",
  "Myeloid",
  "myCAFs",
  
  "Myeloid",
  "CD4+ T-cells",
  "Myeloid",
  "T-Regs",
  
  "Epithelial_Luminal_Mature",
  "Myoepithelial"
  )
Recp <- c("Epithelial_Basal",
          "Myeloid",
          "Endothelial",
          "Epithelial_Basal",
          
          "Myoepithelial",
          "Epithelial_Basal",
          "Epithelial_Luminal_Mature",
          "Epithelial_Basal",
          
          "iCAFs",
          "Myeloid",
          "myCAFs",
          "Myeloid",
          
          "CD4+ T-cells",
          "Myeloid",
          "T-Regs",
          "Myeloid",
          
          "Myoepithelial",
          "Epithelial_Luminal_Mature"
          )
pr <- paste("pair", 1:18, sep = "")
pt <- paste("patient", 1:5, sep = "")
Rep <- 20
 

for (i in 1:length(pr)){
  Lig.cell <- Lig[i]
  Recp.cell <- Recp[i]

  for (j in 1:length(pt)){
    setwd(paste("/Volumes/Extreme SSD/all/folders/", pr[i], pt[j], sep = ""))
    
    load(paste(c(c(paste("Ligand (",Lig.cell,")", sep =""),
                   paste("Receptor (", Recp.cell,")", sep =""),"PMDvalue"),kk, ".RData"), collapse = "_"))
    
    if (!((Lig.n==0)|(Recp.n==0))) {
      pmd.l <- list()
      reps <- 1:Rep

      for (kk in reps){
        load(paste(c(c(paste("Ligand (",Lig.cell,")", sep =""),
                       paste("Receptor (", Recp.cell,")", sep =""),"PMDvalue"),kk, ".RData"), collapse = "_"))
        pmd.l[[kk]] <- pmd
      }

      pmd <- Reduce("+", pmd.l) / length(pmd.l)
      
      pmd.nam <- rownames(pmd)
      table1<-matrix("",lig.P*recp.P,3)
      k <- 1
      for (i1 in 1:lig.P){
        for (j1 in 1:recp.P){
          table1[k,]<-c(pmd.nam[i1],pmd.nam[lig.P + j1],
                        mean(pmd[i1,lig.P + j1],pmd[lig.P + j1, i1]) )
          k <- k + 1
        }
      }
      
      table1 <- data.frame(table1)
      colnames(table1) <- c(paste("Ligand (",Lig.cell,")", sep =""),
                            paste("Receptor (", Recp.cell,")", sep =""),"IMDI")
      
      table1[,1] <- gsub(".Lig", "", table1[,1])
      table1[,2] <- gsub(".Recp", "", table1[,2])
      
      table1 <- table1[which(table1[,3] != 1 ), ]
      
      table1clean <- c()
      for (ii in 1:nrow(table1)){
        if (max((LRdb$ligand %in% table1[ii, 1])+(LRdb$receptor %in% table1[ii,2])) == 2){
          table1clean <- rbind(table1clean,table1[ii,])
        }
      }
      table1clean <- table1clean[order(table1clean$IMDI, decreasing = FALSE),]
      
      
      setwd("/Volumes/Extreme SSD/all/all_ranks")
      nms <- c(paste("Ligand (",Lig.cell,")", sep =""),
               paste("Receptor (", Recp.cell,")", sep =""),"IMDI")
      nms <- paste(c(nms, paste("Patient", j, sep =""), ".csv"), collapse = "_")

      table1clean$Pvalue <- pvalues(pmd, x = as.numeric(pairs$IMDI))

      write.csv(table1clean, file = nms, row.names = FALSE)
    }}}