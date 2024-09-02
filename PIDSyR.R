library(Informeasure)
library(data.table)
library(parallel)
library(dplyr)
library(pbapply)
library(pbmcapply)
library(stringr)

# Load initial miRNA and mRNA data and perform simple processing
cancers <- c("BRCA", "LIHC", "PRAD")
for(i in cancers)
{
  RBPmRNAv27 <-  read.table(file = paste0("inputs/RBPmRNAv27.bed"), header=F, check.name=F)
  RBP <- unique(RBPmRNAv27$V4)  # These are the RBPs to be extracted
  mRNAexpression <- read.table(file = paste0("inputs/",i,"_mRNA_FPKM.bed"), header=TRUE, check.name=F)
  miRNAexpression <- read.table(file = paste0("inputs/",i,"_miRNA_RPM.bed"), header=TRUE, check.name=F)
  rownames(mRNAexpression) <- mRNAexpression$gene_name
  mRNAexpression <- mRNAexpression[,-c(1, 2)]
  rownames(miRNAexpression) <- miRNAexpression$miRNA_id
  miRNAexpression <- miRNAexpression[,-c(1, 2)]
  
  # Handle column names: remove duplicates and keep the first 16 characters for each patient
  micol <- data.frame(cols = colnames(miRNAexpression))
  micol$cols1 <- substr(micol$cols, start = 1, stop = 16)
  micol$unique <- !duplicated(micol$cols1)
  mcol <- data.frame(cols = colnames(mRNAexpression))
  mcol$cols1 <- substr(mcol$cols, start = 1, stop = 16)
  mcol$unique <- !duplicated(mcol$cols1)
  
  # Extract the relevant columns
  miRNAexpression <- miRNAexpression[,micol$unique]
  mRNAexpression <- mRNAexpression[,mcol$unique]
  
  # Set column names
  colnames(miRNAexpression) <- substr(colnames(miRNAexpression), start = 1, stop = 16)
  colnames(mRNAexpression) <- substr(colnames(mRNAexpression), start = 1, stop = 16)
  
  # Find common column names between the two data frames
  common_colnames <- intersect(colnames(miRNAexpression), colnames(mRNAexpression))
  
  # Keep only the common column names in both data frames
  miRNAexpression <- miRNAexpression[, common_colnames]
  mRNAexpression <- mRNAexpression[, common_colnames]
  
  # Create RBP expression data
  RBPexpression <- mRNAexpression[RBP,]
  mRNAexpression <- mRNAexpression[!(rownames(mRNAexpression) %in% RBP), ]
  
  # Verify columns
  table(colnames(miRNAexpression) == colnames(mRNAexpression))
  
  # Save the data
  save(RBPexpression, miRNAexpression, mRNAexpression, file = paste0("outputs/",i,"_data.RData"))
  write.table(miRNAexpression, file = paste0("outputs/", i, "_miRNA_RPM.bed"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  write.table(mRNAexpression, file = paste0("outputs/", i, "_mRNA_FPKM.bed"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  write.table(RBPexpression, file = paste0("outputs/", i, "_RBP_FPKM.bed"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
}

rm(list = ls())

# Compute PID
cancers <- c("BRCA", "LIHC", "PRAD")
triplet <- read.table("./inputs/RBP_m_mi.bed", header=F)
for(i in cancers)
{
  load(file = paste0("outputs/",i,"_data.RData"))
  mRNAexpression  <- log2(mRNAexpression + 1)
  miRNAexpression   <- log2(miRNAexpression + 1)
  RBPexpression <- log2(RBPexpression + 1)
  
  ### Compute PID ####
  colnames(triplet)[1] <- "RBP"
  colnames(triplet)[2] <- "mRNA"
  colnames(triplet)[3] <- "miRNA"
  triplet <- triplet[,c(1:3)]
  
  # Extract row names of expression profiles
  miRNA_rownames <- rownames(miRNAexpression)
  RBP_rownames <- rownames(RBPexpression)
  mRNA_rownames <- rownames(mRNAexpression)
  
  # Filter the triplet data frame
  filtered_triplet <- triplet[
    triplet$miRNA %in% miRNA_rownames &
      triplet$RBP %in% RBP_rownames &
      triplet$mRNA %in% mRNA_rownames,]
  rownames(filtered_triplet) <- 1:nrow(filtered_triplet)
  
  # Function to compute PID
  implement.function = function(RBP, mRNA, miRNA, 
                                RBPexpression,
                                mRNAexpression,
                                miRNAexpression){
    
    miRNAindex <- which(rownames(miRNAexpression) == miRNA, arr.ind = F)
    RBPindex <- which(rownames(RBPexpression) == RBP, arr.ind = F)
    mRNAindex <- which(rownames(mRNAexpression) == mRNA, arr.ind = F)
    
    vector <- c()
    
    if(length(miRNAindex) & length(RBPindex) & length(mRNAindex)){
      
      xyz <- discretize3D(as.numeric(miRNAexpression[miRNAindex, ]),
                          as.numeric(mRNAexpression[mRNAindex, ]),
                          as.numeric(RBPexpression[RBPindex, ]))
      
      pid <- PID.measure(xyz)
      Synergy <- pid[1,1]
      Unique_X <- pid[1,2]
      Unique_Y <- pid[1,3]
      Redundancy <- pid[1,4]
      PID <- pid[1,5]
      
      vector <- c(RBP, miRNA, mRNA, Synergy, Unique_X, Unique_Y, Redundancy, PID)
      
    }
    else{
      vector <- NULL
    }
    
    return(vector)
  }
  
  # Parallel processing
  cores <- detectCores(logical = FALSE)
  cat("The number of cores is ", cores, "\n")
  pboptions(type = "tk")
  system.time(
    infoResults <- pbmclapply(1:(dim(filtered_triplet)[1]), 
                              FUN = function(x) { 
                                implement.function(RBP = filtered_triplet[x,1], mRNA = filtered_triplet[x,2], miRNA = filtered_triplet[x,3],
                                                   RBPexpression = RBPexpression,
                                                   mRNAexpression = mRNAexpression,
                                                   miRNAexpression = miRNAexpression)
                              }, 
                              mc.cores = cores
    )
  )
  
  # Combine results and save
  # RBP, miRNA, mRNA, Synergy, Unique_X, Unique_Y, Redundancy, PID
  infoResults <- do.call("rbind", infoResults) 
  colnames(infoResults) <- c("RBP","miRNA", "mRNA",  "Synergy", "Unique_X", "Unique_Y", "Redundancy", "PID")
  write.table(infoResults, paste0("outputs/triplets",i,"_PID.bed"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

# Compute Pearson correlations
rm(list = ls())
cancers <- c("BRCA", "LIHC", "PRAD")
triplet <- read.table("./inputs/RBP_m_mi.bed", header=F)
for(i in cancers)
{
  load(file = paste0("outputs/",i,"_data.RData"))
  mRNAexpression  <- log2(mRNAexpression + 1)
  miRNAexpression   <- log2(miRNAexpression + 1)
  RBPexpression <- log2(RBPexpression + 1)
  
  ### Compute Pearson correlations ####
  colnames(triplet)[1] <- "RBP"
  colnames(triplet)[2] <- "mRNA"
  colnames(triplet)[3] <- "miRNA"
  triplet <- triplet[,c(1:3)]
  
  # Extract row names of expression profiles
  miRNA_rownames <- rownames(miRNAexpression)
  RBP_rownames <- rownames(RBPexpression)
  mRNA_rownames <- rownames(mRNAexpression)
  
  # Filter the triplet data frame
  filtered_triplet <- triplet[
    triplet$miRNA %in% miRNA_rownames &
      triplet$RBP %in% RBP_rownames &
      triplet$mRNA %in% mRNA_rownames,]
  rownames(filtered_triplet) <- 1:nrow(filtered_triplet)
  
  # Directly compute Pearson correlations between matrices
  mi <- t(miRNAexpression)
  m <- t(mRNAexpression)
  rbp <- t(RBPexpression)
  
  # Calculate Pearson correlations between expression profiles
  mi_m_matrix <- cor(mi, m)
  mi_rbp_matrix <- cor(mi,rbp)
  rbp_m_matrix <- cor(rbp,m)
  
  implementp.function = function(RBP, mRNA, miRNA, 
                                 mi_m_matrix,
                                 mi_rbp_matrix,
                                 rbp_m_matrix){
    pearson1 <- mi_m_matrix[miRNA,mRNA]
    pearson2 <- mi_rbp_matrix[miRNA,RBP]
    pearson3 <- rbp_m_matrix[RBP,mRNA]
    vector <- c(RBP, miRNA, mRNA, pearson1,pearson2,pearson3)
    return(vector)
  }
  
  # Parallel processing
  cores <- detectCores(logical = FALSE)
  cat("The number of cores is ", cores, "\n")
  pboptions(type = "tk")
  system.time(
    infoResults <- pbmclapply(1:(dim(filtered_triplet)[1]), 
                              FUN = function(x) { 
                                implementp.function(RBP = filtered_triplet[x,1], mRNA = filtered_triplet[x,2], miRNA = filtered_triplet[x,3],
                                                    mi_m_matrix,
                                                    mi_rbp_matrix,
                                                    rbp_m_matrix)
                              }, 
                              mc.cores = cores
    )
  )
  
  # Combine results and save
  # RBP, miRNA, mRNA, Synergy, Unique_X, Unique_Y, Redundancy, PID
  infoResults <- do.call("rbind", infoResults) 
  colnames(infoResults) <- c("RBP","miRNA", "mRNA", "pearson_mi_m", "pearson_mi_RBP", "pearson_RBP_m")
  write.table(infoResults, file = paste0("outputs/triplets",i,"_pearson.bed"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  dim(infoResults)
}

#### Merge results ####
rm(list = ls())
cancers <- c("BRCA", "LIHC", "PRAD")
triplet <- read.table("./inputs/RBP_m_mi.bed", header=F)
colnames(triplet)[1] <- "RBP"
colnames(triplet)[2] <- "mRNA"
colnames(triplet)[3] <- "miRNA"
triplet <- triplet[,c(1:3)]
for(i in cancers)
{
  PID <- read.table(file = paste0("outputs/triplets",i,"_PID.bed"),header=TRUE, check.name=F)
  person <-  read.table(file = paste0("outputs/triplets",i,"_pearson.bed"),header=TRUE, check.name=F)
  triplet_cancer <- merge(PID, person, by = c("RBP", "mRNA", "miRNA"))
  write.csv(triplet_cancer, file = paste0("outputs/triplets",i,".csv"), row.names = FALSE)
}
