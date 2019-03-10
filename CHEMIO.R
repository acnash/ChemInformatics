loadCIDFile <- function(fileNameString) {
  if(!file.exists(fileNameString)) {
    print(paste(fileNameString,"does note exist"))
    return(NULL)
  } else {
    asDF <- read.csv(fileNameString, sep="\t",header=F)
    names(asDF) <- c("Name","CID")
    if(nrow(asDF) == 0) {
      print("There are no values in the CID drug file. Returning NULL.")
      return(NULL)
    }
    asDF$CID <- as.numeric(asDF$CID)
    asDF$Name <- as.character(asDF$Name)
    result <- asDF$CID
    if(sum(is.na(result)) > 0) { 
      print("There are NA values in the CID list. Returning poorly formatted DF.") 
      return(asDF)
    }
    print(paste("Loading",length(result),"CID entries"))
    return(asDF)
  }
}

#=====================================================

importSDFFromCID <- function(CIDvector) {
  if(!is.numeric(CIDvector)) {
    print(paste("CIDvector is not numeric. Of type:",class(CIDvector),"Returning NULL"))
    return(NULL)
  }
  if(length(CIDvector) == 0) {
    print("There are no CID entries in the CIDvector. Returning NULL.")
    return(NULL)
  }

  concat <- NULL
  for(i in 1:length(CIDvector)) {
    compounds <- pubchemCidToSDF(CIDvector[i])
    if(length(concat) == 0) {
      concat <- compounds
    } else {
      concat <- c(concat,compounds)
    }
    print(paste("Retrieving",i,"SDF entry from CID",CIDvector[i]))
  }
  print("Making SDF IDs unique")
  uniqueIDs <- makeUnique(sdfid(concat))
  cid(concat) <- uniqueIDs
  return(concat)
}

#======================================================

saveCHEMDF <- function(chemDF, fileName) {
  saveRDS(chemDF, fileName)
}

#======================================================

loadCHEMDF <- function(fileName) {
  if(!file.exists(fileName)) {
    print(paste("File",fileName,"does not exist. Returning NULL."))
    return(NULL)
  }
  df<-readRDS(fileName)
  return(df)
}

#=======================================================

saveSDFFile <- function(SDFset, fileName) {
  if(!(endsWith(fileName,".sdf"))) {
    print(paste("The file",fileName,"does not end with an .sdf extension. Will not save."))
    return(NULL)
  }
  write.SDF(SDFset, file=fileName, cid=TRUE)
}

#=======================================================

loadSDFFile <- function(fileName) {
  if(!file.exists(fileName)) {
    print(paste("File",fileName,"does not exist. Returning NULL."))
    return(NULL)
  }
  sdf <- read.SDFset(sdfstr=fileName)
  return(sdf)
}

#=======================================================

loadClusterMatrix <- function(fileName) {
  if(!file.exists(fileName)) {
    print(paste("File",fileName,"does not exist. Returning NULL."))
    return(NULL)
  }
  simMatrix <- as.matrix(read.table(fileName,header=TRUE))
  names <- colnames(simMatrix)
  row.names(simMatrix) <- names
  return(simMatrix)
}

#=======================================================

saveClusterMatrix <- function(simMatrix, fileName) {
  write.matrix(simMatrix, fileName, sep="\t")
}

