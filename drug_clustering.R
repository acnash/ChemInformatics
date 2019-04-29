library(cmapR)
library(httr)
library(jsonlite)
library(ggplot2)
library(dendextend)
library("energy")
setwd("/Users/anthonynash/Documents/CMAP")
source("../CHEM/CHEMIO.R")

L1000FWD_URL <- 'http://amp.pharm.mssm.edu/L1000FWD/'

#
calculateEuclidDistanceDiffGenes <- function(completeGeneVector, diffGenesAUP, diffGenesBUP, diffGenesADOWN, diffGenesBDOWN) {
  compareVectorA <- as.numeric(rep(0,length=length(completeGeneVector)))
  names(compareVectorA) <- completeGeneVector
  compareVectorB <- compareVectorA
  
  foundUpA <- 0
  for(gene in diffGenesAUP) {
    compareVectorA[names(compareVectorA)==gene] <- 1
    foundUpA <- foundUpA + 1
  }
  stopifnot(foundUpA == length(diffGenesAUP))

  foundUpB <- 0
  for(gene in diffGenesBUP) {
    compareVectorB[names(compareVectorB)==gene] <- 1
    foundUpB <- foundUpB + 1
  }
  stopifnot(foundUpB == length(diffGenesBUP))
  
  foundDownA <- 0
  for(gene in diffGenesADOWN) {
    compareVectorA[names(compareVectorA)==gene] <- -1
    foundDownA <- foundDownA + 1
  }
  stopifnot(foundDownA == length(diffGenesADOWN))
  
  foundDownB <- 0
  for(gene in diffGenesBDOWN) {
    compareVectorB[names(compareVectorB)==gene] <- -1
    foundDownB <- foundDownB + 1
  }
  stopifnot(foundDownB == length(diffGenesBDOWN))
  
  
  euclidDistance<-dist(rbind(compareVectorA, compareVectorB))
  
  correlationDistance <- dcor(compareVectorA, compareVectorB)
  
  resultList <- list(euclidDistance, correlationDistance)
  return(resultList)
}





#load in the genes
geneTable <- read.table("available_genes.txt",header=TRUE, sep="\n", strip.white = TRUE, stringsAsFactors = FALSE)
availableGeneVector <- geneTable$Symbol

#load in Drugs
drugTable <- read.table("drugs.txt",header=FALSE, sep = "\n")
colnames(drugTable) <- c("DrugName")

noMissing <- 0
noFound <- 0
responseList <- list()
drugNamesFoundList <- list()
for(i in 1:nrow(drugTable)) {
  drug <- as.character(drugTable$DrugName[i])
  drugStartName <- strsplit(drug,perl = TRUE, split = " ")
  response <- GET(paste0(L1000FWD_URL, 'synonyms/', drugStartName[[1]][1]))
  if (response$status_code == 200){
    noFound <- noFound + 1
    response <- fromJSON(httr::content(response, 'text'))
    print(paste("Found match for:",drugStartName[[1]][1]))
    responseList[[noFound]] <- response
    drugNamesFoundList[[noFound]] <- drugStartName[[1]][1]
  } else {
    print(paste("Missing:",drugStartName[[1]][1]))
    noMissing <- noMissing + 1
  }
}

#This code pulls out all the pert_id from a list of drug names
print("================================================")
print(paste("Found",noFound,"potential drugs"))
print(paste("Missing",noMissing,"drugs"))

drugDF <- data.frame(Name=character(), pert_id=character(), stringsAsFactors = FALSE) #Name & #pert_id
for(i in 1:length(responseList)) {
  similarDrugNames <- responseList[i]
  
  entry <- similarDrugNames[[1]]
  if(length(entry)==0) {
    next()
  }
  
  lookingFor <- drugNamesFoundList[[i]]
  if(nrow(entry) == 1) {
    tempDF <- data.frame(Name=similarDrugNames[[1]][1], pert_id=similarDrugNames[[1]][2], stringsAsFactors = FALSE)
    drugDF <- rbind(drugDF,tempDF)
  } else {
    for(j in 1:nrow(entry)) {
      print(entry$Name[j])
      if(entry$Name[j]==lookingFor) {
        tempDF <- data.frame(Name=entry$Name[j], pert_id=entry$pert_id[j], stringsAsFactors = FALSE)
        drugDF <- rbind(drugDF, tempDF)
      }
    }
  }
}



#=============================================================================
#dsPath <- "GSE70138_Broad_LINCS_Level3_INF_mlr12k_n78980x22268_2015-06-30.gct"
#gctxPath <- "GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx"
#takes up far too much RAM - DO NOT RUN
#myDS <- parse.gctx(dsPath)

#col_meta_path <- "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt"
#col_meta <- read.delim(col_meta_path, sep="\t", stringsAsFactors=F)

# figure out which signatures correspond to vorinostat by searching the 'pert_iname' column
#idx <- which(col_meta$pert_iname=="vorinostat")

# and get the corresponding sig_ids
#sig_ids <- col_meta$sig_id[idx]

# read only those columns from the GCTX file by using the 'cid' parameter
#vorinostat_ds <- parse.gctx(gctxPath, cid=sig_ids)

#cid - vector or a path to file containing column indices
#rid - vector or a path to file containing row indices
#=============================================================================

#fullMatrixGCTXPath <- "CD_signatures_full_42809x22268.gctx"
fullMatrixGCTXPath <- "CD_signatures_LM_42809x978.gctx"
drugMetadataPath <- "Drugs_metadata.csv"
signatureMetadataPath <- "CD_signature_metadata.csv"

drugMetadata <- read.delim(drugMetadataPath, sep="\t", stringsAsFactors=F)
signatureMetadata <- read.delim(signatureMetadataPath, sep="\t", stringsAsFactors=F)
#I can pull out the up/down genes from the drug for a particular cell type using the signatureMetadata

#This only needs to be run once
n<-dim(signatureMetadata)[1]
signatureMetadataDF <- data.frame(sig_id=character(n),SCS_centered_by_batch=character(n), batch=character(n), cell_id=character(n), mean_cosine_dist_centered_by_batch=character(n), pert_desc=character(n), pert_dose=character(n), pert_id=character(n), pert_time=character(n), stringsAsFactors = FALSE)
for(i in 1:nrow(signatureMetadata)) {
  entries <- strsplit(signatureMetadata[i,1],split=",")
  signatureMetadataDF$sig_id[i] <- entries[[1]][1]
  signatureMetadataDF$SCS_centered_by_batch[i] <- entries[[1]][2]
  signatureMetadataDF$batch[i] <- entries[[1]][3]
  signatureMetadataDF$cell_id[i] <- entries[[1]][4]
  signatureMetadataDF$mean_cosine_dist_centered_by_batch[i] <- entries[[1]][5]
  signatureMetadataDF$pert_desc[i] <- entries[[1]][6]
  signatureMetadataDF$pert_dose[i] <- entries[[1]][7]
  signatureMetadataDF$pert_id[i] <- entries[[1]][8]
  signatureMetadataDF$pert_time[i] <- entries[[1]][9]
}

#entries[[1]][1]

#Go through each pertid
validDrugNameList <- list()
upgenesList <- list()
downgeneList <- list()
timeList <- list()
doseList <- list()
validCount <- 1
cellTypeList <- list()
mean_cosine_dist_centered_by_batchList <- list()
batchList <- list()

for(i in 1:nrow(drugDF)) {
  name <- drugDF$Name[i]
  pertID <- drugDF$pert_id[i]
  
  pertIDSubset <- signatureMetadataDF[signatureMetadataDF$pert_id == pertID,]
  if(nrow(pertIDSubset) == 0) {
    print("Not found")
    next()
  } 
  sigIDVector <- pertIDSubset$sig_id
  errorFound <- FALSE
  for(j in 1:length(sigIDVector)) {
    sig_id <- sigIDVector[j]
    print(sig_id)
    
      response <- GET(paste0(L1000FWD_URL, 'sig/', sig_id))
      if(response$status_code == 200) {
        tryCatch(
          response <- fromJSON(httr::content(response, 'text')),
          warning=function(x) {errorFound<-TRUE},
          error=function(x) {errorFound<-TRUE}
        )
        if(errorFound == TRUE) {
          print("Exception caught")
          errorFound <- FALSE
          next()
        }
        upgenesList[[validCount]] <- response$up_genes
        downgeneList[[validCount]] <- response$down_genes
        cellTypeList[[validCount]] <- response$cell_id
        timeList[[validCount]] <- response$pert_time
        doseList[[validCount]] <- response$pert_dose
        mean_cosine_dist_centered_by_batchList[[validCount]] <- response$mean_cosine_dist_centered_by_batch
        validDrugNameList[[validCount]] <- name
        batchList[[validCount]] <- response$batch
        validCount <- validCount + 1
      } 
  }
}

uniqueCellTypes <- unique(cellTypeList)

xVector <- integer()
yVector <- integer()
xDrug <- character()
yDrug <- character()
timeVectorX <- integer()
doseVectorX <- integer()
timeVectorY <- integer()
doseVectorY <- integer()
cellVector <- character()
#diffUpVector <- integer()
#diffDownVector <- integer()
diffDistance <- integer()
diffDistanceCDis <- integer()
#mean_cosine_dist_centered_by_batchListVectorX <- integer()
#mean_cosine_dist_centered_by_batchListVectorY <- integer()
#batchVectorX <- character()
#batchVectorY <- character()

currentStore <- 1
posX <- 1
posY <- 1
#clusterDF <- data.frame(x=integer(),y=integer(),cell=character(),diffGenesScoreUp<-integer(),diffGenesScoreDown<-integer(),stringsAsFactors = FALSE)

#this singleDF will contain all cell types
#loop through all the cell types 
#this loop with inner loop calculate the diff difference between drugs across multiple cell types
for(i in 1:length(uniqueCellTypes)) {
  #clusterMatrix <- matrix(nrow=,ncol=)
  currentCell <- uniqueCellTypes[i]
  print(currentCell[[1]])
  if(is.null(currentCell[[1]])) {
    next()
  }
  
  for(x in 1:length(validDrugNameList)) {
    cellToExamineX <- cellTypeList[x]
    if(is.null(cellToExamineX[[1]])) {
      next
    }
    if(!(currentCell[[1]] == cellToExamineX[[1]])) {
      next()
    }
    currentUpX <- upgenesList[[x]]
    currentDownX <- downgeneList[[x]]
    
    for(y in 1:length(validDrugNameList)) {
      cellToExamineY <- cellTypeList[y]
      if(is.null(cellToExamineY[[1]])) {
        next
      }
      if(!(currentCell[[1]] == cellToExamineY[[1]])) {
        next()
      }
      currentUpY <- upgenesList[[y]]
      currentDownY <- downgeneList[[y]]
      compareUpVector <- setdiff(currentUpX, currentUpY)
      compareDownVector <- setdiff(currentDownX, currentDownY)
      
      xVector[currentStore] <- posX
      yVector[currentStore] <- posY
      xDrug[currentStore] <- validDrugNameList[[x]]
      yDrug[currentStore] <- validDrugNameList[[y]]
      cellVector[currentStore] <- currentCell[[1]]
      
      difference <- calculateEuclidDistanceDiffGenes(availableGeneVector, currentUpX, currentUpY, currentDownX, currentDownY)
      diffDistance[currentStore] <- difference[[1]]
      diffDistanceCDis[currentStore] <- difference[[2]]
      #diffUpVector[currentStore] <- length(compareUpVector)/974  #this is not right
      #diffDownVector[currentStore] <- length(compareDownVector)/974  #this is not right
      timeVectorX[currentStore] <- timeList[[x]]
      doseVectorX[currentStore] <- doseList[[x]]
      timeVectorY[currentStore] <- timeList[[y]]
      doseVectorY[currentStore] <- doseList[[y]]
      #mean_cosine_dist_centered_by_batchListVectorX[[currentStore]] <- mean_cosine_dist_centered_by_batchList[[x]]
      #mean_cosine_dist_centered_by_batchListVectorY[[currentStore]] <- mean_cosine_dist_centered_by_batchList[[y]]
      #batchVectorX[[currentStore]] <- batchList[[x]]
      #batchVectorY[[currentStore]] <- batchList[[y]]
      
      
      #clusterTempDF <- data.frame(x=posX,y=posY,cell=currentCell[[1]],diffGenesScoreUp<-length(compareUpVector),diffGenesScoreDown<-length(compareDownVector),stringsAsFactors = FALSE)
      #clusterDF <- rbind(clusterDF,clusterTempDF)
      currentStore <- currentStore + 1
      posY<-posY+1
    }
    posX<-posX+1
  }
}


#The higher the score the fewer genes they share i.e., the score highlights the different genes
clusterDF <- data.frame(x=xVector,
                        y=yVector,
                        xDrug=xDrug,
                        yDrug=yDrug,
                        cell=cellVector,
                        timeX=timeVectorX,
                        timeY=timeVectorY,
                        doseX=doseVectorX,
                        doseY=doseVectorY,
                        #batchX=batchVectorX,
                        #batchY=batchVectorY,
                        #meanDistX=mean_cosine_dist_centered_by_batchListVectorX,
                        #meanDistY=mean_cosine_dist_centered_by_batchListVectorY,
                        distanceScore = diffDistance,
                        diffDistanceCDis = diffDistanceCDis,
                        #diffGenesScoreUp=diffUpVector,
                        #diffGenesScoreDown=diffDownVector,
                        stringsAsFactors = FALSE)

#refine the data by restricting to matching drug dosages and times.
dose <- 10
time <- 24
doseModDF <- clusterDF[clusterDF$doseX==dose & clusterDF$doseY==dose,]
doseTimeModDF <- doseModDF[doseModDF$timeX==time & doseModDF$timeY==time,]

#Go through each cell type to generate a matrix
cellResultList <- list()
for(i in 1:length(uniqueCellTypes)) {
  print(uniqueCellTypes[[i]])
  cellDF <- subset(doseTimeModDF, doseTimeModDF$cell == uniqueCellTypes[[i]])
  numSize <- dim(cellDF)[1]
  xName <- character(numSize)
  yName <- character(numSize)
  dist <- numeric(numSize)
  distCor <- numeric(numSize)
  combinedName <- character(numSize)
  
  for(j in 1:numSize) {
    xName[j] <- cellDF$xDrug[j]
    yName[j] <- cellDF$yDrug[j]
    dist[j] <- cellDF$distanceScore[j]
    distCor[j] <- cellDF$diffDistanceCDis[j]
    #dist[j] <- sqrt((cellDF$diffGenesScoreDown[j]^2) + (cellDF$diffGenesScoreUp[j]^2))
    combinedName[j] <- paste(cellDF$xDrug[j],"-",cellDF$yDrug[j])
  }
  
  tempDF <- data.frame(combinedName=combinedName,
                       xName=xName,
                       yName=yName,
                       dist=dist,
                       distCor=distCor,
                       stringsAsFactors=FALSE)
  
  cellResultList[[i]] <- tempDF
}

anchorDrugs <- c("CANDERSARTAN",
                 "LISINOPRIL",
                 "ATENOLOL",
                 "NADOLOL",
                 "VENLAFAXINE",
                 "AMITRIPTYLINE",
                 "PROPRANOLOL",
                 "METOPROLOL",
                 "GABAPENTIN",
                 "TOPIRAMATE",
                 "TIMOLOL",
                 "METHYSERGIDE",
                 "DIVALPROEX"
                 )

for(cellIndex in 1:length(uniqueCellTypes)) {
  currentCellDF <- cellResultList[[cellIndex]]
  for(i in 1:length(anchorDrugs)) {
  #for(i in 1:length(drugNamesFoundList)) {
  
    anchorDF <- subset(currentCellDF, currentCellDF$xName==anchorDrugs[i])
    #anchorDF <- subset(currentCellDF, currentCellDF$xName==drugNamesFoundList[[i]])
    if(nrow(anchorDF)==0) {
      next()
    }
    tempanchorDF <- anchorDF
  
    print(paste("Building plot for",anchorDrugs[i]))
  
    g <- ggplot(data=tempanchorDF, aes(x=tempanchorDF$yName, y=tempanchorDF$distCor)) +
    geom_segment( aes(x=tempanchorDF$yName, xend=tempanchorDF$yName, y=0, yend=tempanchorDF$distCor)) +
    geom_point( size=2, color="red", fill=alpha("blue", 0.2), alpha=0.7, shape=21, stroke=1) + 
    theme(legend.position = "none") +
    theme(axis.text.y=element_text(size=5)) + 
    theme(axis.text.x=element_text(size=5)) + 
    theme_classic() + xlab("Drug") + 
    ylab("Euclidean distance from anchor drug") + 
    ggtitle(paste(anchorDrugs[i],"10 mmol, 24 hours,",uniqueCellTypes[[cellIndex]])) + 
    theme(axis.text.x=element_text(angle=90,hjust=1)) +
    theme(axis.text.y=element_text(angle=45,hjust=1)) 
  
    fileName <- paste0("IMAGES_ANCHORS/",anchorDrugs[i],"10 mmol,", time, "hours,",uniqueCellTypes[[cellIndex]],".jpeg")
    print(fileName)
    jpeg(fileName, width=3000,height=2000,res=300)
    plot(g)
    dev.off()
  }
}

#Build a heatmap to represent a similarity matrix for the pearson correlation distances
#1 x running along the top that will require subsetting
#produce a matrix
#tempD <- NULL
for(n in 1:length(cellResultList)) {
  tempDF <- cellResultList[[n]]
  fileName<- paste0("IMAGES_ANCHORS/",uniqueCellTypes[[n]],"_",time,"_",dose,".jpeg")
  title <- paste0("Cell: ",uniqueCellTypes[[n]],". ","Time: ",time,". ","Dose: ",dose)
  
  compareDrugNames <- unlist(drugNamesFoundList)
  simMatrix <- matrix(data = 0, nrow = length(unique(tempDF$xName)), ncol = length(unique(tempDF$yName)))
  colnames(simMatrix) <- unique(tempDF$yName)
  row.names(simMatrix) <- unique(tempDF$xName)

  for(i in 1:length(compareDrugNames)) {
    currentAnchor <- compareDrugNames[i]
    if(sum(tempDF$xName %in% currentAnchor)== 0) {
      next()
    }
    indDF <- subset(tempDF, tempDF$xName == currentAnchor)
    print(currentAnchor)
    if(nrow(indDF)==0) {
      simMatrix[currentAnchor,] <- 0
      #set everything to zero
    } else {
      for(j in 1:nrow(indDF)) {
        yDrug <- indDF$yName[j]
        value <- indDF$distCor[j]
        simMatrix[currentAnchor,yDrug] <- value
      }
    }
  }
  d<-dist(simMatrix,method="manhattan")
  if(length(d) < 2) {
    next()
  }
  #tempD <- d
  #h<-hclust(d)
  dend <- d %>% hclust %>% as.dendrogram(method = "average")
  dend <- dend %>% set("labels_cex",0.9)
  dend <- dend %>% set("branches_k_color", k = 8)
  jpeg(fileName,width=4000,height=2000,res=300)
  errorFound <- FALSE
  tryCatch(
    #plot(h,main=title, xlab="Drug", ylab="Manhattan distance"),
    plot(dend,main=title),#, ylab="Manhattan distance"),
    warning=function(x) {errorFound<-TRUE},
    error=function(x) {errorFound<-TRUE}
  )
  if(errorFound == TRUE) {
    print("Exception caught")
    errorFound <- FALSE
    next()
  }
  dev.off()
}

#dend <- tempD %>% hclust %>% as.dendrogram(method = "average")

#d<-dist(simMatrix,method = "manhattan")
#h<-hclust(d)
#heat<-heatmap(simMatrix)

#uniqueDrugNames <- currentCellDF$xName
#propranololDF <- subset(currentCellDF, currentCellDF$xName=="PROPRANOLOL")
#dim(propranololDF)

#g<- ggplot(propranololDF, aes(x=propranololDF$yName, y=propranololDF$dist)) +
#  geom_segment( aes(x=propranololDF$yName, xend=propranololDF$yName, y=0, yend=propranololDF$dist)) +
#  geom_point( size=5, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=2)

#g <- g + xlab("Drug") + 
#  ylab("Diff gene difference") + 
#  theme(axis.text.x=element_text(angle=45,hjust=1)) +
#  theme(axis.text.y=element_text(angle=45,hjust=1)) +
#  theme(legend.position = "none") +
#  theme(axis.text.y=element_text(size=5)) + 
#  theme(axis.text.x=element_text(size=5)) + 
#  ggtitle(paste("Propranolol 10 mmol, 6 hours, A375"))




