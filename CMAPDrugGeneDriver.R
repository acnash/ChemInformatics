getDrugGeneExpCMAP <- function(geneUpFile, geneDownFile, gctDF, colNumber) {
  geneListUp <- read.csv(geneUpFile,header = F,col.names = c("gene","ExpValue"))
  geneListDown <- read.csv(geneDownFile,header = F,col.names = c("gene","ExpValue"))
  geneListUp$gene <- tolower(geneListUp$gene)
  geneListDown$gene <- tolower(geneListDown$gene)
  
  cell <- dplyr::select(gctDF,colNumber)
  cellName <- as.character(cell["cell_id",])
  dose <- as.character(cell["pert_idose",])
  time <- as.character(cell["pert_itime",])
  
  cell <- cell[grepl("^[0-9]+.*", rownames(cell)),]
  namesToUse <- names(cell)
  
  cellDF <- data.frame(gene=namesToUse, value=as.character(cell),stringsAsFactors = FALSE)
  
  geneListUpFromDrug <- cellDF[cellDF$value >= 0,]
  geneListDownFromDrug <- cellDF[cellDF$value < 0,]
  geneListUpFromDrug <- geneListUpFromDrug[base::order(geneListUpFromDrug$value,decreasing = TRUE),]
  geneListDownFromDrug <- geneListDownFromDrug[base::order(geneListDownFromDrug$value,decreasing = FALSE),]
  
  
  #WHAT IS THIS and where the fuck did I get it from??????
  #I need to replace these with a entrezgene and gene (sumbol) column
  ###geneIDsDrugUp <- loadCPRDDataframe("~/Documents/ALEX/geneIDsDrugUp.Rda")
  ###geneIDsDrugDown <- loadCPRDDataframe("~/Documents/ALEX/geneIDsDrugDown.Rda")
  
  #using the entrezID from the "gene" column of geneListUpFromDrug, retrieve the true "gene" symbol
  #and store as entrez and gene in geneIDsDrugUp
  geneIDsDrugUp<-getEntrezGeneSymbolDF(geneListUpFromDrug)
  geneIDsDrugDown<-getEntrezGeneSymbolDF(geneListDownFromDrug)
  geneIDsDrugUp$gene_symbol <- tolower(geneIDsDrugUp$gene_symbol)
  geneIDsDrugDown$gene_symbol <- tolower(geneIDsDrugDown$gene_symbol)
  
  names(geneIDsDrugUp) <- c("entrezgene","gene")
  names(geneIDsDrugDown) <- c("entrezgene","gene")
  names(geneListUpFromDrug) <-c("entrezgene","value")
  names(geneListDownFromDrug) <-c("entrezgene","value")
  
  mergedRNAGenes <- rbind(geneListUp,geneListDown) #the UP and DOWN genes merged from our RNA seq exp
  
  
  
  mergedUpGenes <- merge(geneIDsDrugUp,geneListUpFromDrug,by="entrezgene",all.x=TRUE)
  mergedUpGenes <- merge(mergedUpGenes,mergedRNAGenes,by="gene",all.x=TRUE)
  mergedUpGenes <- mergedUpGenes[base::order(mergedUpGenes$ExpValue,decreasing = TRUE),] 
  
  mergedDownGenes <- merge(geneIDsDrugDown,geneListDownFromDrug,by="entrezgene",all.x=TRUE)
  mergedDownGenes <- merge(mergedDownGenes,mergedRNAGenes,by="gene",all.x=TRUE)
  mergedDownGenes <- mergedDownGenes[base::order(mergedDownGenes$ExpValue,decreasing = TRUE),] 
  
  #remove all genes that I do not have our up/down genes for (not the drug ones)
  mergedDownGenes <- mergedDownGenes[complete.cases(mergedDownGenes),]
  mergedUpGenes <- mergedUpGenes[complete.cases(mergedUpGenes),]
  
  if(length(intersect(mergedDownGenes$gene,mergedUpGenes$gene)) > 0) {
    print("A drug at a particular time and with a particular dose should not have the same gene in up and down space.")
    stop()
  }
  
  #mergedUpGenes as UP genes in the DRUG
  #mergedDownGenes are DOWN genes in the DRUG
  resultList <- list(cellName, mergedUpGenes, mergedDownGenes, dose, time)
  return(resultList)
}

#==============================

getEntrezGeneSymbolDF <- function(entrezDF) {
  entrezDF$gene_symbol <- mapIds(keys=entrezDF$gene,org.Hs.eg.db,column="SYMBOL",keytype="ENTREZID",multiVals="first")
  entrezDF <- dplyr::select(entrezDF,c(gene,gene_symbol))
  return(entrezDF)
}

#============================
#geneDF must look like:

#gene entrezgene     value ExpValue
#1      A375          0 6 h|10 µM    0.000
#1282  fosl1       8061     -8.89    2.090
#2931   rela       5970     -1.61    1.050
#3455    src       6714     -0.32    0.810
#2280 nfkbia       4792     -6.11    0.605

#This is generated from the function getDrugGeneExpCMAP (above in this file)
calculateSimilarityMap <- function(geneDF) {
  geneNames <- as.character(unique(geneDF$gene)) #this becomes row names and number of rows
  cellTimeDoseDF <- subset(geneDF, grepl(".*h.*",geneDF$value)==TRUE)
  
  possibleTimeDoses <- (as.character(unique(cellTimeDoseDF$value)))
  splitTimeDoses <- unlist(strsplit(possibleTimeDoses,"\\|"))
  
  possibleCellTypes <- (as.character(unique(cellTimeDoseDF$gene)))
  possibleDoses <- unique(splitTimeDoses[!grepl(".*h",splitTimeDoses)])  
  possibleTime <- unique(splitTimeDoses[grepl(".*h",splitTimeDoses)])  
  
  geneNames <- subset(geneNames, !(geneNames %in% possibleCellTypes))
  
  numMatrices <- length(possibleDoses)
  numCols <- length(possibleTime)*length(possibleCellTypes)
  numRows <- length(geneNames)
  
  resultList <- list()
  
  #you make one matrix at a time, each time you parse through the whole dataframe looking for matches
  for(n in 1:numMatrices) {
    colNames <- character(0)
    doseToConsider <- possibleDoses[n]
    count <- 1
    for(i in 1:length(possibleTime)) {
      for(j in 1:length(possibleCellTypes)) {
        colNames[count] <- paste0(possibleCellTypes[j],"|",possibleTime[i],"|",possibleDoses[n])
        count <- count + 1
      }
    }
  
    mat <- matrix(data = 0, nrow = numRows, ncol = numCols)
    row.names(mat) <- geneNames
    colnames(mat) <- colNames
    
    found <- FALSE
    for(i in 1:nrow(geneDF)) {
      entry<-geneDF[i,]
      
      if(entry$gene %in% possibleCellTypes) {
        timeDoseData <- as.character(entry$value[1])
        timeFound <- unlist(strsplit(timeDoseData,"\\|"))[1]
        doseFound <- unlist(strsplit(timeDoseData,"\\|"))[2]
        cellTypeToUse <- entry$gene[1]
        if(doseFound == doseToConsider) { 
          found <- TRUE
          next()
        } else {
          found <- FALSE
          next()
        }
      }
      
      if(found == TRUE) {
        #use the entry
        
        gene <- entry$gene[1]
        #print(class(entry))
        #print(str(entry))
        #print("=====")
        #print(entry$value)
        #print(entry$value[1])
        #print(as.numeric(as.character(entry$value[1])))
        #stop()
        
        value <- as.numeric(as.character(entry$value[1]))
        expValue <- as.numeric(entry$ExpValue[1])
        
        diffValue <- 0
        if(value >= 0 & expValue >=0) {
          diffValue <- -1
        } else if(value >= 0 & expValue < 0) {
          diffValue <- 1
        } else if(value < 0 & expValue >= 0) {
          diffValue <- 1
        } else {
          diffValue <- -1
        }
        
        #print(paste(value,expValue,diffValue))
        
        colEntry <- paste0(cellTypeToUse,"|",timeFound,"|",doseFound)
        rowEntry <- as.character(gene)
        mat[c(rowEntry),c(colEntry)] <- diffValue
        
      }
    } #end of geneDF loop
    resultList[[n]] <- mat
  } #end of matrix loop
  
  return(resultList)
}


#=======================================================
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))








