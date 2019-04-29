
getBioassayDatabase <- function(DBLocation) {
  if(!file.exists(DBLocation)) {
    stop(paste("The database file location", DBLocation, "does not exist."))
  }

  pubChemDB <- connectBioassayDB(DBLocation)
  print("Connecting to protein DB")
  print(pubChemDB)
  return(pubChemDB)
}

#===============================================================
#return a dataframe (within a list) showing the drug targets (protein), the number of total assay screens and the total fraction of activity
getProteinDrugTargets <- function(db, CID) {
  #returns a matrix
  drugTargets <- activeTargets(db, CID)
  if(is.logical(drugTargets)) {
    print("CID was not recognised. Returning NULL.")
    return(NULL)
  }
  print("Returning drug **protein** NCBI Protein identifiers Matrix.")
  print(paste("Number of protein entries", nrow(drugTargets)))
  g <- getNCBIFractionActivity(drugTargets)
  resultList <- list(drugTargets,g)
  return(resultList)
}

#===============================================================

getUniProtIDs <- function(NCBIMatrix, db) {
  #run translateTargetId on each target identifier
  uniProtIDs <- lapply(row.names(NCBIMatrix), translateTargetId, database=db, category="UniProt")
  #if any targets had more than one UniProt ID, keep only the first one. The annotation details we are looking for are all likely to the the same. 
  uniProtIDs <- sapply(uniProtIDs, function(x) x[1])
  return(uniProtIDs)
}

#===============================================================

getEnsemblProteinDetails <- function(uniProtIDs, attributeCharVector, filtersCharVector, drugTargets) {
  if(length(attributeCharVector) == 0) {
    print("Missing attribute character vector entries. Returning NULL.")
    return(NULL)
  }
  ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")
  proteinDetails <- getBM(attributes=attributeCharVector, filters=filtersCharVector, mart=ensembl, values=uniProtIDs)

  return(proteinDetails)
}

#===============================================================

clusterCompoundsByActivityProfile <- function(db, compoundCIDs) {
 
  #retrieve activity data for a given list of compounds into a bioassaySet object
  assayData <- getBioassaySetByCids(db, compoundCIDs)
  print(assayData)

  #convert assayData into a matrix of "targets(rows) vs. compounds(cols)"
  #data from multiple assays hitting the same target are summarized into a single value

  #in the matrix a (2) represents an active compound vs target combination, a (1) represents an inactive combindation, and a (0) represents an untested or inconclusive combination. Inactive results are filtered out with inactive=FALSE.
  activityMatrix <- perTargetMatrix(assayData, inactive=FALSE, summarizeReplicates = "mode")
  #print(activityMatrix)

  assayTargets <- assaySetTargets(assayData)
  customMerge <- sapply(assayTargets, translateTargetId, database=db, category="kClust")

  #the merged matrix is smaller because several protein targets have been collapsed into single clusters.
  mergedActivityMatrix <- perTargetMatrix(assayData, inactive=FALSE, assayTargets=customMerge)
  print(paste("Before merging targets by similar sequence:",dim(activityMatrix)))
  print(paste("After merging", dim(mergedActivityMatrix)))

  binaryMatrix <- 1*(mergedActivityMatrix > 1)

  transposedMatrix <- t(binaryMatrix)
  distanceMatrix <- dist(transposedMatrix)

  clusterResults <- hclust(distanceMatrix, method="average")
  return(clusterResults)

  #plot(clusterResults)
}


