
clusterWithFMCSAromatic <- function(SDFset, au=2, bu=1, overlapCoefficient=TRUE) {
  if(overlapCoefficient==TRUE) {
    print("Running batch aromatic FMCS and returning the overlap coefficient")
    d <- sapply(cid(SDFset), function(x) fmcsBatch(SDFset[x], SDFset, au=au, bu=bu, matching.mode="aromatic", numParallel=12)[,"Overlap_Coefficient"])
  } else {
    print("Running batch aromatic FMCS and returning the Tanimoto coefficient")
    d <- sapply(cid(SDFset), function(x) fmcsBatch(SDFset[x], SDFset, au=au, bu=bu, matching.mode="aromatic", numParallel=12)[,"Tanimoto_Coefficient"])
  }
  return(d) 
}

#=========================================

clusterWithFMCSStatic <- function(SDFset, au=2, bu=1, overlapCoefficient=TRUE) {
  if(overlapCoefficient==TRUE) {
    print("Running batch static FMCS and returning the overlap coefficient")
    d <- sapply(cid(SDFset), function(x) fmcsBatch(SDFset[x], SDFset, au=au, bu=bu, matching.mode="static", numParallel=12)[,"Overlap_Coefficient"])
  } else {
     print("Running static aromatic FMCS and returning the Tanimoto coefficient")
    d <- sapply(cid(SDFset), function(x) fmcsBatch(SDFset[x], SDFset, au=au, bu=bu, matching.mode="static", numParallel=12)[,"Tanimoto_Coefficient"])
  }
  return(d)
}

#=========================================

clusterWithFMCSRing <- function(SDFset, au=2, bu=1, overlapCoefficient=TRUE) {
   if(overlapCoefficient==TRUE) {
    print("Running batch ring FMCS and returning the overlap coefficient")
    d <- sapply(cid(SDFset), function(x) fmcsBatch(SDFset[x], SDFset, au=au, bu=bu, matching.mode="ring", numParallel=12)[,"Overlap_Coefficient"])
  } else {
     print("Running batch ring FMCS and returning the Tanimoto coefficient")
    d <- sapply(cid(SDFset), function(x) fmcsBatch(SDFset[x], SDFset, au=au, bu=bu, matching.mode="ring", numParallel=12)[,"Tanimoto_Coefficient"])
  }
  return(d)
}

#=========================================

calculateFMCS <- function(sdfObject1, sdfObject2, au=2, bu=1) {
  mcsa <- fmcs(sdfObject1, sdfObject2, fast=FALSE, au=au, bu=bu)
  print("Running Flexible Maximum Common Substruction (FMCS) calculation on two SDF objects")
  print("execute plotMCS(msca, regenCoords=TRUE) to see a visual comparison")
  print(mcsa)
  return(mcsa)
}

#=========================================

calculateNumAtomsToMCSReference <- function(referenceSDF, ignoreIndex, SDFset, drugNameVector, au=2, bu=1) {
  counter <- 1
  drugsVector <- as.character()
  sizeVector <- as.integer()
  for(i in 1:length(SDFset)) {
    if(i == ignoreIndex) {
      next
    }
    mcsa <- calculateFMCS(referenceSDF, SDFset[[i]], au, bu)
    #MCS size in atoms
    size <- mcsa@stats[3]
    sizeVector[counter] <- size
    drugName <- drugNameVector[i]
    drugsVector[counter] <- drugName
    counter <- counter + 1
  }

  df <- data.frame(drug=drugsVector, size=sizeVector, stringsAsFactors=TRUE)
  return(df)
}
