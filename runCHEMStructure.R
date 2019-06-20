source("CHEMAPI.R")

#I don't think I can call external to the server i.e., I need to have net access for this part
#drugFile <- loadCIDFile("drugs_names.txt")
compounds <- loadSDFFile("compounds.sdf") 

aromaticMatrix <- clusterWithFMCSAromatic(compounds)
saveClusterMatrix(aromaticMatrix,"aromaticMatrix.txt")

ringMatrix <- clusterWithFMCSRing(compounds)
saveClusterMatrix(ringMatrix,"ringMatrix.txt")

staticMatrix <- clusterWithFMCSStatic(compounds)
saveClusterMatrix(staticMatrix,"staticMatrix.txt")

stop("STOPPING")

print("Loading CID entries.")
cidDF <- loadCHEMDF("cidDF.rda")
cidVector <- cidDF$CID
print(paste("A total of",length(cidVector),"CID entries."))

db <- getBioassayDatabase("pubchem_protein_only.sqlite")
print("Building bioassay cluster matrix.")
clusterMatrix <- clusterCompoundsByActivityProfile(db, as.character(cidVector))

#print("Saving bioassay cluster matrix.")
#df <- as.data.frame(clusterMatrix)
#saveCHEMDF(df,"clusterMatrix.rda")
#saveClusterMatrix(clusterMatrix,"bioactivityCluster.txt")


disconnectBioassayDB(db)
