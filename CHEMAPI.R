library(ChemmineR)
library(fmcsR)
library(dendextend)
library(RFLPtools)
library(viridis)
library(gplots)
library(MASS)
library(ggplot2)
library(bioassayR)
library(biomaRt)

source("CHEMIO.R")
source("CHEMAnalysis.R")
source("CHEMDisplay.R")
source("CHEMBioassay.R")


listCHEMMethods <- function() {
  print("CHEMIO.R")
  print("importSDFFromCID(CIDvector) : returns SDFset")
  print("loadCHEMDF(fileName) : returns data.frame")
  print("loadCIDFile(fileName) : returns data.frame (drugName, CID)")
  print("loadSDFFile(fieName) : returns SDFset")
  print("loadClusterMatrix(fileName) : returns matrix")
  print("saveCHEMDF(df,fileName)")
  print("saveClusterMatrix(simMatrix, fileName)")
  print("saveSDFFile(SDFset, fileName)")
  print(" ")

  print("CHEMAnalysis.R")
  print("clusterWithFMCSAromatic(SDFset) : returns a matrix of similarities")
  print("clusterWithFMCSStatic(SDFset) : returns a matrix of similarities")
  print("clusterWithFMCSRing(SDFset) : returns a matrix of similarities")
  print("")

  print("CHEMDisplay.R")
  print("displayClusterDend(simMatrix)")

}
