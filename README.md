# ChemInformatics
Cheminformatics R wrapper code and work flows. There is no actual 'new' science here, just a list of functions for the direct purpose of comparing and clustering drugs by their substructure and CMAP transcriptional response. 

This is a 5% project, i.e., just something I am interested in and I'll do what I can to build on this work but it isn't a priority.

Please cite this github repo if this code is used in your work. 

Files you will need
-----------
pubchem_protein_only.sqlite</br>
available_genes.txt</br>
drugs.txt</br>
CD_signatures_LM_42809x978.gctx</br>
Drugs_metadata.csv</br>
CD_signature_metadata.csv</br>

IMPORTANT
-----------
The L1000 gene transcriptional response to drug code is in drug_clustering.R. However, this file is very much tied into my machine. Until I get around to building this into an API, you'd need to hack this file apart.

List of wrapper functions and load into the R environment
-----------
You are going to need a Net connection.</br>
Load the R files into the R environment using:</br>
source("CHEMAPI.R")

For a list of the most recent available <b>(this list is old and outdated)</b> functions execute:</br>
listChemMethods()

<b>CHEMIO.R</b></br>
importSDFFromCID(CIDvector) : returns SDFset</br>
loadCHEMDF(fileName) : returns data.frame</br>
loadCIDFile(fileName) : returns data.frame (drugName, CID)</br>
loadSDFFile(fieName) : returns SDFset</br>
loadClusterMatrix(fileName) : returns matrix</br>
saveCHEMDF(df,fileName)</br>
saveClusterMatrix(simMatrix, fileName)</br>
saveSDFFile(SDFset, fileName)</br>

<b>CHEMAnalysis.R</b></br>
calculateFMCS(sdfObject1, sdfObject2, au=2, bu=1) : returns an fmcs object</br>
calculateNumAtomsToMCSReference(referenceSDF, ignoreIndex, SDFset, drugNameVector, au=2, bu=1) : returns data.frame for plotting</br>
clusterWithFMCSAromatic(SDFset, au=2, bu=1, overlapCoefficient=TRUE) : returns a matrix of similarities</br>
clusterWithFMCSStatic(SDFset, au=2, bu=1, overlapCoefficient=TRUE) : returns a matrix of similarities</br>
clusterWithFMCSRing(SDFset, au=2, bu=1, overlapCoefficient=TRUE) : returns a matrix of similarities</br>

<b>CHEMDisplay.R</b></br>
displayClusterDend(simMatrix) : returns dendrogram objecct for plot()</br>
displayClusterHeatMap(simMatrix, CIDvector, title=NULL) : returns heatmap object for plot()</br>
displayFMCSSize(df, fileNameJPEG) : returns ggplot object</br>
getNCBIFractionActivity(drugTargetsMat, title=NULL) : returns ggplot object</br>

<b>CHEMBioassay.R</b></br>
getBioassayDatabase(DBLocation) : returns a db object, don't forget to close it</br>
getEnsemblProteinDetails(uniProtIDs, attributeCharVector, filtersCharVector, drugTargets) : returns list of protein details</br>
getProteinDrugTargets(db, CID) : returns a dataframe (within a list) showing the drug targets (protein), the number of total assay screens and the total fraction of activity</br>
getUniProtIDs(NCBIMatrix, db) : returns a list of uniprot IDs</br>
clusterCompoundsByActivityProfile(db, compoundCIDs) : returns hcluster object</br>


Tools
-----------
-ChemmineR

-BioassayR

-BiomaRt

-cmapR


DB Access
-----------
-PubChem

-Ensembl

-UniProt

-CLUE

-Drugbank.ca
