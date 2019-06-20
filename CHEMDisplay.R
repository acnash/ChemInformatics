displayClusterDend <- function(simMatrix, Namevector) {
  colnames(simMatrix) <- Namevector
  row.names(simMatrix) <- Namevector
  hCluster <- hclust(as.dist(1-simMatrix), method="complete") 
  dend <- as.dendrogram(hCluster)
  dend %>% dendextend::set('labels_cex',0.5) %>% highlight_branches_col(viridis(100)) 
  return(dend)
}

#========================================================================

displayClusterHeatMap <- function(simMatrix, compoundName) {
  colnames(simMatrix) <- compoundName
  row.names(simMatrix) <- compoundName
  hc <- hclust(as.dist(1-simMatrix), method="single")
  heatmap <- heatmap.2(1-simMatrix, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "darkred", "yellow", "white"), density.info="none", trace="none")
  return(heatmap)
}

#========================================================================

displayFMCSSize <- function(df, fileNameJPEG) {
  #jpeg(fileNameJPEG, units="in", width=24, height=8, res=300)
  g <- ggplot(df, aes(df$drug, df$size)) + geom_point(shape=21, colour="#999999", fill="#56B4E9", size=4, stroke=1) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text = element_text(size=8)) + ylab("MCS Size (number of atoms)") + xlab("Drug")
  #g
  #dev.off()
  return(g)
}

#========================================================================

getNCBIFractionActivity <- function(drugTargetsMat, title=NULL) {
  df <- as.data.frame(drugTargetsMat)
  g <- ggplot(data=df, aes(y=df$fraction_active, x=row.names(df))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text = element_text(size=8)) + ylab("Total Fraction of Activity") + xlab("NCBI Target Identifier") + geom_bar(color="#999999", fill="#56B4E9", stat="identity")
  if(!is.null(title)) {
    g <- g + ggtitle(title)
  }
  return(g)
}
