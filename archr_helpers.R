# Helper functions for ArchR

suppressPackageStartupMessages({
  library(ArchR)
  library(magrittr)
})


getMatrixValuesFromProj <- function(proj, matrixName="GeneScoreMatrix", names=NULL, imputeMatrix=FALSE){
  # Return the imputed matrix from an ArchR project
  # Must have already added imputedWeights, etc.
  # Names is a vector of feature names to return. If not provided, will use all available features
  # Warning though: imputing the matrix for all features may take a very long time

  # Returns a summarized experiment:
  se <- getMatrixFromProject(proj, useMatrix = matrixName, binarize = FALSE)

  # Get mat with cell names and row names
  mat <- assays(se)[[matrixName]]
  colnames(mat) <- rownames(colData(se))
  rownames(mat) <- rowData(se)$name # All matrix rowData has a name field

  # Subset by provided names
  if(!is.null(names)){
    # Check if any names are invalid
    validNames <- names[names %in% rownames(mat)]
    if(any(!names %in% validNames)){
      invalidNames <- names[!names %in% validNames]
      message(sprintf("Warning! name(s) %s are not present in matrix!", paste(invalidNames, collapse=',')))
    }
    mat <- mat[validNames,]
  }

  # Impute matrix values
  if(imputeMatrix){
    message("Imputing matrix...")
    imputeWeights <- getImputeWeights(proj)
    mat <- ArchR::imputeMatrix(mat = as.matrix(mat), imputeWeights = imputeWeights)
  }
  mat
}


getClusterPeaks <- function(proj, clusterNames, peakGR=NULL, replicateScoreQuantileCutoff=0, originalScore=FALSE){
  # This function will return the subset of peaks from the full ArchR project that were 
  # initially called on the clusters provided in clusterNames.
  ######################################################################################
  # proj = ArchR project
  # clusterNames = name or names of clusters to pull peaks from. These cluster names must
  #   match the cluster names originally used to call peaks
  # peakGR = the ArchR peak genomic range obtained using 'getPeakSet'. Will return a subset of
  #   these peaks that overlap the peaks originally called using the provided clusters
  # replicateScoreQuantileCutoff = A numeric quantile cutoff for selecting peaks. 
  if(is.null(peakGR)){
    peakGR <- getPeakSet(proj)
  }
  peakDir <- paste0(proj@projectMetadata$outputDirectory, "/PeakCalls")
  calledPeaks <- lapply(clusterNames, function(x){
      readRDS(paste0(peakDir, sprintf("/%s-reproduciblePeaks.gr.rds", x)))
    }) %>% as(., "GRangesList") %>% unlist()
  calledPeaks <- calledPeaks[calledPeaks$replicateScoreQuantile >= replicateScoreQuantileCutoff]
  peakGR <- peakGR[overlapsAny(peakGR, calledPeaks)]
  if(originalScore){
    message(sprintf("Getting original scores from clusters..."))
    # Replace 'score' column with the score of the original peak call in this cluster
    # (if multiple clusters, replaces with the maximum score)
    ol <- findOverlaps(peakGR, calledPeaks, type="any", maxgap=0, ignore.strand=TRUE)
    odf <- as.data.frame(ol)
    odf$og_score <- calledPeaks$score[odf$subjectHits]
    score_df <- odf %>% group_by(queryHits) %>% summarize(max_score=max(og_score)) %>% as.data.frame()
    peakGR$score[score_df$queryHits] <- score_df$max_score
  }
  peakGR
}


buildUMAPdfFromArchR <- function(proj, cellColData=NULL, embeddingName="UMAP", 
  useCells=NULL, dataMat=NULL, featureName=NULL, shuffle=TRUE, 
  lowerPctLim=NULL, upperPctLim=NULL){
  # Return a three column UMAP df from an ArchR project
  # If cellColData is not null, return the indicated column
  # dataMat is a pre-populated cell x feature matrix of values to plot. 
  # The featureName indicates which one

  # Get UMAP coordinates first:
  df <- proj@embeddings[[embeddingName]]$df
  if(is.null(useCells)){
    useCells <- rownames(df)
  }
  colnames(df) <- c("UMAP1", "UMAP2")
  df <- df[useCells,] %>% as.data.frame()
  if(!is.null(cellColData)){
    df[,3] <- proj@cellColData[useCells,cellColData] %>% as.vector()
    colnames(df) <- c("UMAP1", "UMAP2", cellColData)
  }
  if(!is.null(dataMat) & !is.null(featureName)){
    df <- merge(df, dataMat, by=0, all=TRUE)
    df <- df[,c("UMAP1", "UMAP2", featureName)] 
  }
  if(shuffle){
    df <- df[sample(nrow(df), replace=FALSE),]
  }
  # Force limits if indicated
  if(!is.null(lowerPctLim)){
    lowerLim <- quantile(df[,3], probs=c(lowerPctLim))
    df[,3][df[,3] <= lowerLim] <- lowerLim
  }
  if(!is.null(upperPctLim)){
    upperLim <- quantile(df[,3], probs=c(upperPctLim))
    df[,3][df[,3] >= upperLim] <- upperLim
  }
  df
}


scoreGeneSet <- function(expr, geneSet){
  # Generate scores for each cell in expr matrix (log2TP10K, genes x cells)
  # See: Smillie et al. Cell 2019

  # Subset expr matrix by genes in geneSet:
  validGenes <- geneSet[geneSet %in% rownames(expr)]
  subExpr <- expr[validGenes,]

  # Remove any genes that have no expression in any cells
  subExpr <- subExpr[rowSums(subExpr) > 0,]

  # Prevent highly expressed genes from dominating gene score signature by
  # scaling each gene by its root mean squared expression
  scaledSubExpr <- subExpr %>% t() %>% scale(., center=FALSE) %>% t()

  # Signature score is the mean scaled expression across all genes in signature
  scores <- colMeans(scaledSubExpr)
  return(scores)
}


# Functions for creating 'low-overlapping aggregates' of cells

computeKNN <- function(data=NULL, query=NULL, k=50, includeSelf=FALSE, ...){
  # Compute KNN for query points (usually a reduced dims matrix)
  # This returns a matrix of indices mapping query to neighbors in data
  # If query has n cells (rows) and k = 50, will be a n x 50 matrix
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}


getLowOverlapAggregates <- function(proj, target.agg=500, k=100, overlapCutoff=0.8, dimReduc="IterativeLSI", seed=1){
  # Generate low-overlapping aggregates of cells
  ##############################################
  # proj = ArchR project
  # target.agg = number of target aggregates (before filtering)
  # k = number of cells per aggreagate
  # overlapCutoff = Maximum allowable overlap between aggregates
  set.seed(seed)

  # Get reduced dims:
  rD <- getReducedDims(proj, reducedDims=dimReduc)

  # Subsample
  idx <- sample(seq_len(nrow(rD)), target.agg, replace = !nrow(rD) >= target.agg)

  # Get KNN Matrix:
  knnObj <- computeKNN(data=rD, query=rD[idx,], k=k)

  # Check whether aggregates pass overlap cutoff
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]

  # Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  # Name aggregates and return as a df of cell ids x aggs
  names(knnObj) <- paste0("agg", seq_len(length(knnObj)))
  knnDF <- data.frame(knnObj)[,c(3,2)]
  colnames(knnDF) <- c("cell_name", "group")
  knnDF$cell_name <- as.character(knnDF$cell_name)
  knnDF
}


# Cluster visualization helpers

relabelClusters <- function(proj, clusterName="Clusters"){
  # Relabel clusters to be ordered by cluster size

  ogClusters <- getCellColData(proj)[[clusterName]]
  tabDF <- base::table(ogClusters) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")
  tabDF["NewClusters"] <- rank(-tabDF$count)
  swapVec <- paste0("C", tabDF$NewClusters)
  names(swapVec) <- tabDF$Clusters

  # Now replace cluster names
  newClust <- sapply(ogClusters, function(x) swapVec[x]) %>% unname()
  proj <- addCellColData(proj, data=newClust, name=clusterName, cells=getCellNames(proj), force=TRUE)
  return(proj)
}


visualizeClustering <- function(proj, pointSize=0.75, prefix="", clusterName="Clusters", sampleName="Sample2", embedding="UMAP", 
  sampleCmap=NULL, diseaseCmap=NULL, barwidth=0.9){
  # Plot various clustering results
  # Set colormap
  qualcmap <- cmaps_BOR$stallion
  quantcmap <- cmaps_BOR$solarExtra
  namedSampCmap <- TRUE
  namedDiseaseCmap <- TRUE

  if(is.null(sampleCmap)){
    sampleCmap <- qualcmap
    namedSampCmap <- FALSE
  }

  if(is.null(diseaseCmap)){
    diseaseCmap <- qualcmap
    namedDiseaseCmap <- FALSE
  }

  # Plot the UMAPs by Sample and Cluster:
  p1 <- plotEmbedding(proj, colorBy="cellColData", name=sampleName, embedding=embedding, plotAs="points", size=pointSize, pal=sampleCmap, labelMeans=FALSE)
  p2 <- plotEmbedding(proj, colorBy="cellColData", name=clusterName, embedding= embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  p3 <- plotEmbedding(proj, colorBy="cellColData", name="diseaseStatus", embedding=embedding, plotAs="points", size=pointSize, pal=diseaseCmap, labelMeans=FALSE)
  proj@cellColData$log10nFrags <- log10(proj@cellColData$nFrags)
  p4 <- plotEmbedding(proj, colorBy="cellColData", name="log10nFrags", embedding=embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  p5 <- plotEmbedding(proj, colorBy="cellColData", name="TSSEnrichment", embedding=embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  p6 <- plotEmbedding(proj, colorBy="cellColData", name="DoubletScore", 
    embedding = embedding, plotAs="points", size=pointSize, labelMeans=FALSE, imputeWeights=getImputeWeights(proj))
  p7 <- plotEmbedding(proj, colorBy = "cellColData", name="cellCallUncertainty", 
    embedding = embedding, plotAs="points", size=pointSize, labelMeans=FALSE, imputeWeights=getImputeWeights(proj))
  ggAlignPlots(p1,p2,p3,p4,p5,p6,p7, type="h")
  plotPDF(p1,p2,p3,p4,p5,p6,p7, name = paste0(prefix,"Plot-UMAP-Sample-Clusters.pdf"), ArchRProj=proj, addDOC=FALSE, width=5, height=5)

  # Non-ArchR plots:
  plotDir <- paste0(proj@projectMetadata$outputDirectory, "/Plots")

  # Bar plot cluster counts
  clustVec <- getCellColData(proj)[[clusterName]] %>% gsub("[^[:digit:].]", "", .) %>% as.numeric()
  tabDF <- base::table(clustVec) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")

  pdf(paste0(plotDir,sprintf("/%sclusterBarPlot.pdf", prefix)))
  print(qcBarPlot(tabDF, cmap=qualcmap, barwidth=barwidth))
  dev.off()

  # Stacked bar plot fraction samples in clusters
  clustBySamp <- fractionXbyY(clustVec, proj$Sample2, add_total=TRUE, xname="Cluster", yname="Sample")

  pdf(paste0(plotDir, sprintf("/%sclustBySampleBarPlot.pdf", prefix)))
  print(stackedBarPlot(clustBySamp, cmap=sampleCmap, namedColors=namedSampCmap, barwidth=barwidth))
  dev.off()

  # Stacked bar plot fraction disease in clusters
  diseaseBySamp <- fractionXbyY(clustVec, proj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="diseaseStatus")

  pdf(paste0(plotDir, sprintf("/%sclustByDiseaseBarPlot.pdf", prefix)))
  print(stackedBarPlot(diseaseBySamp, cmap=diseaseCmap, namedColors=namedDiseaseCmap, barwidth=barwidth))
  dev.off()

  return(proj)
}


# Functions for working with peak to gene linkages

getP2G_GR <- function(proj, corrCutoff=NULL, varCutoffATAC=0.25, varCutoffRNA=0.25, filtNA=TRUE){
  # Function to get peaks and genes involved in peak to gene links
  # (See: https://github.com/GreenleafLab/ArchR/issues/368)
  ############################################################
  # proj: ArchR project that alreayd has Peak2GeneLinks
  # corrCutoff: minimum numeric peak-to-gene correlation to return
  # varCutoffATAC: minimum variance quantile of the ATAC peak accessibility when selecting links
  # varCutoffRNA: minimum variance quantile of the RNA gene expression when selecting links
  p2gDF <- metadata(proj@peakSet)$Peak2GeneLinks
  p2gDF$symbol <- mcols(metadata(p2gDF)$geneSet)$name[p2gDF$idxRNA] %>% as.character()
  p2gDF$peakName <- (metadata(p2gDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2gDF$idxATAC]
  # Remove peaks with 'NA' correlation values
  if(filtNA){
    p2gDF <- p2gDF[!is.na(p2gDF$Correlation),]
  }
  if(!is.null(corrCutoff)){
    p2gDF <- p2gDF[(p2gDF$Correlation > corrCutoff),]
  }
  # Filter by variance quantile
  p2gDF <- p2gDF[which(p2gDF$VarQATAC > varCutoffATAC & p2gDF$VarQRNA > varCutoffRNA),]
  # The genomic range contains just the peak ranges:
  p2gGR <- metadata(p2gDF)$peakSet[p2gDF$idxATAC]
  mcols(p2gGR) <- p2gDF
  p2gGR
}


grLims <- function(gr){
  # Get the minimum and maximum range from a GR
  if(length(gr) == 0){
    return(NA)
  }
  starts <- start(gr)
  ends <- end(gr)
  c(min(starts, ends), max(starts, ends))
}


getP2Gregions <- function(proj, genes, p2gGR=NULL, corrCutoff=0.4, buffer_space=0.05, min_width=25000, ...) {
  # Function to get regions containing entire peak to gene region,
  # i.e. a GR that contains all peak to gene links
  ###############################################################
  # p2gGR: genomic range containing all peak to gene links
  # genes: vector of genes to look up
  # buffer_space: fraction of total length to expand on each side of region

  # Get gene GR from ArchR project
  geneGR <- promoters(getGenes(proj)) # Promoters gets 2kb upstream and 200bp downstream
  geneGR <- geneGR[!is.na(geneGR$symbol)]

  # if p2gGR not provided, pull it from ArchR project
  if(is.null(p2gGR)){
    p2gGR <- getP2G_GR(proj, corrCutoff=corrCutoff, ...)
  }

  # Now for each gene, construct GR of all loops and gene TSS
  resultGR <- geneGR[match(genes, geneGR$symbol)]
  start(resultGR) <- sapply(resultGR$symbol, function(x){
      min(grLims(resultGR[resultGR$symbol == x]), grLims(p2gGR[p2gGR$symbol == x]), na.rm=TRUE)
    })
  end(resultGR) <- sapply(resultGR$symbol, function(x){
      max(grLims(resultGR[resultGR$symbol == x]), grLims(p2gGR[p2gGR$symbol == x]), na.rm=TRUE)
    })

  # Finally, resize by buffer space
  resultGR <- resize(resultGR, width=width(resultGR) + buffer_space*width(resultGR), fix="center")
  resultGR <- resize(resultGR, width=ifelse(width(resultGR) > min_width, width(resultGR), min_width), fix="center")
  resultGR
}


