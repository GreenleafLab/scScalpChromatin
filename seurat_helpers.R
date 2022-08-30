# Helper functions for Seurat

suppressPackageStartupMessages({
  library(Seurat)
  library(magrittr)
  library(DoubletFinder)
  library(celda)
  library(dplyr)
  library(tidyr)
  library(ggrastr)
})


scRNAdataPreProcessing <- function(
  obj, objPrefix, plotDir, # Seurat object, prefix, and plot dir for qc plots
  minFeatures=200, maxFeatures=Inf, minCounts=1000, maxCounts=Inf, maxPctMito=10, # Basic quality filters
  nfeatures=2500, dims=1:15, res=0.5, # Seurat clustering parameters
  runDoubletFinder=TRUE, # Should DoubletFinder be run?
  estDubRate=0.075, # DoubletFinder parameters
  runDecontX=TRUE, # Should DecontX be run?
  assays="RNA", # Which assays should be kept?
  ncores=1, use_logfile=TRUE
  ){

  # Perform some basic filtering on a seurat object.
  # Seurat object should be created from a single sample (e.g. one 10x run)
  # Optionally run other pre-processing tools:
  #
  # - DoubletFinder to estimate likely doublets 
  #   https://github.com/chris-mcginnis-ucsf/DoubletFinder
  #   https://www-cell-com.stanford.idm.oclc.org/cell-systems/fulltext/S2405-4712(19)30073-0
  #
  # - DecontX to reduce ambient RNA contamination
  #   https://github.com/campbio/celda
  #   https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6
  # 
  ################################
  # obj = seurat object
  # objPrefix = prefix for raq_qc plots
  # plotDir = directory for plotting
  # minFeatures = minimum number of features (genes) per cell
  # maxFeatures = maximum number of features (genes) per cell
  # minCounts = minimum number of UMIs per cell
  # maxCounts = maximum number of UMIs per cell
  # maxPctMito = maximum mitochondrial read percentage
  # nfeatures = number of variable features to be used in DoubletFinder
  # dims = which PCA dimensions will be used in DoubletFinder
  # estDubRate = estimated percent doublets (DF.classify will identify exactly this percent as doublets!)

  if(use_logfile){
    logfile <- paste0(plotDir, sprintf("/%s_preprocess_log_%s.txt", objPrefix, format(Sys.time(), "%Y%m%d-%H%M%S")))
    con <- file(logfile, open = "wt")
    sink(con, type="output")
    sink(con, type="message")
  }

  # Add percent.mt 
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  cellsBeforeFiltering <- dim(obj)[2]
  
  # Save some quick plots of raw qc
  histBreaks <- 100

  pdf(paste0(plotDir, sprintf("/%s_nCountHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), log10nCount=log10(obj$nCount_RNA))
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", lower_lim=log10(minCounts), upper_lim=log10(maxCounts))
  print(p)
  dev.off()

  pdf(paste0(plotDir, sprintf("/%s_nFeatureHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), log10nFeatures=log10(obj$nFeature_RNA))
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", lower_lim=log10(minFeatures), upper_lim=log10(maxFeatures))
  print(p)
  dev.off()

  pdf(paste0(plotDir, sprintf("/%s_pctMitoHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), PctMito=obj$percent.mt)
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", upper_lim=maxPctMito)
  print(p)
  dev.off()

  # Perform basic hard threshold filters
  message("Will filter based on:")
  message(sprintf("%s < unique genes < %s", minFeatures, maxFeatures))
  message(sprintf("%s < UMIs (counts) < %s", minCounts, maxCounts))
  message(sprintf("Percent mitochondrial < %s", maxPctMito))

  obj <- subset(obj, 
    subset = (
      nFeature_RNA > minFeatures & 
      nFeature_RNA < maxFeatures & 
      nCount_RNA > minCounts & 
      nCount_RNA < maxCounts & 
      percent.mt < maxPctMito
    )
  )
  cellsAfterFiltering <- dim(obj)[2]
  message(sprintf("%s filtered down to %s (%s%% remaining)", 
    cellsBeforeFiltering, cellsAfterFiltering, 
    round(100*(cellsAfterFiltering/cellsBeforeFiltering), 2)))

  # Perform standard Seurat pre-processing:
  obj <- seuratPreProcess(obj, selectionMethod="vst", nFeatures=nfeatures, dims=dims)

  # Run DoubletFinder, if indicated
  if(runDoubletFinder){
  	message("Running DoubletFinder...")
    obj <- runDoubletFinder(obj, dims, estDubRate=estDubRate, ncores=ncores)
  }

  # Run DecontX, if indicated
  if(runDecontX){
  	message("Running DecontX...")
  	obj <- runDecontX(obj)
  	assays <- c(assays, "origCounts")
  }

  # Return filtered Seurat object
  obj <- DietSeurat(obj, counts=TRUE, data=TRUE, scale.data=FALSE, assays=assays)

  # Close connections
  message("Finished preprocessing...")
  if(use_logfile){
    on.exit({ sink(type = "message"); sink(type = "output"); close(con) })
  }

  # Return processed Seurat object
  return(obj)
}


seuratPreProcess <- function(obj, selectionMethod="vst", nFeatures=2000, dims=1:10, res=0.5, seed=1){
  # perform standard Seurat preprocessing 
  # (Normalization, Scaling, Dimensionality Reduction, Clustering)
  obj <- NormalizeData(obj)
  obj <- ScaleData(obj, verbose=FALSE)
  obj <- FindVariableFeatures(obj, selection.method=selectionMethod, nfeatures=nFeatures)
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims=dims)
  obj <- FindNeighbors(obj, dims=dims)
  obj <- FindClusters(obj, resolution=res, random.seed=1)
  return(obj)
}


runDoubletFinder <- function(obj, dims, estDubRate=0.075, ncores=1){
  # Run DoubletFinder on a provided (preprocessed) Seurat object
  # Return the seurat object with the selected pANN parameter and the 
  # DoubletFinder doublet classifications

  ### pK Identification (parameter-sweep) ###
  # "pK ~ This defines the PC neighborhood size used to compute pANN (proportion of artificial nearest neighbors), 
  # expressed as a proportion of the merged real-artificial data. 
  # No default is set, as pK should be adjusted for each scRNA-seq dataset"

  sweep.res.list <- paramSweep_v3(obj, PCs=dims, sct=FALSE, num.cores=ncores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))

  # Get expected doublets (DF.classify will identify exactly this percent as doublets!)
  nExp_poi <- round(estDubRate * length(Cells(obj)))

  # DoubletFinder:
  obj <- doubletFinder_v3(obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  # Rename results into more useful annotations
  pann <- grep(pattern="^pANN", x=names(obj@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(obj@meta.data), value=TRUE)
  obj$pANN <- obj[[pann]]
  obj$DF.classify <- obj[[classify]]
  obj[[pann]] <- NULL
  obj[[classify]] <- NULL

  return(obj)
}


runDecontX <- function(obj, seed=1){
  # Run DecontX on a provided Seurat object
  # From the DecontX vignette: 
  # "**Only the expression profile of *"real"* cells after cell calling are required to run DecontX. 
  # Empty cell droplet information (low expression cell barcodes before cell calling) are not needed.**"

  # DecontX can take either `SingleCellExperiment` object... or a single counts matrix as input. 
  # `decontX` will attempt to convert any input matrix to class `dgCMatrix` before beginning any analyses.
  counts <- GetAssayData(object = obj, slot = "counts")
  clusters <- Idents(obj) %>% as.numeric()

  # Run on only expressed genes
  x <- counts[rowSums(counts)>0,]
  message(sprintf("Running decontX on %s cells with %s non-zero genes...", dim(x)[2], dim(x)[1]))
  decon <- decontX(x, z=clusters, verbose=TRUE, seed=seed)

  # Save desired information back to Seurat Object
  # We will place the estimated 'decontaminated counts' in place of the original counts ('RNA')
  # and keep the original counts as a separate assay called 'origCounts'
  obj[["origCounts"]] <- CreateAssayObject(counts = counts)
  newCounts <- decon$decontXcounts
  # Add back unexpressed genes and sort according to original counts
  newCounts <- rbind(newCounts, counts[rowSums(counts)==0,])[rownames(counts),]
  obj[["RNA"]]@counts <- as(round(newCounts), "sparseMatrix")
  obj$estConp <- decon$contamination # Estimated 'contamination proportion, 0 to 1'

  return(obj)
}



plotClusterQC <- function(obj, subgroup, plotDir, pointSize=1.0, barwidth=0.9, sampleCmap=NULL, diseaseCmap=NULL){
  
  # Plot basic clustering plots
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

  ### Bar plot cluster counts ###
  tabDF <- base::table(obj$Clusters) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")

  pdf(paste0(plotDir, sprintf("/clusterBarPlot_%s.pdf", subgroup)))
  print(qcBarPlot(tabDF, cmap=qualcmap, barwidth=barwidth))
  dev.off()

  clustBySamp <- fractionXbyY(obj$Clusters, obj$Sample, add_total=TRUE, xname="Cluster", yname="Sample")

  pdf(paste0(plotDir, sprintf("/clustBySampleBarPlot_%s.pdf", subgroup)))
  print(stackedBarPlot(clustBySamp, cmap=sampleCmap, namedColors=namedSampCmap, barwidth=barwidth))
  dev.off()

  ### Stacked bar plot fraction disease in clusters ###
  clustByDisease <- fractionXbyY(obj$Clusters, obj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="Disease")

  pdf(paste0(plotDir, sprintf("/clustByDiseaseBarPlot_%s.pdf", subgroup)))
  print(stackedBarPlot(clustByDisease, cmap=diseaseCmap, namedColors=namedDiseaseCmap, barwidth=barwidth))
  dev.off()

  ### Cluster UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$Clusters)
  # Randomize cells before plotting UMAP
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/cluster_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=qualcmap, point_size=pointSize))
  dev.off()

  ### Sample UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$Sample)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/sample_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=sampleCmap, namedColors=namedSampCmap, point_size=pointSize))
  dev.off()

  ### Disease UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$diseaseStatus)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/disease_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=diseaseCmap, namedColors=namedDiseaseCmap, point_size=pointSize))
  dev.off()

  ### Percent Mitochondrial UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$percent.mt)
  pdf(paste0(plotDir, sprintf("/pctMito_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="quantitative", cmap=quantcmap, point_size=pointSize))
  dev.off()

  ### nCounts (UMIs) UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), log10(obj$nCount_RNA))
  pdf(paste0(plotDir, sprintf("/log10nCount_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="quantitative", cmap=quantcmap, point_size=pointSize))
  dev.off()

  ### nFeatures (unique genes) UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$nFeature_RNA)
  pdf(paste0(plotDir, sprintf("/nFeature_RNA_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="quantitative", cmap=quantcmap, point_size=pointSize))
  dev.off()

  ### Cell Cycle Phase UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$Phase)
  pdf(paste0(plotDir, sprintf("/ccPhase_RNA_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=qualcmap, point_size=pointSize))
  dev.off()
}



