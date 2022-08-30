#!/usr/bin/env Rscript

#############################
# Processing scalp data
#############################

#Load ArchR (and associated libraries)
library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 0.25
barwidth <- 0.9

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_proj <- readRDS("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scRNA_preprocessing/preprocessing_output/scalp.rds")
plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

# Color Maps
broadClustCmap <- readRDS(paste0(scriptPath, "/scalpClusterColors.rds")) %>% unlist()
atacNamedClustCmap <- readRDS(paste0(scriptPath, "/scATAC_NamedClust_cmap.rds")) %>% unlist()
rnaNamedClustCmap <- readRDS(paste0(scriptPath, "/scRNA_NamedClust_cmap.rds")) %>% unlist()
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
rna_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(rna_proj$Sample)] %>% unlist()
atac_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(atac_proj$Sample2)] %>% unlist()

# Get label cmaps
source(paste0(scriptPath, "/cluster_labels.R"))
atacLabelClustCmap <- atacNamedClustCmap
names(atacLabelClustCmap) <- unlist(atac.NamedClust)[names(atacNamedClustCmap)]
rnaLabelClustCmap <- rnaNamedClustCmap
names(rnaLabelClustCmap) <- unlist(rna.NamedClust)[names(rnaNamedClustCmap)]

disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

##########################################################################################
# Integrate RNA and ATAC data
##########################################################################################

## WARNING! Very memory intensive

# Ensure Impute Weights are up to date
atac_proj <- addImputeWeights(atac_proj)

# Unconstrained integration with RNA:
set.seed(1)
atac_proj <- addGeneIntegrationMatrix(
    ArchRProj = atac_proj, 
    force = TRUE,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    sampleCellsATAC = 25000, # Default for both was 10000
    sampleCellsRNA = 25000,
    nGenes = 3000, # Default was 2000
    seRNA = rna_proj, # Can be a seurat object
    addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
    groupRNA = "NamedClust", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this function.
    nameCell = "RNA_paired_cell", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "NamedClust_RNA", #Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore" #Name of column where prediction score from scRNA
)

cM <- as.matrix(confusionMatrix(atac_proj$NamedClust, atac_proj$NamedClust_RNA))

# To see the most likely ATAC label:
lcM <- cM
colnames(lcM) <- unlist(rna.NamedClust)[colnames(lcM)]
apply(lcM,1, function(x)colnames(lcM)[which.max(x)])

atacOrder <- c(
    "aTc2", # "CD8.Tc"
    "aTc1", # "CD4.Tc"
    "aTc3", # "Tregs"
    "aMy1", # "DCs_1"
    "aMy2", # "Macs_1"
    "aMy3", # "CLEC9a.DC"
    "aKc1", # "Basal.Kc_1"
    "aKc2", # "Spinous.Kc_1"
    "aKc3", # "Spinous.Kc_2"
    "aKc4", # "HF.Kc_1"
    "aKc5", # "HF.Kc_2"
    "aKc6", # "HF.Kc_3",
    "aKc7", # "HF.Kc_4",
    "aFb1", # "D.Fib" # Papillary/Reticular dermal fibroblasts
    "aFb2", # "D.Sheath" # Dermal sheath
    "aMu1", # "Muscle_1"
    "aMu2", # "Muscle_2" 
    "aVe1", # "Vas.Endo_1"
    "aVe2", # "Vas.Endo_2"
    "aLe1", # "Lymph.Endo"
    "aMe1", # "Melanocytes"
    "aBc1" # "B.cells"
)

rnaOrder <- c(
    "rTc3", # "CD8.Tc",
    "rTc2", # "CD4.Tc",
    "rTc1", # "Tregs",
    "rMy1", # "DCs_1",
    "rMy2", # "Macs_1",
    "rMy3", # "CLEC9a.DC",
    "rMy4", # "M1_Macs",
    #"rMa1", # "Mast", # Mast cells
    "rKc5", # "Basal.Kc_2",
    "rKc4", # "Basal.Kc_1",
    "rKc1", # "Spinous.Kc_1",
    "rKc2", # "Spinous.Kc_2",
    "rKc3", # "HF.Kc_1",
    "rFb1", # "D.Fib", # Papillary/Reticular dermal fibroblasts
    "rFb2", # "D.Sheath", # Dermal sheath
    "rMu1", # "Muscle_1",
    "rMu2", # "Muscle_2", 
    "rVe1", # "Vas.Endo",
    "rLe1", # "Lymph.Endo",
    "rMe1", # "Melanocytes",
    "rBc1" # "Bc",
    #"rMe2" # "McSC"
)

atacOrder <- atacOrder[atacOrder %in% rownames(cM)]
rnaOrder <- rnaOrder[rnaOrder %in% colnames(cM)]

cM <- cM[atacOrder, rnaOrder]

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}
cM <- new_cM

# cM <- prettyOrderMat(t(cM),clusterCols=FALSE)$mat %>% t()

rownames(cM) <- unlist(atac.NamedClust)[rownames(cM)]
colnames(cM) <- unlist(rna.NamedClust)[colnames(cM)]

pdf(paste0(plotDir, "/ATAC-RNA-integration-cM-heatmap.pdf"), width=8, height=8)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <- HeatmapAnnotation(atac_cluster=rownames(cM),col=list(atac_cluster=atacLabelClustCmap), 
  which="row", show_legend=c("atac_cluster"=FALSE))
ta <- HeatmapAnnotation(rna_cluster=colnames(cM),col=list(rna_cluster=rnaLabelClustCmap), 
  show_legend=c("rna_cluster"=FALSE))
hm <- BORHeatmap(
  cM, 
  limits=c(0,1), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$whitePurple,
  left_annotation = ra,
  top_annotation = ta,
  row_names_side = "left",
  width = ncol(cM)*unit(0.5, "cm"),
  height = nrow(cM)*unit(0.5, "cm"),
  border_gp=gpar(col="black"), # Add a black border to entire heatmap
  legendTitle="Jaccard Index"
  )
draw(hm)
dev.off()

# Add label named clust
atac_proj$NLabelClust <- unlist(atac.NamedClust)[atac_proj$NamedClust]
atac_proj$NLabelClust_RNA <- unlist(rna.NamedClust)[atac_proj$NamedClust_RNA]

# Plot predicted score by original ATAC and RNA cluster and by Sample
pScoreLim <- c(0.0,1.0)
plotList <- list()
plotList[[1]] <- plotGroups(ArchRProj = atac_proj, 
  groupBy = "Sample2", 
  colorBy = "colData", 
  name = "predictedScore",
  ylim = pScoreLim,
  pal = atac_sample_cmap
)
plotList[[2]] <- plotGroups(ArchRProj = atac_proj, 
  groupBy = "NLabelClust", 
  colorBy = "colData", 
  name = "predictedScore",
  ylim = pScoreLim,
  pal = atacLabelClustCmap
)

plotList[[3]] <- plotGroups(ArchRProj = atac_proj, 
  groupBy = "NLabelClust_RNA", 
  colorBy = "colData", 
  name = "predictedScore",
  ylim = pScoreLim,
  pal = rnaLabelClustCmap
)
plotPDF(plotList = plotList, name = "prediction_accuracy_by_group", width = 4, height = 4,  ArchRProj = atac_proj, addDOC = FALSE)

# Plot the UMAPs by Sample and Cluster:
p1 <- plotEmbedding(atac_proj, colorBy="cellColData", name="Sample2", embedding="UMAP", labelMeans=FALSE,
  size = pointSize, plotAs = "points", pal=atac_sample_cmap)
p2 <- plotEmbedding(atac_proj, colorBy="cellColData", name="diseaseStatus", embedding="UMAP", labelMeans=FALSE,
  size = pointSize, plotAs = "points", pal=disease_cmap)
p3 <- plotEmbedding(atac_proj, colorBy="cellColData", name="NLabelClust", embedding="UMAP", labelMeans = FALSE,
  size = pointSize, plotAs = "points", pal=atacLabelClustCmap)
p4 <- plotEmbedding(atac_proj, colorBy="cellColData", name="NLabelClust_RNA", embedding="UMAP", labelMeans=FALSE,
  size = pointSize, plotAs = "points", pal=rnaLabelClustCmap)
p5 <- plotEmbedding(atac_proj, colorBy = "cellColData", name = "predictedScore", embedding = "UMAP", labelMeans = FALSE,
  size = pointSize, plotAs = "points", colorLimits=c(0.0,1.0))

ggAlignPlots(p1, p2, p3, p4, p5, type = "h")
plotPDF(p1, p2, p3, p4, p5,
    name = "Plot-UMAP-ATACvRNAclusters.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)


# Plot scATAC with ATAC cluster labels
umapPlots <- list()

# ATAC NamedClust on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="NLabelClust")
umapPlots[["namedATAC_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=atacLabelClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="ATAC_labels", useRaster=TRUE)

# ATAC BroadClust on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="BroadClust")
umapPlots[["broadATAC_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=broadClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="ATAC_labels", useRaster=TRUE)

# RNA NamedClust on RNA clustering:
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), unlist(rna.NamedClust)[rna_proj$NamedClust])
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
umapPlots[["namedRNA_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=rnaLabelClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="RNA_labels", useRaster=TRUE)

# RNA BroadClust on RNA clustering:
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), rna_proj$BroadClust)
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
umapPlots[["broadRNA_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=broadClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="RNA_labels", useRaster=TRUE)

# RNA NamedClust on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="NLabelClust_RNA")
umapPlots[["RNA_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=rnaLabelClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="RNA_labels", useRaster=TRUE)

# ATAC Samples on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="Sample2")
umapPlots[["Samples_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=atac_sample_cmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="Samples", useRaster=TRUE)

# ATAC diseaseStatus on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="diseaseStatus")
umapPlots[["diseaseStatus_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=disease_cmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="diseaseStatus", useRaster=TRUE)

# RNA Samples on RNA clustering:
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), rna_proj$Sample)
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
umapPlots[["Samples_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=rna_sample_cmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="Samples", useRaster=TRUE)

# RNA diseaseStatus on RNA clustering:
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), rna_proj$diseaseStatus)
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
umapPlots[["diseaseStatus_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=disease_cmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="diseaseStatus", useRaster=TRUE)

# Predicted Score UMAP
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="predictedScore", lowerPctLim=0.01)
umapPlots[["predScore"]] <- plotUMAP(umapDF, dataType = "quantitative", cmap = cmaps_BOR$solarExtra, 
  point_size=pointSize, covarLabel="predictedScore", colorLims=c(0,1))

pdf(paste0(plotDir,"/clusters_UMAP_ATACvRNA.pdf"), width=10, height=8)
umapPlots
dev.off()

# Save version of each umap with only only major cluster colored
atac_proj$pseudoBC <- atac_proj$BroadClust
rna_proj$pseudoBC <- rna_proj$BroadClust
pList <- list()
for(clust in unique(rna_proj$pseudoBC)){
  # ATAC NamedClust on ATAC clustering:
  umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="NLabelClust", shuffle=FALSE)
  umapDF[atac_proj$pseudoBC != clust, 3] <- NA
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]
  pList[[sprintf("namedATAC_on_ATAC_%s", clust)]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=atacLabelClustCmap, 
    namedColors=TRUE, point_size=pointSize, covarLabel="ATAC_labels", na.value="grey80", useRaster=TRUE)

  # RNA NamedClust on RNA clustering:
  umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), unlist(rna.NamedClust)[rna_proj$NamedClust])
  umapDF[rna_proj$pseudoBC != clust, 3] <- NA
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]
  pList[[sprintf("namedRNA_on_RNA_%s", clust)]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=rnaLabelClustCmap, 
    namedColors=TRUE, point_size=pointSize, covarLabel="RNA_labels", na.value="grey80", useRaster=TRUE)

  # RNA diseaseStatus on RNA clustering:
  umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), rna_proj$diseaseStatus)
  umapDF[rna_proj$pseudoBC != clust, 3] <- NA
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]
  pList[[sprintf("diseaseStatus_on_RNA_%s", clust)]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=disease_cmap, 
    namedColors=TRUE, point_size=pointSize, covarLabel="RNA_labels", na.value="grey80", useRaster=TRUE)
}
pdf(paste0(plotDir,"/diseaseStatus_loo_clusters_UMAP_ATACvRNA.pdf"), width=10, height=8)
pList
dev.off()

# NamedClust by Sample barplot
clustBySamp <- fractionXbyY(atac_proj$NLabelClust, atac_proj$Sample2, add_total=TRUE, xname="Cluster", yname="Sample")
pdf(paste0(plotDir, "/sampleByManualLabelsBarPlot.pdf"))
stackedBarPlot(clustBySamp, cmap=atac_sample_cmap, namedColors=TRUE, barwidth=barwidth)
dev.off()

# BroadClust by Sample barplot
clustBySamp <- fractionXbyY(atac_proj$BroadClust, atac_proj$Sample2, add_total=TRUE, xname="Cluster", yname="Sample")
pdf(paste0(plotDir, "/sampleByBroadClustBarPlot.pdf"))
stackedBarPlot(clustBySamp, cmap=atac_sample_cmap, namedColors=TRUE, barwidth=barwidth)
dev.off()

# NamedClust by diseaseStatus barplot
clustByDisease <- fractionXbyY(atac_proj$NLabelClust, atac_proj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="Disease")
pdf(paste0(plotDir, "/diseaseByManualLabelsBarPlot.pdf"))
stackedBarPlot(clustByDisease, cmap=disease_cmap, namedColors=TRUE, barwidth=barwidth)
dev.off()

# BroadClust by diseaseStatus barplot
clustByDisease <- fractionXbyY(atac_proj$BroadClust, atac_proj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="Disease")
pdf(paste0(plotDir, "/diseaseByBroadClustBarPlot.pdf"))
stackedBarPlot(clustByDisease, cmap=cmaps_BOR$stallion, namedColors=FALSE, barwidth=barwidth)
dev.off()

##########################################################################################
# Calculate peak-to-gene links using the mapped pseudoRNA
##########################################################################################

atac_proj <- addCoAccessibility(
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI",
    k = 100 # Default is 100
)

atac_proj <- addPeak2GeneLinks( 
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI",
    k = 100 # Default is 100
)

# Plot Peak2Gene heatmap
p <- plotPeak2GeneHeatmap(atac_proj, groupBy="NLabelClust", palGroup=atacLabelClustCmap, k=25)
pdf(paste0(plotDir, "/peakToGeneHeatmap.pdf"), width=16, height=10)
p
dev.off()

# Save intermediate output
saveArchRProject(atac_proj)

##########################################################################################
# Re-identify marker genes using integrated datasets
##########################################################################################

# Marker genes we want to highlight for labeling broad clusters:
markerGenes  <- c(
  "KRT1", "KRT10", "KRT14", "KRT15", # Keratinocytes
  "ITGB8", "SOX9", "LGR5", "LHX2", "KRT75", "CD200", # Hair follicle
  "THY1", "COL1A1", "COL1A2", "COL11A1", # Fibroblasts
  "CD3D", "CD8A", "CD4", "FOXP3", "IKZF2", "CCL5", # T-cells
  "CD14", "CD86", "CD74", "CD163", #Monocytes / macrophages
  "VWF", "PECAM1", "SELE", "SOX17",  # Endothelial
  "ACKR2", "FLT4", "LYVE1", "PROX1",  # Lymphatic endothelial (FLT4 = VEGFR-3)
  "MITF", "TYR", "SOX10", "MLANA", # Melanocyte markers
  "ITGAX", "CD1C", "CD1A", "CLEC9A", # Dendritic/Langerhans cells? (ITGAX = cd11c)
  "FCGR1A", "FCER1A", "ITGAM", # Antigen presentation / Ig interacting 
  "TPM1", "TPM2", "TAGLN", "MYL9", # Muscle
  "ICOS", "CD84", "COL6A3",  "KRT17", "IL10RA", "KRT16", "IKZF1",  # Highly Regulated genes
  "TWIST2", "PRRX1", "POU2F3", "BCL11B", # Highly Regulated TFs
  "TEAD1", "TEAD2", "TEAD3", "TEAD4", "ETS1", "ETS2",
  "RUNX1", "RUNX2", "RUNX3"
) %>% unique()

# Replot GeneScores heatmap with named clusters
markersGS <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "NLabelClust",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# Lists of GeneScores 'marker genes'
markerListGS <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  pal = cmaps_BOR$horizonExtra,
  cutOff = "FDR <= 0.05 & Log2FC >= 1.0", 
  labelMarkers = markerGenes,
  transpose = FALSE
)

draw(heatmapGS, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapGS, name="GeneScores-Marker-Heatmap", width=6, height=10, ArchRProj=atac_proj, addDOC=FALSE)

#Plot the UMAP Embedding with pseudoRNA Marker Genes Overlayed
p <- plotEmbedding(
    ArchRProj = atac_proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes,
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(atac_proj), 
    plotAs="points", size = pointSize
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-GSM_Marker-Genes-W-Imputation.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)


# Now Imputed RNA:

markersPRNA <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "GeneIntegrationMatrix", 
    groupBy = "NLabelClust",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# Lists of pseudo-RNA 'marker genes'
markerListPRNA <- getMarkers(markersPRNA, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")

heatmapPRNA <- plotMarkerHeatmap(
  seMarker = markersPRNA,
  pal = cmaps_BOR$sunrise,
  cutOff = "FDR <= 0.01 & Log2FC >= 2", 
  labelMarkers = markerGenes,
  transpose = FALSE
)

draw(heatmapPRNA, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapPRNA, name="pseudoRNA-Marker-Heatmap", width=6, height=10, ArchRProj=atac_proj, addDOC=FALSE)

#Plot the UMAP Embedding with pseudoRNA Marker Genes Overlayed
p <- plotEmbedding(
    ArchRProj = atac_proj, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes,
    pal = cmaps_BOR$sunrise,
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL, 
    plotAs="points", size = pointSize
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-pseudoRNA_Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)


# Tracks of genes:
p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "NLabelClust", 
    useGroups = unlist(atac.NamedClust)[atacOrder],
    pal = atacLabelClustCmap,
    plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
    sizes = c(8, 0.2, 1.25, 0.5),
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(atac_proj),
    tileSize=250,
    minCells=100
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 6, height = 7)


# Save intermediate output
saveArchRProject(atac_proj)


#############################
# Add motif information
#############################

atac_proj <- addMotifAnnotations(atac_proj, motifSet="cisbp", name="Motif", force=TRUE)

# Add background peaks
atac_proj <- addBgdPeaks(atac_proj, force = TRUE)

atac_proj <- addDeviationsMatrix(
  ArchRProj = atac_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

# Save intermediate output
saveArchRProject(atac_proj)

#############################
# Motif Footprinting
#############################

motifPositions <- getPositions(atac_proj)

motifs <- c(
  "SPIB", "RORA", "NR3C1", "BATF", "AR_", "TP63", "SOX9", "SOX10", "LHX2",
  "BLC11A", "NFIC", "TCF21", "EOMES", "KLF2",
  "TWIST2", "PRRX1", "POU2F3", "RUNX3", "ETS1", "BCL11B", "ZFHX3",
  "TEAD1", "TEAD2", "TEAD3", "TEAD4", "ETS1", "ETS2",
  "RUNX1", "RUNX2", "RUNX3"
  ) %>% unique()
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

# Compute group coverages (Need to recalculate for each 'groupBy' you want to use for footprinting)
atac_proj <- addGroupCoverages(
  ArchRProj=atac_proj, 
  groupBy="NLabelClust", 
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=FALSE
  )

seFoot <- getFootprints(
  ArchRProj = atac_proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "NLabelClust"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = atac_proj, 
  pal = atacLabelClustCmap,
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

# Save intermediate output
saveArchRProject(atac_proj)

##########################################################################################
# Identifying Marker Peaks
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
atac_proj$LFineClust <- unlist(atac.FineClust)[atac_proj$FineClust]

# Exclude certain groups
all_LFineClust <- unique(atac_proj$LFineClust)
exclude <- c("B.cells", "Unknown", "Th_3", "M2.macs_3")
use_groups <- all_LFineClust[all_LFineClust %ni% exclude]

# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "LFineClust",
    useGroups = use_groups,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

# Save marker peaks for later analyses
saveRDS(markerList, file=paste0(wd, "/LFineClust_markerPeakList.rds"))

#Visualize Markers as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  nLabel = 1,
  binaryClusterRows = TRUE,
  transpose = FALSE
)
draw(heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapPeaks, name="Peak-Marker-Heatmap", width=10, height=15, ArchRProj=atac_proj, addDOC=FALSE)

##########################################################################################
# Motif Enrichments
##########################################################################################

#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = atac_proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Rename motifs for more aesthetic plotting:
rownames(enrichMotifs) <- lapply(rownames(enrichMotifs), function(x) strsplit(x, "_")[[1]][1]) %>% unlist()

# Subset to clusters that have at least some enrichment
log10pCut <- 5
keepClust <- colnames(enrichMotifs)[apply(assays(enrichMotifs)$mlog10Padj, 2, max) > log10pCut]

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(enrichMotifs[,keepClust], n=5, transpose=TRUE, cutOff=log10pCut)
draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapEM, name="Motifs-Enriched-Marker-Heatmap", width=12, height=8, ArchRProj=atac_proj, addDOC=FALSE)

# Let's plot a different style heatmap
plot_mat <- plotEnrichHeatmap(enrichMotifs[,keepClust], n=5, transpose=TRUE, 
  cutOff=log10pCut, returnMatrix=TRUE)
plot_mat <- plot_mat[,!grepl("^ENSG", colnames(plot_mat))]

# Get colors for cluster annotation
# (Link FineClusts to BroadClust cmap)
cM <- as.matrix(confusionMatrix(atac_proj$LFineClust, atac_proj$BroadClust))
map_colors <- apply(cM, 1, function(x)colnames(cM)[which.max(x)])
plotColors <- map_colors[rownames(plot_mat)]

bc_to_fc_map <- invertList(plotColors)
bc_to_fc_map <- bc_to_fc_map[c("Tc", "My", "Kc", "Fb", "Mu", "Ve", "Le", "Me")]

col_order <- lapply(bc_to_fc_map, function(bc){
  if(length(bc) == 1){
    bc
  }else{
    colnames(prettyOrderMat(t(plot_mat[bc,]), clusterCols=TRUE, cutOff=1)$mat)
  }
  }) %>% do.call(c,.)

plot_mat <- prettyOrderMat(t(plot_mat[col_order,]), clusterCols=FALSE, cutOff=1)$mat

pdf(paste0(plotDir, "/Custom-Motifs-Enriched-Marker-Heatmap.pdf"), width=15, height=12)
fontsize <- 6
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=plotColors[colnames(plot_mat)],col=list(atac_cluster=broadClustCmap), 
  show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))
hm <- BORHeatmap(
  plot_mat, 
  limits=c(0.0,100.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  top_annotation = ta,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = fontsize),
  column_names_gp = gpar(fontsize = fontsize),
  width = ncol(plot_mat)*unit(0.4, "cm"),
  height = nrow(plot_mat)*unit(0.25, "cm"),
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]",
  border_gp = gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()

##########################################################################################
# ChromVAR UMAPs
##########################################################################################


# Motifs to focus in on:
motifs <- c(
  #"GATA1", "EBF1", "IRF4", "TBX21", "PAX5", # Default immune things
  "TP63", "TCF3", "TCF4", "LHX2", "RBPJ", # Keratinocytes / HFSCs
  "NFIC", "NFIB", "CBFB", "RELA", "KLF4",
  "SOX9", "SOX10", "SOX17", "SOX18", # SOX factors
  "SPIB", "SPI1", "BCL11A", "BCL11B", "RUNX1", "RUNX2", "RUNX3", "NFKB1", "NFKB2", # More immune
  "ELF2", "STAT6", "NFIX", "TCF12", # Endothelial
  "JUND", "FOSL1", "EBF1", # Stress / activation
  "ZNF238", "TCF21", "CEBPA", "AR_", "TBX15", "ARNT", # Fibroblasts
  "TEAD1", "TEAD2", "TEAD3", "TEAD4", "ETS1", "ETS2", # Model hits
  "RUNX1", "RUNX2", "RUNX3"
  ) %>% unique()

markerMotifs <- getFeatures(atac_proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

# Plot over UMAP:
p <- plotEmbedding(
    ArchRProj = atac_proj, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(atac_proj),
    plotAs="points", size = pointSize
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-ChromVar.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)


# End of preprocessing
saveArchRProject(atac_proj)

#############################################################################