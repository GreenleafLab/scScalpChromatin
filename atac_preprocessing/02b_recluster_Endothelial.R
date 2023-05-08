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
source(paste0(scriptPath, "/cluster_labels.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory
subgroup <- "Endothelial"
wd <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/subclustered_%s", subgroup)

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 1.5

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_proj <- readRDS(sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/harmonized_subclustering/%s/%s.rds", subgroup, subgroup))

plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

# Colormaps
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
atac_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(atac_proj$Sample2)] %>% unlist()
disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

rna_sub_cmap <- readRDS(paste0(scriptPath, sprintf("/rna_cmap_%s.rds", subgroup)))

##########################################################################################
# Re-cluster subclustered ArchR project
##########################################################################################

atac_proj <- addIterativeLSI(
    ArchRProj = atac_proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    clusterParams = list(resolution=c(2), sampleCells=30000, maxClusters=6, n.start=10),
    sampleCellsPre = 30000,
    varFeatures = 25000,
    dimsToUse = 1:25,
    force = TRUE
)

# (Batch correct sample effects)
atac_proj <- addHarmony(
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)

# Identify Clusters from Iterative LSI
atac_proj <- addClusters(
    input = atac_proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.2,
    force = TRUE
)

set.seed(1)
atac_proj <- addUMAP(
    ArchRProj = atac_proj, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

# Relabel clusters so they are sorted by cluster size
atac_proj <- relabelClusters(atac_proj)

# Make various cluster plots:
atac_proj <- addImputeWeights(atac_proj)
atac_proj <- visualizeClustering(atac_proj, pointSize=pointSize, sampleCmap=atac_sample_cmap, diseaseCmap=disease_cmap)

# Relabel clusters:
nclust <- length(unique(atac_proj$Clusters))
fineClust <- c(
    "aVe1",
    "aVe2",
    "aVe3",
    "aLe1",
    "aVe4"
    )
names(fineClust) <- atac_proj$Clusters %>% getFreqs() %>% names()
atac_proj$FineClust <- fineClust[atac_proj$Clusters] %>% unname

# Save again and return project
saveArchRProject(atac_proj)

##########################################################################################
# Integrate RNA and ATAC data
##########################################################################################

# integration with scRNA:
set.seed(1)
atac_proj <- addGeneIntegrationMatrix(
    ArchRProj = atac_proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = rna_proj, # Can be a seurat object
    addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
    force = TRUE,
    groupRNA = "FineClust", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this function.
    nameCell = "RNA_paired_cell", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "FineClust_RNA", #Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore" #Name of column where prediction score from scRNA
)
cM <- as.matrix(confusionMatrix(atac_proj$FineClust, atac_proj$FineClust_RNA))

# rna cmap
rna_label_cmap <- rna_sub_cmap
names(rna_label_cmap) <- unlist(rna.FineClust)[names(rna_label_cmap)]

# atac cmap
# (First need to get an RNA to ATAC mapping)
labels_to_rna <- invertList(unlist(rna.FineClust)[names(rna_sub_cmap)])
rna_to_atac <- invertList(atac.FineClust)[names(labels_to_rna)]
names(rna_to_atac) <- unlist(labels_to_rna)
atac_to_rna <- invertList(rna_to_atac) # Removes non-matching clusters

atac_sub_cmap <- rna_sub_cmap[unlist(atac_to_rna)]
names(atac_sub_cmap) <- names(atac_to_rna)

# Re-use rna labels:
remainingColors <- cmaps_BOR$stallion[cmaps_BOR$stallion %ni% unique(atac_sub_cmap, rna_sub_cmap)]
leftout <- fineClust[fineClust %ni% names(atac_sub_cmap)]
catchup <- remainingColors[1:length(leftout)]
names(catchup) <- leftout
atac_sub_cmap <- c(atac_sub_cmap, catchup)

saveRDS(atac_sub_cmap, file=paste0(scriptPath, sprintf("/atac_cmap_%s.rds", subgroup)))

atac_label_cmap <- atac_sub_cmap
names(atac_label_cmap) <- unlist(atac.FineClust)[names(atac_label_cmap)]

# Add labels to project
atac_proj$LabelClust <- unlist(atac.FineClust)[atac_proj$FineClust]
atac_proj$LabelClust_RNA <- unlist(rna.FineClust)[atac_proj$FineClust_RNA]

# Plot confusion matrix heatmap
cM <- as.matrix(confusionMatrix(atac_proj$LabelClust, atac_proj$LabelClust_RNA))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}
cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_sub_cmap <- readRDS(paste0(scriptPath, sprintf("/rna_cmap_%s.rds", subgroup)))
rna_label_cmap <- rna_sub_cmap
names(rna_label_cmap) <- unlist(rna.FineClust)[names(rna_label_cmap)]
atac_sub_cmap <- readRDS(paste0(scriptPath, sprintf("/atac_cmap_%s.rds", subgroup)))
atac_label_cmap <- atac_sub_cmap
names(atac_label_cmap) <- unlist(atac.FineClust)[names(atac_label_cmap)]

rna_order <- c(
  "Lymph.Endo_1",
  "Vas.Endo_1",
  "Vas.Endo_2",
  "Vas.Endo_3",
  "Vas.Endo_4"
)
atac_order <- c(
  "Lymph.Endo",
  "Vas.Endo_1",
  "Vas.Endo_2",
  "Vas.Endo_3",
  "Unknown_2"
)
cM <- cM[atac_order, rna_order]

atac_label_cmap <- atac_label_cmap[names(atac_label_cmap) %in% rownames(cM)]
rna_label_cmap <- rna_label_cmap[names(rna_label_cmap) %in% colnames(cM)]

pdf(paste0(plotDir, sprintf("/ATAC-RNA-integration-cM-heatmap-%s.pdf", subgroup)), width=6, height=6)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <- HeatmapAnnotation(atac_cluster=rownames(cM),col=list(atac_cluster=atac_label_cmap), 
  which="row", show_legend=c("atac_cluster"=FALSE))
ta <- HeatmapAnnotation(rna_cluster=colnames(cM),col=list(rna_cluster=rna_label_cmap), 
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
  border_gp=gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()


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
  groupBy = "LabelClust", 
  colorBy = "colData", 
  name = "predictedScore",
  ylim = pScoreLim,
  pal = atac_label_cmap
)

plotList[[3]] <- plotGroups(ArchRProj = atac_proj, 
  groupBy = "LabelClust_RNA", 
  colorBy = "colData", 
  name = "predictedScore",
  ylim = pScoreLim,
  pal = rna_label_cmap
)
plotPDF(plotList = plotList, name = "prediction_accuracy_by_group", width = 4, height = 4,  ArchRProj = atac_proj, addDOC = FALSE)


# Plot the UMAPs by Sample and Cluster:
p1 <- plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Sample2", embedding = "UMAP", labelMeans = FALSE,
  size = pointSize, plotAs = "points", pal=atac_sample_cmap)
p2 <- plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "diseaseStatus", embedding = "UMAP", labelMeans = FALSE,
  size = pointSize, plotAs = "points", pal=disease_cmap)
p3 <- plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "LabelClust", embedding = "UMAP", labelMeans = FALSE,
  size = pointSize, plotAs = "points", pal=atac_label_cmap)
p4 <- plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "LabelClust_RNA", embedding = "UMAP", labelMeans = FALSE,
  size = pointSize, plotAs = "points", pal=rna_label_cmap)
p5 <- plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "predictedScore", embedding = "UMAP", labelMeans = FALSE,
  size = pointSize, plotAs = "points", colorLimits=c(0.0,1.0))

ggAlignPlots(p1, p2, p3, p4, p5, type = "h")
plotPDF(p1, p2, p3, p4, p5,
    name = "Plot-UMAP-ATACvRNAclusters.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)


# Plot scATAC with ATAC cluster labels
umapPlots <- list()
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="LabelClust")
umapPlots[["ATAC_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=atac_label_cmap, namedColors=TRUE, point_size=pointSize, covarLabel="ATAC_labels")
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), unlist(rna.FineClust)[rna_proj$FineClust])
# Randomize cells before plotting UMAP
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
umapPlots[["RNA_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=rna_label_cmap, namedColors=TRUE, point_size=pointSize, covarLabel="RNA_labels")
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="LabelClust_RNA")
umapPlots[["RNA_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=rna_label_cmap, namedColors=TRUE, point_size=pointSize, covarLabel="RNA_labels")
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="predictedScore", lowerPctLim=0.01)
umapPlots[["predScore"]] <- plotUMAP(umapDF, dataType = "quantitative", cmap = cmaps_BOR$solarExtra, point_size=pointSize, covarLabel="predictedScore")

pdf(paste0(plotDir,sprintf("/clusters_UMAP_%s.pdf", subgroup)), width=10, height=8)
umapPlots
dev.off()

clustBySamp <- fractionXbyY(atac_proj$FineClust, atac_proj$Sample2, add_total=TRUE, xname="FineClust", yname="Sample2")

pdf(paste0(plotDir, "/sampleByClustBarPlot.pdf"))
stackedBarPlot(clustBySamp, cmap = atac_sample_cmap, namedColors=TRUE)
dev.off()

clustByDisease <- fractionXbyY(atac_proj$Clusters, atac_proj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="Disease")

pdf(paste0(plotDir, "/diseaseByClustBarPlot.pdf"))
stackedBarPlot(clustByDisease, cmap = disease_cmap, namedColors=TRUE)
dev.off()

mClustBySample <- fractionXbyY(factor(atac_proj$LabelClust, levels=names(getFreqs(atac_proj$LabelClust)), ordered=TRUE), 
  atac_proj$Sample2, add_total=TRUE, xname="Cluster", yname="Sample2")

pdf(paste0(plotDir, "/sampleByLabeledClustBarPlot.pdf"))
stackedBarPlot(mClustBySample, cmap = atac_sample_cmap, namedColors=TRUE)
dev.off()

mClustByDisease <- fractionXbyY(factor(atac_proj$LabelClust, levels=names(getFreqs(atac_proj$LabelClust)), ordered=TRUE), 
  atac_proj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="Disease")

pdf(paste0(plotDir, "/diseaseByManualClustBarPlot.pdf"))
stackedBarPlot(mClustByDisease, cmap = disease_cmap, namedColors=TRUE)
dev.off()

# Save again and return project
saveArchRProject(atac_proj)

###########################################################################################