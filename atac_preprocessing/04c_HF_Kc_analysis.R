#!/usr/bin/env Rscript

#######################################################################################
# Final analyses on keratinocytes
# (This assumes you have already performed peak assignments, RNA integration, etc.)
#######################################################################################

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

# set working directory
subgroup <- "HF_Kc"
wd <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/subclustered_%s", subgroup)
full_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/fine_clustered"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 2

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

# Load labels from file
source(paste0(scriptPath, "/cluster_labels.R"))

##########################################################################################
# Reassign peak matrix
##########################################################################################

# Get peaks that were called on this subproject's subclusters from full ArchR project
full_proj <- loadArchRProject(full_dir, force=TRUE)
full_peaks <- getPeakSet(full_proj)
peaks <- getClusterPeaks(full_proj, clusterNames=unique(atac_proj$FineClust), peakGR=full_peaks)
rm(full_proj); gc() # Cleanup big project

# Now add these peaks to the subproject and generate peak matrix
atac_proj <- addPeakSet(atac_proj, peakSet=peaks, force=TRUE)
atac_proj <- addPeakMatrix(atac_proj, force=TRUE)
atac_proj <- addMotifAnnotations(atac_proj, motifSet="cisbp", name="Motif", force=TRUE)

##########################################################################################
# Re-cluster subclustered ArchR project
##########################################################################################

atac_proj <- addIterativeLSI(
    ArchRProj = atac_proj,
    useMatrix = "PeakMatrix", 
    name = "IterativeLSI", 
    clusterParams = list(resolution=c(2), sampleCells=30000, maxClusters=6, n.start=10), 
    sampleCellsPre = 30000,
    varFeatures = 25000, # Default is 25000
    dimsToUse = 1:30, # Detault is 1:30
    force = TRUE
)

# Identify Clusters from Iterative LSI
atac_proj <- addClusters(
    input = atac_proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4, 
    force = TRUE
)

set.seed(1)
atac_proj <- addUMAP(
    ArchRProj = atac_proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 20, 
    minDist = 0.1, 
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
fineClust <- sapply(1:nclust, function(x) paste0("aHF", x))
names(fineClust) <- atac_proj$Clusters %>% getFreqs() %>% names()
atac_proj$HFClust <- fineClust[atac_proj$Clusters] %>% unname

# Ensure Impute Weights are up to date
atac_proj <- addImputeWeights(atac_proj)

# Save again and return project
saveArchRProject(atac_proj)

##########################################################################################
# (Re)Integrate RNA and ATAC data
##########################################################################################

set.seed(1)
atac_proj <- addGeneIntegrationMatrix(
    ArchRProj = atac_proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = rna_proj, # Can be a seurat object
    addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
    force = TRUE,
    groupRNA = "HFClust", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this function.
    nameCell = "HF_RNA_paired_cell", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "HFClust_RNA", #Name of column where group from scRNA is matched to each cell
    nameScore = "HFpredictedScore" #Name of column where prediction score from scRNA
)
cM <- as.matrix(confusionMatrix(atac_proj$HFClust, atac_proj$HFClust_RNA))

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
remainingColors <- cmaps_BOR$stallion[cmaps_BOR$stallion %ni% unique(c(atac_sub_cmap, rna_sub_cmap))]
leftout <- fineClust[fineClust %ni% names(atac_sub_cmap)]
catchup <- remainingColors[1:length(leftout)]
names(catchup) <- leftout
atac_sub_cmap <- c(atac_sub_cmap, catchup)

saveRDS(atac_sub_cmap, file=paste0(scriptPath, sprintf("/atac_cmap_%s.rds", subgroup)))

atac_label_cmap <- atac_sub_cmap
names(atac_label_cmap) <- unlist(atac.FineClust)[names(atac_label_cmap)]

# Add labels to project
atac_proj$LabelClust <- unlist(atac.FineClust)[atac_proj$HFClust]
atac_proj$LabelClust_RNA <- unlist(rna.FineClust)[atac_proj$HFClust_RNA]

cM <- as.matrix(confusionMatrix(atac_proj$LabelClust, atac_proj$LabelClust_RNA))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}
cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

pdf(paste0(plotDir, "/ATAC-RNA-integration-cM-heatmap.pdf"), width=6, height=6)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  cM, 
  limits=c(0,1), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$whitePurple,
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
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), unlist(rna.FineClust)[rna_proj$HFClust])
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

##########################################################################################
# Re-compute coaccessibility, P2G links, deviations, etc.
##########################################################################################

# Calculate coaccessibility
atac_proj <- addCoAccessibility(
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI"
)

# Calculate peak-to-gene links using the broad cluster mapped pseudoRNA
atac_proj <- addPeak2GeneLinks(
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI"
)

# Add background peaks
atac_proj <- addBgdPeaks(atac_proj, force = TRUE)

atac_proj <- addDeviationsMatrix(
  ArchRProj = atac_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(atac_proj)

####################################################################################
# Marker Genes
####################################################################################

allGenes <- getGenes(atac_proj)$symbol %>% unname()
allGenes <- allGenes[!is.na(allGenes)]

fzd_genes <- allGenes[grepl("FZD", allGenes)]
fzd_genes <- fzd_genes[fzd_genes %ni% c("FZD10-DT")] #(Not present in GeneIntegrationMatrix)
wnt_genes <- allGenes[grepl("WNT", allGenes)]

markerGenes  <- c(
  # https://www.proteinatlas.org/humanproteome/tissue/skin
  "KRT14", "KRT5", "KRT15", "COL17A1", # Basal epithelia
  "ITGA6", "ITGB1", "CD200", "LGR5","LHX2", "FRZB", "FZD1", "FZD5", "FZD10",  "IL31RA", "OSMR", # HFSCs?
  "CD34", "CDH3", "LGR5", "LGR6", "CDKN2A", "RUNX1", # Hair germ markers (some say HG is CD34 neg?)
  "KRT81", "KRT82", "KRT83", "HOXC13", "LEF1", "WNT3", # Matrix hair keratins/genes
  "KRT75", # ORS keratins / genes
  "ELOVL3", # Sebaceous 
  "AXIN2", "DKK3", "DKK1", "DKK2", "SFRP1", # HFSC WNT signaling??
  "SOX9", "LHX2", "NFATC1", "TCF3", "TCF4", "IRF1", # Key HFSC TFs
  fzd_genes,
  wnt_genes,
  "DVL1", "DVL2", "DVL3", "AXIN1", "AXIN2", "CTNNB1" # CTNNB1 = beta-catenin
) %>% unique()

# GeneScores first:


# Identify Marker Gene through Pairwise Test vs Bias-Matched Background:
markersGS <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "HFClust",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersGS)$name

# Lists of 'marker genes'
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.50")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = FALSE # Can't transpose if we want pretty ordered heatmaps...
)

draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width=6, height=8, ArchRProj=atac_proj, addDOC=FALSE)


p <- plotEmbedding(
    ArchRProj = atac_proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(atac_proj), 
    size = pointSize, plotAs = "points"
)

plotPDF(plotList = p, 
    name = "HF-Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)

# Now Imputed RNA:

markersPRNA <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "GeneIntegrationMatrix", 
    groupBy = "HFClust",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersPRNA)$name

# Lists of pseudo-RNA 'marker genes'
markerList <- getMarkers(markersPRNA, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")

heatmapPRNA <- plotMarkerHeatmap(
  seMarker = markersPRNA, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.5", 
  labelMarkers = markerGenes,
  transpose = FALSE # Can't transpose if we want pretty ordered heatmaps...
)

draw(heatmapPRNA, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPRNA, name = "pseudoRNA-Marker-Heatmap", width = 6, height = 8, ArchRProj = atac_proj, addDOC = FALSE)


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
    name = "HF-Plot-UMAP-pseudoRNA_Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)


atac_sub_cmap <- readRDS(file=paste0(scriptPath, sprintf("/atac_cmap_%s.rds", subgroup)))

# Tracks of genes:
p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "HFClust", 
    useGroups = c("aHF6", "aHF3", "aHF2", "aHF1", "aHF5", "aHF4"),
    pal = atac_sub_cmap,
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

# Plot ChromVAR deviations
motifPositions <- getPositions(atac_proj)

plotVarDev <- getVarDeviations(atac_proj, plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = atac_proj, addDOC = FALSE)

motifs <- c(
  "AR_", "TP63", "SOX4", "SOX9", "JUND", "FOS", "BATF", "FOXC1", "LHX2", "LEF1",
  "PRRX1", "PRRX2", "POU2F3", "RUNX3", "CEBPB", "HOXA1", "HOXC13", # Highly Regulated TFs
  "NFIC", "KLF4", "NFATC1", "TCF3", "TCF4", "IRF1", "EGR1", "EGR2", "EGR3", "EGR4"
  )
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
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
    name = "HF-Plot-UMAP-ChromVar.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)

#saveArchRProject(atac_proj)

####################################################################################
# Milo on scATAC
####################################################################################

library(miloR)
library(SingleCellExperiment)

# Convert ArchR project into a singleCellExperiment

# Subset Milo object to remove samples that are too lowly represented:
ccd <- atac_proj@cellColData
samp_freqs <- getFreqs(ccd$Sample2)
valid_samps <- samp_freqs[samp_freqs > 10] %>% names()
valid_cells <- rownames(ccd[ccd$Sample2 %in% valid_samps,])

# First get just the peak matrix (what we will use as our 'counts' matrix)
counts_mat <- getMatrixFromProject(atac_proj, useMatrix="PeakMatrix") # This doesn't return a matrix with identifyable peak names...
counts_mat <- counts_mat[,valid_cells]

sce_atac <- SingleCellExperiment(counts_mat)
names(assays(sce_atac)) <- names(assays(counts_mat)) # Doesn't carry over assay names in construction for some reason
dim_reduc <- atac_proj@reducedDims$IterativeLSI$matSVD
umap <- atac_proj@embeddings$UMAP$df
reducedDims(sce_atac)[["LSI"]] <- dim_reduc[valid_cells,]
reducedDims(sce_atac)[["UMAP"]] <- umap[valid_cells,]
colData(sce_atac) <- colData(counts_mat)

# Compare AA to C (all controls)
sce_atac$diseaseStatus <- factor(ifelse(sce_atac$diseaseStatus == "AA", "AA", "C"), levels=c("C", "AA"), ordered=TRUE)

# Create Milo object
set.seed(123)
milo_proj <- Milo(sce_atac)

# build KNN graph using LSI reduced dimensions
d <- ncol(reducedDims(sce_atac)$LSI)
k <- 30 # 30 is default
reduced.dims <- "LSI" # Milo capitalizes the dim reduc names
milo_proj <- buildGraph(milo_proj, k=k, d=d, reduced.dim=reduced.dims)

# Defining representative neighborhoods on the KNN graph
prop <- 0.3
milo_proj <- makeNhoods(milo_proj, prop=prop, k=k, d=d, refined=TRUE, reduced_dims=reduced.dims)

# "As a rule of thumb we want to have an average neighbourhood size over 5 x N_samples"
pdf(paste0(plotDir, "/miloR_nhood_size_hist.pdf"), width=8, height=6)
plotNhoodSizeHist(milo_proj)
dev.off()

# Counting cells in neighborhoods
milo_proj <- countCells(milo_proj, meta.data=as.data.frame(colData(milo_proj)), sample="Sample2")

# Defining experimental design
milo_design <- data.frame(colData(milo_proj))[,c("Sample2", "diseaseStatus")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$Sample

# Computing neighborhood connectivity
milo_proj <- calcNhoodDistance(milo_proj, d=d, reduced.dim=reduced.dims)

# Testing differential abundance
da_results <- testNhoods(milo_proj, design = ~ diseaseStatus, design.df=milo_design, reduced.dim=reduced.dims)

# Visualize results:
milo_proj <- buildNhoodGraph(milo_proj)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo_proj, da_results, layout="UMAP", 
    alpha=0.1, size_range=c(1, 4), node_stroke=0.1) 

pdf(paste0(plotDir, sprintf("/miloR_atac_DA_UMAP_%s.pdf", subgroup)), width=7, height=6)
nh_graph_pl
dev.off()

# Save differential abundance results

# First need to map Milo neighborhoods to majority cluster IDs
nhoodmat <- milo_proj@nhoods
cell_to_clust <- atac_proj$HFClust
names(cell_to_clust) <- atac_proj$cellNames

# For each nhood, identify the most common cluster
nhood_to_clust <- apply(nhoodmat, 2, function(x){
  cellnames <- names(x[x>0]) # Which cells are part of this nhood
  names(sort(table(cell_to_clust[cellnames]), decreasing=TRUE))[1] # Most frequent cluster
  })

# Assign majority cluster to nhoods
da_results$majority_cluster <- nhood_to_clust[da_results$Nhood]
da_results$majority_lclust <- unlist(atac.FineClust)[da_results$majority_cluster]
da_results <- da_results[order(da_results$PValue, decreasing=FALSE),]

# Save table
table_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/supplemental_tables"
write.table(da_results, file=paste0(table_dir, "/HF_subclustered_Milo_results.tsv"), quote=FALSE, sep="\t", col.names=NA, row.names=TRUE) 

####################################################################################
# Hair follicle trajectories
####################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)

# trajectory (aHF6 is stem, aHF3 is terminal; don't include aHF5)
atac_proj <- addSlingShotTrajectories(
  ArchRProj = atac_proj,
  name = "sheath",
  useGroups = c("aHF1", "aHF2", "aHF3", "aHF4", "aHF6"),
  principalGroup = "aHF6",
  groupBy = "HFClust",
  embedding = "UMAP",
  force = TRUE,
  seed = 1
)

# Plot trajectories
p1 <- plotTrajectory(
  atac_proj, 
  trajectory = "sheath.Curve1", 
  embedding = "UMAP",
  colorBy = "cellColData", 
  name = "sheath.Curve1", 
  continuousSet = "fireworks2",
  plotAs = "points",
  smoothWindow = 8,
  size = pointSize
  )

plotPDF(p1, name="HF_slingshot_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=5, height=5)

####################################################################################
# More trajectory plots
####################################################################################

# Get list of genes we want to highlight (e.g. genes involved in HF development)
library(org.Hs.eg.db)
library(GO.db)
go_id = GOID(GOTERM[Term(GOTERM) == "hair follicle development"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
hfdev_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

label_genes <- c(hfdev_genes, markerGenes) %>% unique() %>% sort()

# ChromVAR motifs:
sheath.trajMM  <- getTrajectory(ArchRProj=atac_proj, name="sheath.Curve1", 
    useMatrix="MotifMatrix", log2Norm=FALSE, groupEvery=2.0)

var.cutoff <- 0.925 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)

p1 <- plotTrajectoryHeatmap(sheath.trajMM, pal=paletteContinuous(set="solarExtra"),
    varCutOff=var.cutoff, limits=lims)

plotPDF(p1, name="HF_MM_trajectories_shaft.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=7, height=10)

# GeneScoreMatrix :

sheath.trajGSM  <- getTrajectory(ArchRProj=atac_proj, name="sheath.Curve1", 
    useMatrix="GeneScoreMatrix", log2Norm=TRUE, groupEvery=2.0)

rownames(sheath.trajGSM) <- rownames(sheath.trajGSM) %>% strsplit(":") %>% sapply('[',2)

var.cutoff <- 0.9 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)
labeltop <- 25 # Default is 50

p1 <- plotTrajectoryHeatmap(sheath.trajGSM, pal=paletteContinuous(set="horizonExtra"),
    varCutOff=var.cutoff, limits=lims, labelMarkers=label_genes, labelTop=labeltop)

plotPDF(p1, p2, name="HF_GSM_trajectories_shaft.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=7, height=10)

# GeneIntegrationMatrix :
sheath.trajGIM  <- getTrajectory(ArchRProj=atac_proj, name="sheath.Curve1", 
    useMatrix="GeneIntegrationMatrix", log2Norm=TRUE, groupEvery=2.0)

rownames(sheath.trajGIM) <- rownames(sheath.trajGIM) %>% strsplit(":") %>% sapply('[',2)

var.cutoff <- 0.925 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)

p1 <- plotTrajectoryHeatmap(sheath.trajGIM, pal=cmaps_BOR$sunrise,
    varCutOff=var.cutoff, limits=lims, labelMarkers=label_genes, labelTop=labeltop)

plotPDF(p1, name="HF_GIM_trajectories_shaft.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=7, height=10)

####################################################################################
# Integrative pseudo-time analysis
####################################################################################

# GeneScoreMatrix:
corGSM_MM <- correlateTrajectories(shaft.trajGSM, shaft.trajMM)

trajGSM2 <- shaft.trajGSM[corGSM_MM[[1]]$matchname1, ]
trajMM2 <- shaft.trajMM[corGSM_MM[[1]]$name2, ]

trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat=TRUE, varCutOff=0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal=paletteContinuous(set="horizonExtra"),  varCutOff=0, rowOrder=rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal=paletteContinuous(set="solarExtra"), varCutOff=0, rowOrder=rowOrder)
plotPDF(ht1 + ht2, name="shaft_GSM_MM_corr_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=12, height=7)

# GeneIntegrationMatrix:
corGIM_MM <- correlateTrajectories(shaft.trajGIM, shaft.trajMM)

trajGIM2 <- shaft.trajGIM[corGIM_MM[[1]]$matchname1, ]
trajMM2 <- shaft.trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat=TRUE, varCutOff=0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))

ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal=cmaps_BOR$sunrise,  varCutOff=0, rowOrder=rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal=paletteContinuous(set="solarExtra"), varCutOff=0, rowOrder=rowOrder)
plotPDF(ht1 + ht2, name="shaftHG_GIM_MM_corr_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=12, height=7)

# Save project
saveArchRProject(atac_proj)

####################################################################################

# Plot WNT / FZD 'waves'
shaft.trajGIM  <- getTrajectory(ArchRProj=atac_proj, name="shaft.Curve1", 
    useMatrix="GeneIntegrationMatrix", log2Norm=TRUE, groupEvery=3.0)

shaft.trajGSM  <- getTrajectory(ArchRProj=atac_proj, name="shaft.Curve1", 
    useMatrix="GeneScoreMatrix", log2Norm=TRUE, groupEvery=3.0)

rownames(shaft.trajGIM) <- rownames(shaft.trajGIM) %>% strsplit(":") %>% sapply('[',2)
rownames(shaft.trajGSM) <- rownames(shaft.trajGSM) %>% strsplit(":") %>% sapply('[',2)

# Exclude ribosomal genes
blacklist <- c(
  grep(pattern="^RPS", x=rownames(shaft.trajGIM), value=TRUE),
  grep(pattern="^RPL", x=rownames(shaft.trajGIM), value=TRUE)
)
shaft.trajGIM <- shaft.trajGIM[rownames(shaft.trajGIM) %ni% blacklist,]
shaft.trajGSM <- shaft.trajGSM[rownames(shaft.trajGSM) %ni% blacklist,]


# Pre-filter genes by variance cutoff
wnt.pathway <- c(
  wnt_genes, # WNT effectors
  fzd_genes, # Frizzled receptors
  "CTNNB1", # beta-catenin
  "TCF3", "TCF4", "LEF1", # WNT-responsive TFs
  "LRP1", "LRP4", "LRP5", "LRP6",
  "DKK1", "DKK2", "DKK3", "DKK4", "DKKL1", "FRZB", "SFRP1", "SFRP2", "SFRP3", "SFRP4", "SFRP5", # WNT pathway inhibitors
  "DVL1", "DVL2", "DVL3", "AXIN1", "AXIN2" # WNT degradation complex components
  )

# Or using GO term:
go_id = GOID(GOTERM[Term(GOTERM) == "cell-cell signaling by wnt"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
wnt.GO.genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

varCutoff <- 0.5 # Cutoff for variable genes
varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(assays(shaft.trajGIM)$smoothMat))
var.traj.genes <- rownames(shaft.trajGIM[order(varQ, decreasing=TRUE),]) %>% head((1-varCutoff)*nrow(shaft.trajGIM))

wnt.trajGIM <- shaft.trajGIM[var.traj.genes[var.traj.genes %in% wnt.pathway],]
wnt.trajGSM <- shaft.trajGSM[var.traj.genes[var.traj.genes %in% wnt.pathway],]

var.cutoff <- 0 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)
#exclude <- c("VAX", "NFE2")
p1 <- plotTrajectoryHeatmap(wnt.trajGIM, pal=cmaps_BOR$sunrise, 
    varCutOff=var.cutoff, limits=lims, labelTop=50)

plotPDF(p1, name="HF_GIM_shaft_WNTs.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=5, height=7)

var.cutoff <- 0 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)
#exclude <- c("VAX", "NFE2")
p1 <- plotTrajectoryHeatmap(wnt.trajGSM, pal=cmaps_BOR$horizonExtra, 
    varCutOff=var.cutoff, limits=lims, labelTop=50)

plotPDF(p1, name="HF_GSM_shaft_WNTs.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=5, height=7)


# GO terms on trajectory variable genes
source(paste0(scriptPath, "/GO_wrappers.R"))

varCutoff <- 0.9 # Cutoff for variable genes
varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(assays(shaft.trajGIM)$smoothMat))
var.traj.genes <- rownames(shaft.trajGIM[order(varQ, decreasing=TRUE),]) %>% head((1-varCutoff)*nrow(shaft.trajGIM))

upGO <- rbind(
  calcTopGo(rownames(shaft.trajGIM), interestingGenes=var.traj.genes, nodeSize=5, ontology="BP") 
  #calcTopGo(rownames(shaft.trajGIM), interestingGenes=var.traj.genes, nodeSize=5, ontology="MF")
  )
upGO <- upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
up_go_plot <- topGObarPlot(upGO, cmap=cmaps_BOR$comet, nterms=10, border_color="black", 
  barwidth=0.9, title="variable HF trajectory terms", enrichLimits=c(0.0, 3))

pdf(paste0(plotDir, "/HF_shaft_varGene_topGO.pdf"), width=8, height=6)
print(up_go_plot)
dev.off()


####################################################################################
# Explore differences in WNT dynamics between control and AA samples
####################################################################################

