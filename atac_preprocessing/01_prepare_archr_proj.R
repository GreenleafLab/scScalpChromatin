#!/usr/bin/env Rscript

#####################################################################
# Build ArchR project and perform basic pre-processing and subsetting
#####################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(mclust)
})

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin/"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/sample_metadata.R"))

# Set Threads to be used
addArchRThreads(threads = 10)

# set working directory
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/baseline_preprocessing"

# color palettes
sample_cmap <- readRDS("/home/users/boberrey/git_clones/hairATAC/sample_cmap.rds")

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
pointSize <- 0.5

##########################################################################################
# Preparing Data
##########################################################################################

#Get scalp input files
inputFiles <- getInputFiles("/oak/stanford/groups/wjg/boberrey/hairATAC/scATAC/fragment_files")

# Create Arrow Files (~30 minutes)
# recommend you use as many threads as samples.
# This step will for each sample :
# 1. Read Accessible Fragments.
# 2. Identify Cells QC Information (TSS Enrichment, Nucleosome info).
# 3. Filter Cells based on QC parameters.
# 4. Create a TileMatrix 500-bp genome-wide.
# 5. Create a GeneScoreMatrix.
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  geneAnno = geneAnno,
  genomeAnno = genomeAnno,
  minTSS = 0, # Don't filter at this point
  minFrags = 1000, # Default is 1000.
  addTileMat = FALSE, # Don't add tile or geneScore matrices yet. Will add them after we filter
  addGeneScoreMat = FALSE
)

#Create ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno,
  outputDirectory = "unfiltered_output"
)

# Remove Arrow files after copying
unlink(paste0(wd, "/*.arrow"))

# Now, identify likely cells:
identifyCells <- function(df, TSS_cutoff=6, nFrags_cutoff=2000, minTSS=5, minFrags=1000, maxG=4){
    # Identify likely cells based on gaussian mixture modelling.
    # Assumes that cells, chromatin debris, and other contaminants are derived from
    # distinct gaussians in the TSS x log10 nFrags space. Fit a mixture model to each sample
    # and retain only cells that are derived from a population with mean TSS and nFrags passing
    # cutoffs
    ####################################################################
    # df = data.frame of a single sample with columns of log10nFrags and TSSEnrichment
    # TSS_cutoff = the TSS cutoff that the mean of a generating gaussian must exceed
    # nFrags_cutoff = the log10nFrags cutoff that the mean of a generating gaussian must exceed
    # minTSS = a hard cutoff of minimum TSS for keeping cells, regardless of their generating gaussian
    # maxG = maximum number of generating gaussians allowed

    cellLabel <- "cell"
    notCellLabel <- "not_cell"

    if(nFrags_cutoff > 100){
        nFrags_cutoff <- log10(nFrags_cutoff)
        minFrags <- log10(minFrags)
    } 
    
    # Fit model
    set.seed(1)
    mod <- Mclust(df, G=2:maxG, modelNames="VVV")

    # Identify classifications that are likely cells
    means <- mod$parameters$mean

    # Identify the gaussian with the maximum TSS cutoff
    idents <- rep(notCellLabel, ncol(means))
    idents[which.max(means["TSSEnrichment",])] <- cellLabel

    names(idents) <- 1:ncol(means)

    # Now return classifications and uncertainties
    df$classification <- idents[mod$classification]
    df$classification[df$TSSEnrichment < minTSS] <- notCellLabel
    df$classification[df$nFrags < minFrags] <- notCellLabel
    df$cell_uncertainty <- NA
    df$cell_uncertainty[df$classification == cellLabel] <- mod$uncertainty[df$classification == cellLabel]
    return(list(results=df, model=mod))
}

# Run classification on all samples
minTSS <- 5
samples <- unique(proj$Sample)
cellData <- getCellColData(proj)
cellResults <- lapply(samples, function(x){
  df <- cellData[cellData$Sample == x,c("nFrags","TSSEnrichment")]
  df$log10nFrags <- log10(df$nFrags)
  df <- df[,c("log10nFrags","TSSEnrichment")]
  identifyCells(df, minTSS=minTSS)
  })
names(cellResults) <- samples

# Save models for future reference
saveRDS(cellResults, file = paste0(proj@projectMetadata$outputDirectory, "/cellFiltering.rds"))

# Plot filtering results
for(samp in samples){
    df <- as.data.frame(cellResults[[samp]]$results)
    cell_df <- df[df$classification == "cell",]
    non_cell_df <- df[df$classification != "cell",]

    xlims <- c(log10(500), log10(100000))
    ylims <- c(0, 18)
    # QC Fragments by TSS plot w/ filtered cells removed:
    p <- ggPoint(
        x = cell_df[,1], 
        y = cell_df[,2], 
        size = 1.5,
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = xlims,
        ylim = ylims,
        title = sprintf("%s droplets plotted", nrow(cell_df)),
        rastr = TRUE
    )
    # Add grey dots for non-cells
    p <- p + geom_point_rast(data=non_cell_df, aes(x=log10nFrags, y=TSSEnrichment), color="light grey", size=0.5)
    p <- p + geom_hline(yintercept = minTSS, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
    plotPDF(p, name = paste0(samp,"_EM_model_filtered_cells_TSS-vs-Frags.pdf"), ArchRProj = proj, addDOC = FALSE)
}


# Now, filter ATAC project to contain only cells
finalCellCalls <- lapply(cellResults, function(x) x$results) %>% do.call(rbind, .)
proj <- addCellColData(proj, data=finalCellCalls$classification, name="cellCall", cells=rownames(finalCellCalls), force=TRUE)
proj <- addCellColData(proj, data=finalCellCalls$cell_uncertainty, name="cellCallUncertainty", cells=rownames(finalCellCalls), force=TRUE)

# Add Demuxlet results to C_SD_POOL sample
demuxResults <- c("/oak/stanford/groups/wjg/boberrey/hairATAC/bulkATAC/20200108_C_SD_Demux/demuxlet_output/C_SD_POOL.best")
proj <- addDemuxletResults(proj, bestFiles=demuxResults, sampleNames="C_SD_POOL")

# Relabel demuxlet samples to match existing sample formatting
demuxConvert <- c(
    "C_SD_01_S15" = "C_SD4",
    "C_SD_06_S16" = "C_SD5",
    "C_SD_08_S17" = "C_SD6",
    "C_SD_10_S18" = "C_SD7"
    )
proj$Sample2 <- ifelse(proj$DemuxletBest == "NotClassified", proj$Sample, demuxConvert[proj$DemuxletBest])

# Real cells pass QC filter and for C_SD_POOL are classified singlets
realCells <- getCellNames(proj)[(proj$cellCall == "cell") & (proj$DemuxletClassify %ni% c("AMB", "DBL")) & (proj$Sample2 != "C_SD_POOL")]
subProj <- subsetArchRProject(proj, cells=realCells, 
    outputDirectory="filtered_output", dropCells=TRUE, force=TRUE)

# Add sample metadata
subProj$preservation <- samp.preservation[subProj$Sample2] %>% unlist() %>% as.factor()
subProj$sex <- samp.sex[subProj$Sample2] %>% unlist() %>% as.factor()
subProj$age <- samp.age[subProj$Sample2] %>% unlist()

# Now, add tile matrix and gene score matrix to ArchR project
subProj <- addTileMatrix(subProj, force=TRUE)
subProj <- addGeneScoreMatrix(subProj, force=TRUE)

# Add Infered Doublet Scores to ArchR project (~5-10 minutes)
subProj <- addDoubletScores(subProj, dimsToUse=1:20, scaleDims=TRUE, LSIMethod=2)

sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(subProj$Sample2)]
samp_cmap <- unlist(sample_cmap)

# Visualize numeric metadata per grouping with a violin plot now that we have created an ArchR Project.
plotList <- list()
plotList[[1]] <- plotGroups(ArchRProj = subProj, 
  groupBy = "Sample2", 
  colorBy = "colData", 
  name = "TSSEnrichment",
  pal=samp_cmap
)
plotList[[2]] <- plotGroups(ArchRProj = subProj, 
  groupBy = "Sample2", 
  colorBy = "colData", 
  name = "DoubletEnrichment",
  pal=samp_cmap
)
plotPDF(plotList = plotList, name = "TSS-Doublet-Enrichment", width = 4, height = 4,  ArchRProj = subProj, addDOC = FALSE)

# Filter doublets:
subProj <- filterDoublets(subProj, filterRatio = 1)

# Save filtered ArchR project
saveArchRProject(subProj)

##########################################################################################
# Reduced Dimensions and Clustering
##########################################################################################

# Add info about whether cell is 'control' or 'disease' (C_SD or AA)
subProj$diseaseStatus <- NA
subProj$diseaseStatus <- ifelse(grepl("C_SD", subProj$Sample), "C_SD", subProj$diseaseStatus)
subProj$diseaseStatus <- ifelse(grepl("C_PB", subProj$Sample), "C_PB", subProj$diseaseStatus)
subProj$diseaseStatus <- ifelse(grepl("AA", subProj$Sample), "AA", subProj$diseaseStatus)
disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

# Reduce Dimensions with Iterative LSI (<5 minutes)
set.seed(1)
subProj <- addIterativeLSI(
    ArchRProj = subProj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    sampleCellsPre = 15000,
    varFeatures = 50000, 
    dimsToUse = 1:25,
    force = TRUE
)

# Identify Clusters from Iterative LSI
subProj <- addClusters(
    input = subProj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.6,
    force = TRUE
)

##########################################################################################
# Visualize Data
##########################################################################################

set.seed(1)
subProj <- addUMAP(
    ArchRProj = subProj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 50, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

# Relabel clusters so they are sorted by cluster size
subProj <- relabelClusters(subProj)

subProj <- addImputeWeights(subProj)

# Make various cluster plots:
subProj <- visualizeClustering(subProj, pointSize=pointSize, sampleCmap=samp_cmap, diseaseCmap=disease_cmap)

# Save unfiltered ArchR project
saveArchRProject(subProj)

##########################################################################################
# Remove doublet clusters
##########################################################################################

# There are a few clusters that appear to be mostly doublets / poor quality cells
# (Poor TSS enrichment, few marker peaks / GS in downstream analysis, higher cellCallUncertainty)
# and were not filtered by automated doublet removal or by basic QC filtering
# They will be manually filtered here.
proj <- loadArchRProject(paste0(wd, "/filtered_output/"), force=TRUE)

# Identify Marker Gene through Pairwise Test vs Bias-Matched Background:
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersGS)$name

# Lists of 'marker genes'
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")

# Marker genes we want to highlight for labeling broad clusters:
markerGenes  <- c(
  "KRT1", "KRT10", "KRT5", "KRT14", "KRT15", "DSP", # Keratinocytes
  "THY1", "COL1A1", "COL1A2", "COL3A1", "DCN", "MGP", "COL6A2", # Fibroblasts
  "CD3D", "CD8A", "CD4", "PTPRC", "FOXP3", "IKZF2", "CCL5", # T-cells
  "CD19", "MS4A1", # B-cells
  "CD14", "CD86", "CD74", "CD163", #Monocytes / macrophages
  "VWF", "PECAM1", "SELE", # Endothelial
  "MITF", "TYR", # Melanocyte markers
  "ITGAX", "CD1C", "CD1A", "CLEC1A", "CD207", # Dendritic cells (ITGAX = cd11c)
  "TPM1", "TPM2", "TAGLN" # Muscle
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.00", 
  labelMarkers = markerGenes,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE
)

hm <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(hm, name = "filtered-GeneScores-Marker-Heatmap", width = 6, height = 10, ArchRProj = proj, addDOC = FALSE)

nonMultipletCells <- getCellNames(proj)[proj$Clusters %ni% c("C7", "C13", "C15", "C18")]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = nonMultipletCells,
  outputDirectory = "multiplets_removed_output",
  dropCells=TRUE, force=TRUE
)
saveArchRProject(proj)

# Now, redo clustering and visualization:

# Reduce Dimensions with Iterative LSI
set.seed(1)
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI",
    sampleCellsPre = 20000,
    varFeatures = 50000, 
    dimsToUse = 1:50,
    force = TRUE
)

# Identify Clusters from Iterative LSI
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.7,
    force = TRUE
)

set.seed(1)
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 60, 
    minDist = 0.6, 
    metric = "cosine",
    force = TRUE
)

# Relabel clusters so they are sorted by cluster size
proj <- relabelClusters(proj)
proj <- addImputeWeights(proj)
# Make various cluster plots:
proj <- visualizeClustering(proj, pointSize=pointSize, sampleCmap=sample_cmap, diseaseCmap=disease_cmap)

# Save filtered ArchR project
saveArchRProject(proj)

########################################################################