#!/usr/bin/env Rscript

#####################################################################
# Build ArchR project and perform basic pre-processing and subsetting
#####################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
})

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

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
# Load Previously Prepared ArchR project
##########################################################################################

proj <- loadArchRProject(paste0(wd, "/multiplets_removed_output"), force=TRUE)
proj <- addImputeWeights(proj)

sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(proj$Sample2)] %>% unlist()

###############################
# Marker Genes on full project
###############################

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
  "KRT1", "KRT10", "KRT14", "KRT15", # Keratinocytes
  "ITGB8", "SOX9", "LGR5", "LHX2", "KRT16", "KRT75", # Hair follicle
  "THY1", "COL1A1", "COL11A1", # Fibroblasts
  "CD3D", "CD8A", "CD4", "FOXP3", "IKZF2", # T-cells
  "MS4A1", "IGLL5", # B-cells
  "CD14", "CD86", "CD74", "CCR7", "CD163", #Monocytes / macrophages
  "TPSB2", "FCER1A", "KIT", "HPGD", # Mast cells? (KIT and MITF also melanocytes)
  "VWF", "PECAM1", "SELE", # Endothelial
  "MITF", "TYR", "SOX10", # Melanocyte markers
  "ITGAX", "CD1C", "CD1A", "CD207", # Dendritic cells (ITGAX = cd11c)
  "TPM1", "TPM2", "MYL9", # Muscle
  "FOXE1", "SMIM23", "GJB6" # HF keratinocytes
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
plotPDF(hm, name = "GeneScores-Marker-Heatmap", width = 6, height = 10, ArchRProj = proj, addDOC = FALSE)

# Marker gene imputation with Magic:
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj), 
    plotAs="points", size = pointSize
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

# Tracks of genes:
p <- ArchRBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    tileSize=250,
    minCells=250
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

# Assign cluster names based on marker genes:

# Fb = Fibroblast
# Tc = T-cell / Lymphoid
# Ve = Vascular endothelial
# Le = Lymphoic endothelial
# Kc = Keratinocyte
# My = Myeloid
# Mu = Muscle / myofibroblast
# Me = Melanocytes
# Bc = B cells / Plasma
# Ma = Mast cells

clustNames <- list(
  "C1" = "aFb1",
  "C2" = "aTc1", 
  "C3" = "aKc1", 
  "C4" = "aMy1",
  "C5" = "aVe1", 
  "C6" = "aKc2",
  "C7" = "aMu1",
  "C8" = "aKc3",
  "C9" = "aFb2",
  "C10" = "aKc4",
  "C11" = "aKc5",
  "C12" = "aTc2",
  "C13" = "aMu2",
  "C14" = "aTc3",
  "C15" = "aMe1",
  "C16" = "aMy2",
  "C17" = "aKc6",
  "C18" = "aVe2",
  "C19" = "aLe1",
  "C20" = "aMy3",
  "C21" = "aBc1",
  "C22" = "aKc7"
  )
proj$NamedClust <- clustNames[proj$Clusters] %>% unlist()

# Assign 'broad clusters' as well
proj$BroadClust <- proj$NamedClust %>% gsub('[0-9]+', '', .) %>% sub('.', '', .)

# Determine colors to use for plots

# Reuse RNA colors for clusters with the same label
source(paste0(scriptPath, "/cluster_labels.R")) # Load labels from file
colorPal <- readRDS("/home/users/boberrey/git_clones/hairATAC/scalpClusterColors.rds") # Load pre-set color palate for all broad clusters
rnaNamedClustCmap <- readRDS(paste0(scriptPath, "/scRNA_NamedClust_cmap.rds")) %>% unlist()
rnaLabelClustCmap <- rnaNamedClustCmap
names(rnaLabelClustCmap) <- unlist(rna.NamedClust)[names(rnaLabelClustCmap)]

# atac cmap
# (First need to get an RNA to ATAC mapping)
labels_to_rna <- invertList(unlist(rna.NamedClust)[names(rnaNamedClustCmap)])
rna_to_atac <- invertList(atac.NamedClust)[names(labels_to_rna)]
names(rna_to_atac) <- unlist(labels_to_rna)
atac_to_rna <- invertList(rna_to_atac) # Removes non-matching clusters

atacNamedClustCmap <- rnaNamedClustCmap[unlist(atac_to_rna)]
names(atacNamedClustCmap) <- names(atac_to_rna)

# Get colors for clusters that do no have an equivalent RNA label
set.seed(1)
expandedColors <- getColorMap(cmaps_BOR$stallion, n=50)
expandedColors <- expandedColors[expandedColors %ni% c(colorPal, rnaNamedClustCmap)]

broadClust <- unique(atac_proj$BroadClust)
namedClust <- unique(atac_proj$NamedClust)
labelHierarchy <- lapply(broadClust, function(x) namedClust[grepl(x, namedClust)])
names(labelHierarchy) <- broadClust
labelMap <- invertList(labelHierarchy)
leftout <- namedClust[namedClust %ni% names(atacNamedClustCmap)]

leftoutColors <- c()
for(x in leftout){
    ccolor <- mostSimilarColors(colorPal[[labelMap[[x]]]], colorOptions=expandedColors, n=1)
    expandedColors <- expandedColors[expandedColors != ccolor]
    leftoutColors <- c(leftoutColors, ccolor)
}
names(leftoutColors) <- leftout

atacNamedClustCmap <- c(atacNamedClustCmap, leftoutColors)

saveRDS(atacNamedClustCmap, file="/home/users/boberrey/git_clones/hairATAC/scATAC_NamedClust_cmap.rds")

barwidth=0.9

# Plot the UMAPs by Sample and Cluster:
p1 <- plotEmbedding(ArchRProj=proj, colorBy="cellColData", name="Clusters", embedding="UMAP", plotAs="points", size=pointSize, labelMeans=FALSE)
p2 <- plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="NamedClust", pal=unlist(narrowColors),
    embedding = "UMAP", plotAs="points", size=pointSize, labelMeans = FALSE)
p3 <- plotEmbedding(ArchRProj = proj, colorBy="cellColData", name="BroadClust", pal=unlist(colorPal),
    embedding="UMAP", plotAs="points", size=pointSize, labelMeans=FALSE)
ggAlignPlots(p1,p2,p3, type="h")
plotPDF(p1,p2,p3, name = "Plot-UMAP-Named-Clusters.pdf", ArchRProj=proj, addDOC=FALSE, width=5, height=5)

### Stacked bar plot fraction sample in named and broad clusters ###
plotDir <- paste0(proj@projectMetadata$outputDirectory, "/Plots")

namedClustBySamp <- fractionXbyY(proj$NamedClust, proj$Sample2, add_total=TRUE, xname="NamedClust", yname="Sample")
pdf(paste0(plotDir, "/sampleByNamedClustBarPlot.pdf"))
print(stackedBarPlot(namedClustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

broadClustBySamp <- fractionXbyY(proj$BroadClust, proj$Sample2, add_total=TRUE, xname="BroadClust", yname="Sample")
pdf(paste0(plotDir, "/sampleByBroadClustBarPlot.pdf"))
print(stackedBarPlot(broadClustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

# Save labeled ArchR project
saveArchRProject(proj)

##########################################################################################
# Subcluster all major cell groups
##########################################################################################

# All major cell clusters will be subclustered for downstream analysis
# Peaks will also be called on subclustered groups

proj <- loadArchRProject(paste0(wd, "/multiplets_removed_output"), force=TRUE)

subClusterGroups <- list(
  "Lymphoid" = c("Tc"), 
  "Keratinocytes" = c("Kc"),
  "Myeloid" = c("Ma", "My"),
  "Fibroblasts" = c("Fb"),
  "Endothelial" = c("Ve", "Le"),
  "Other" = c("Me", "Bc", "Mu", "Other")
  )

# subgroups are now not non-overlapping
subClusterCells <- lapply(subClusterGroups, function(x){
  getCellNames(proj)[as.character(proj@cellColData[["BroadClust"]]) %in% x]
  })

subClusterArchR <- function(proj, subCells, outdir){
  # Subset an ArchR project for focused analysis

  message(sprintf("Subgroup has %s cells.", length(subCells)))
  sub_proj <- subsetArchRProject(
      ArchRProj = proj,
      cells = subCells,
      outputDirectory = outdir,
      dropCells = TRUE,
      force = TRUE
  )
  saveArchRProject(sub_proj)

  # Delete stuff we don't want to copy...
  unlink(paste0(outdir, "/Plots/*"), recursive = TRUE)
  unlink(paste0(outdir, "/IterativeLSI/*"), recursive = TRUE)
  unlink(paste0(outdir, "/Embeddings/*"), recursive = TRUE)
  unlink(paste0(outdir, "/ImputeWeights/*"), recursive = TRUE)
}

subgroups <- c("Lymphoid", "Myeloid", "Keratinocytes", "Fibroblasts", "Endothelial")

# Generate and cluster each of the subprojects
sub_proj_list <- lapply(subgroups, function(sg){
  message(sprintf("Subsetting %s...", sg))
  outdir <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/subclustered_%s", sg)
  subClusterArchR(proj, subCells=subClusterCells[[sg]], outdir=outdir)
})
names(sub_proj_list) <- subgroups

# Save project
saveArchRProject(proj)
