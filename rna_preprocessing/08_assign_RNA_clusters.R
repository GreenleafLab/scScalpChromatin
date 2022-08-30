#!/usr/bin/env Rscript

#####################################################
# Subcluster scRNA to match subclustered scATAC data
#####################################################

library(dplyr)
library(tidyr)
library(Seurat)
library(ggrastr)
library(future) # For parallelization

# change the current plan to access parallelization (for Seurat)
nThreads <- 8
plan("multicore", workers = nThreads)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scRNA_preprocessing/preprocessing_output"
plotDir <- paste0(wd,"/expression_plots_scalp")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# color palettes
sample_cmap <- readRDS("/home/users/boberrey/git_clones/hairATAC/sample_cmap.rds")

##########################################
# Read in previously created Seurat object
##########################################
message("Reading in data...")
obj <- readRDS(paste0(wd, '/scalp.rds'))

sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)]

allGenes <- rownames(obj)

pointSize <- 0.15

##########################################
# Assigning cell type identity to clusters
##########################################

message("Identifying broad clusters and saving new objects for subclustering...")

# Fb = Fibroblast
# Tc = T-cell / Lymphoid
# Ve = Vascular endothelial
# Le = Lymphoic endothelial
# Kc = Keratinocyte
# My = Myeloid
# Mu = Muscle / myofibroblast
# Me = Melanocytes
# Bc = B cells
# Ma = Mast cells

clustNames <- list(
    "0" = "rFb1", 
    "1" = "rVe1",
    "2" = "rTc1",
    "3" = "rTc2",
    "4" = "rKc1",
    "5" = "rKc2",
    "6" = "rMu1",
    "7" = "rKc3",
    "8" = "rKc4",
    "9" = "rMy1", 
    "10" = "rTc3", 
    "11" = "rFb2", 
    "12" = "rMy2", 
    "13" = "rMu2", 
    "14" = "rKc5", 
    "15" = "rMa1", 
    "16" = "rMe1", 
    "17" = "rMy3", 
    "18" = "rLe1",
    "19" = "rMy4",
    "20" = "rBc1",
    "21" = "rMe2" 
)

obj$NamedClust <- clustNames[obj$Clusters] %>% unlist() %>% unname

# Assign 'broad clusters' as well
# (Remove digits and first character)
obj$BroadClust <- obj$NamedClust %>% gsub('[0-9]+', '', .) %>% sub('.', '', .)

# Plot the UMAPs by Sample and Cluster:

subClusterGroups <- list(
  "Lymphoid" = c("Tc"), 
  "Keratinocytes" = c("Kc"),
  "Myeloid" = c("My"),
  "Fibroblasts" = c("Fb"),
  "Endothelial" = c("Ve", "Le"),
  "Other" = c("Ma", "Me", "Bc", "Mu", "Other")
  ) %>% invertList()

obj$SubGroup <- subClusterGroups[obj$BroadClust] %>% unlist()

compartments <- list(
  "immune" = c("Tc", "My", "Ma", "Bc"), 
  "epidermal" = c("Kc", "Me"),
  "stromal" = c("Fb", "Mu", "Ve", "Le"),
  "other" = c("Other")
  ) %>% invertList()

obj$compartment <- compartments[obj$BroadClust] %>% unlist()

# Create color palette that will be used for all broad NamedClust
scalpClusterColors <- c(
  "Tc" = "#D51F26", # red
  "My" = "#8A9FD1", # sky blue 
  "Ma" = "#C06CAB", # light purple
  "Bc" = "#E6C2DC", # pale pink
  "Fb" = "#272E6A", # dark blue
  "Kc" = "#208A42", # green
  "Ve" = "#F47D2B", # orange
  "Le" = "#FEE500", # yellow
  "Mu" = "#89288F", # purple
  "Me" = "#D24B27", # brick
  "Other" = "#A9A9A9" # grey
  )

# Save these colors for use in other scripts
saveRDS(scalpClusterColors, file = "/home/users/boberrey/git_clones/hairATAC/scalpClusterColors.rds")
#scalpClusterColors <- readRDS("/home/users/boberrey/git_clones/hairATAC/scalpClusterColors.rds")
colorPal <- scalpClusterColors
expandedColors <- getColorMap(cmaps_BOR$stallion, n=50)
expandedColors <- expandedColors[expandedColors %ni% colorPal]

# Now identify colormap for subclusters
broadLabels <- unique(obj$BroadClust)
narrowLabels <- unique(obj$NamedClust)

labelHierarchy <- lapply(broadLabels, function(x) narrowLabels[grepl(x, narrowLabels)])
names(labelHierarchy) <- broadLabels
narrowColors <- list()
takenColors <- colorPal
for(bl in broadLabels){
    subLabs <- labelHierarchy[[bl]]
    if(length(subLabs) == 1){
        # If only a single label, use broad label
        narrowColors[[subLabs[1]]] <- colorPal[[bl]]
    }else{
        validColors <- expandedColors[expandedColors %ni% takenColors]
        subCols <- mostSimilarColors(colorPal[[bl]], colorOptions=validColors, n=length(subLabs))
        for(i in seq_along(subLabs)){
            narrowColors[[subLabs[i]]] <- subCols[i]
            takenColors <- c(takenColors, subCols[i])
        }
    }
}
# Save color palette for 'NamedClust'
saveRDS(narrowColors, file = "/home/users/boberrey/git_clones/hairATAC/scRNA_NamedClust_cmap.rds")

# Load colormaps for plotting:
barwidth=0.9
sample_cmap <- readRDS("/home/users/boberrey/git_clones/hairATAC/sample_cmap.rds")
sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)] %>% unlist()
scalpClusterColors <- readRDS(paste0(scriptPath, "/scalpClusterColors.rds")) %>% unlist()
narrowColors <- readRDS(paste0(scriptPath, "/scRNA_NamedClust_cmap.rds")) %>% unlist()

qualcmap <- cmaps_BOR$stallion

### Named cluster UMAP ###
umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), obj$NamedClust)
# Randomize cells before plotting
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir,"/NamedClusters_UMAP.pdf"), width=10, height=10)
plotUMAP(umapDF, dataType="qualitative", cmap=narrowColors, namedColors=TRUE, point_size=pointSize)
dev.off()

clustBySamp <- fractionXbyY(obj$NamedClust, obj$Sample, add_total=TRUE, xname="NamedClust", yname="Sample")
pdf(paste0(plotDir, "/clustBySampleBarPlot_NamedClusters.pdf"))
print(stackedBarPlot(clustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

clustByDisease <- fractionXbyY(obj$NamedClust, obj$diseaseStatus, add_total=TRUE, xname="NamedClust", yname="diseaseStatus")
disease_cmap <- head(cmaps_BOR$stallion, 3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB") 
pdf(paste0(plotDir, "/clustByDiseaseBarPlot_NamedClusters.pdf"))
print(stackedBarPlot(clustByDisease, cmap=disease_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()


### Broad cluster UMAP ###
umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$BroadClust)
# Randomize cells before plotting
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir,"/BroadClusters_UMAP.pdf"))
plotUMAP(umapDF, dataType="qualitative", cmap=scalpClusterColors, namedColors=TRUE, point_size=pointSize)
dev.off()

clustBySamp <- fractionXbyY(obj$BroadClust, obj$Sample, add_total=TRUE, xname="BroadClust", yname="Sample")
qualcmap <- cmaps_BOR$stallion
pdf(paste0(plotDir, "/clustBySampleBarPlot_BroadClusters.pdf"))
print(stackedBarPlot(clustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

clustByDisease <- fractionXbyY(obj$BroadClust, obj$diseaseStatus, add_total=TRUE, xname="BroadClust", yname="diseaseStatus")
qualcmap <- cmaps_BOR$stallion
pdf(paste0(plotDir, "/clustByDiseaseBarPlot_BroadClusters.pdf"))
print(stackedBarPlot(clustByDisease, cmap=disease_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

# Save whole project with all cluster information:
saveRDS(obj, file = paste0(wd, "/scalp.rds"))

########################################################################################################
# Subcluster groups of interest
########################################################################################################

# Make new Seurat objects for each sub-clustered group

makeSubClusts <- function(obj, ident, subgroups, outdir){
  Idents(obj) <- ident
  for(subg in subgroups){
    subsubdir <- paste0(outdir, sprintf("/%s", subg))
    dir.create(subsubdir, showWarnings = FALSE, recursive = TRUE)
    subObj <- subset(obj, idents = c(subg))
    counts <- GetAssayData(object = subObj, slot = "counts")
    newObj <- CreateSeuratObject(counts = counts, project = subg, min.cells = 0, min.features = 200)
    old.meta <- subObj@meta.data
    # Drop selected columns from old metadata
    old.cols <- colnames(old.meta)
    drop.cols <- old.cols[grepl("^RNA_snn", old.cols)]
    newObj@meta.data <- old.meta[,old.cols %ni% drop.cols]
    message(sprintf("Subcluster %s has %s cells", subg, dim(newObj)[2]))
    saveRDS(newObj, file = paste0(subsubdir, "/", subg, ".rds"))
  }
}

subclustDir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scRNA_preprocessing/harmonized_subclustering"
dir.create(subclustDir, showWarnings = FALSE, recursive = TRUE)

makeSubClusts(
  obj, 
  ident="SubGroup", 
  subgroups=c("Lymphoid", "Myeloid", "Keratinocytes", "Fibroblasts", "Endothelial"),
  outdir=subclustDir
)

