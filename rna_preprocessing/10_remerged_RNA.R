#!/usr/bin/env Rscript

######################################################
# Add subclustered information back to full project
######################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(future)
  library(Matrix)
})

#### Parameters ####

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/preprocessing_output"
setwd(wd)
plotDir <- paste0(wd,"/expression_plots")

# change the current plan to access parallelization (for Seurat)
nThreads <- 8
plan("multicore", workers = nThreads)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))

#######################################################################
# Read in previously created Seurat object
#######################################################################
message("Reading in data...")
obj <- readRDS(paste0(wd, '/scalp.rds'))

# color palettes
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)] %>% unlist()

pointSize <- 0.15

#######################################################################
# Add subclustered labels back to original Seurat project
#######################################################################

# Subgroups to retrieve data from
subgroups <- c("Lymphoid", "Myeloid", "Keratinocytes", "Fibroblasts", "Endothelial")

FineClustLabels <- obj$NamedClust
names(FineClustLabels) <- Cells(obj)

for(subgroup in subgroups){
    message(sprintf("Reading in subcluster %s", subgroup))

    # Read in subclustered object
    sub_dir <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/harmonized_subclustering/%s", subgroup)
    sub_obj <- readRDS(paste0(sub_dir, sprintf('/%s.rds', subgroup)))

    # Add ManualLabels to full Seurat object
    FineClustLabels[Cells(sub_obj)] <- sub_obj$FineClust
}

# Add FineCluster information to full Seurat object
obj$FineClust <- FineClustLabels[Cells(obj)]

### FineCluster UMAP ###
umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), obj$FineClust)
# Randomize cells before plotting
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir,"/FineClusters_UMAP.pdf"), width=10, height=10)
plotUMAP(umapDF, dataType="qualitative", cmap=cmaps_BOR$stallion, point_size=pointSize)
dev.off()

# Save whole project with all cluster information:
saveRDS(obj, file = paste0(wd, "/scalp.rds"))

#######################################################################