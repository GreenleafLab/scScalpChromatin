#!/usr/bin/env Rscript

#####################################
# Cluster scRNA using iterative LSI
#####################################

library(dplyr)
library(tidyr)
library(Seurat)
library(ggrastr)

pointSize <- 0.2

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin/"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/seurat_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/preprocessing_output"
plotDir <- paste0(wd,"/clustering_qc")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# color palettes
sample_cmap <- readRDS("/home/users/boberrey/git_clones/scScalpChromatin/sample_cmap.rds")

##########################################
# Read in previously created Seurat object
##########################################

obj <- readRDS(paste0(wd,'/scalp.rds'))

sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)]

disease_cmap <- head(cmaps_BOR$stallion, 3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB") 

# Plot clustering results:
message("Plotting clustering results...")
plotClusterQC(obj, subgroup="scalp", plotDir=plotDir, pointSize=pointSize, sampleCmap=sample_cmap, diseaseCmap=disease_cmap)

message("Done.")

