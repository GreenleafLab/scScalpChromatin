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

# set working directory
subgroup <- "Fibroblasts"
wd <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/subclustered_%s", subgroup)
full_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 1.0

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_proj <- readRDS(sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scRNA_preprocessing/harmonized_subclustering/%s/%s.rds", subgroup, subgroup))

plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

# Colormaps
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
atac_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(atac_proj$Sample2)] %>% unlist()
disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

rna_sub_cmap <- readRDS(sprintf("/home/users/boberrey/git_clones/hairATAC/rna_cmap_%s.rds", subgroup))
atac_sub_cmap <- readRDS(sprintf("/home/users/boberrey/git_clones/hairATAC/atac_cmap_%s.rds", subgroup))

###########################################################################################
# Do not proceed prior to calling peaks
###########################################################################################

# Compute group coverages
atac_proj <- addGroupCoverages(
  ArchRProj=atac_proj, 
  groupBy="FineClust", 
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
  )

# Get peaks that were called on this subproject's subclusters from full ArchR project
full_proj <- loadArchRProject(full_dir, force=TRUE)
full_peaks <- getPeakSet(full_proj)
peaks <- getClusterPeaks(full_proj, clusterNames=unique(atac_proj$FineClust), peakGR=full_peaks)
rm(full_proj); gc()

# Now add these peaks to the subproject and generate peak matrix
atac_proj <- addPeakSet(atac_proj, peakSet=peaks, force=TRUE)
atac_proj <- addPeakMatrix(atac_proj, force=TRUE)

atac_proj <- addMotifAnnotations(atac_proj, motifSet="cisbp", name="Motif", force=TRUE)

# Calculate coaccessibility
atac_proj <- addCoAccessibility(
    ArchRProj = atac_proj,
    reducedDims = "Harmony"
)

# Calculate peak-to-gene links
atac_proj <- addPeak2GeneLinks(
    ArchRProj = atac_proj,
    reducedDims = "Harmony"
)

# Add background peaks
atac_proj <- addBgdPeaks(atac_proj, force = TRUE)

# Add deviations matrix
atac_proj <- addDeviationsMatrix(
  ArchRProj = atac_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

# Save intermediate output
saveArchRProject(atac_proj)

##########################################################################################
