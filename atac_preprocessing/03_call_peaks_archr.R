#!/usr/bin/env Rscript

#####################################################################
# Call peaks on high-resolution sub-clusters
#####################################################################

# ***NOTE***: Use condaMacs2 conda environment to ensure correct macs2 path

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
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# Working directory
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/baseline_preprocessing/multiplets_removed_output"

# New directory
outdir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
pointSize <- 0.25

##########################################################################################
# Load Previously Prepared ArchR project
##########################################################################################

orig_proj <- loadArchRProject(wd, force=TRUE)

##########################################################################################
# Copy ArchR project to new directory for adding subcluster information
##########################################################################################

atac_proj <- subsetArchRProject(
  ArchRProj = orig_proj,
  cells = getCellNames(orig_proj),
  outputDirectory = outdir,
  dropCells = TRUE,
  force = TRUE
)
saveArchRProject(atac_proj)
rm(orig_proj); gc()

# Delete stuff we don't need to duplicate
unlink(paste0(outdir, "/Plots/*"), recursive = TRUE)
unlink(paste0(outdir, "/GroupCoverages"), recursive = TRUE)
unlink(paste0(outdir, "/PeakCalls"), recursive = TRUE)
unlink(paste0(outdir, "/*_filtered_barcodes.txt"), recursive = TRUE)

##########################################################################################
# Identify cluster labels from subclustered ArchR projects
##########################################################################################

atac_proj <- loadArchRProject(outdir, force=TRUE)
# color palettes
sample_cmap <- readRDS("/home/users/boberrey/git_clones/hairATAC/sample_cmap.rds")
sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(atac_proj$Sample2)] %>% unlist()
atacNamedClustCmap <- readRDS(paste0(scriptPath, "/scATAC_NamedClust_cmap.rds")) %>% unlist()

subclustered_projects <- c("Lymphoid", "Myeloid", "Keratinocytes", "Fibroblasts", "Endothelial")

FineClustLabels <- atac_proj$NamedClust # Default to NamedClust where no subgroup label exists
names(FineClustLabels) <- getCellNames(atac_proj)

for(subgroup in subclustered_projects){
    message(sprintf("Reading in subcluster %s", subgroup))
    # Read in subclustered project
    sub_dir <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/subclustered_%s", subgroup)
    sub_proj <- loadArchRProject(sub_dir, force=TRUE)

    # Add FineClust to full ArchR project
    FineClustLabels[getCellNames(sub_proj)] <- sub_proj$FineClust
}

# Add FineCluster information to full ArchR project
atac_proj$FineClust <- FineClustLabels[getCellNames(atac_proj)]

# Plot the UMAPs by Sample and Cluster:
p1 <- plotEmbedding(ArchRProj=atac_proj, colorBy = "cellColData", name="NamedClust", pal=atacNamedClustCmap,
    embedding = "UMAP", plotAs="points", size=pointSize, labelMeans = FALSE)
p2 <- plotEmbedding(ArchRProj = atac_proj, colorBy="cellColData", name="FineClust",
    embedding="UMAP", plotAs="points", size=pointSize, labelMeans=FALSE)
ggAlignPlots(p1,p2, type="h")
plotPDF(p1,p2, name = "Plot-UMAP-FineClusters.pdf", ArchRProj=atac_proj, addDOC=FALSE, width=5, height=5)

saveArchRProject(atac_proj)

##########################################################################################
# Call Peaks
##########################################################################################

atac_proj <- loadArchRProject(outdir, force=TRUE)

# Create Group Coverage Files that can be used for downstream analysis
atac_proj <- addGroupCoverages(
  ArchRProj=atac_proj, 
  groupBy="FineClust", 
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
  )

# Find Path to Macs2 binary
# ***NOTE***: Use condaMacs2 conda environment to ensure correct macs2 path
pathToMacs2 <- findMacs2()

# Call Reproducible Peaks w/ Macs2
atac_proj <- addReproduciblePeakSet(
    ArchRProj = atac_proj, 
    groupBy = "FineClust", 
    peaksPerCell = 500, # The upper limit of the number of peaks that can be identified per cell-grouping in groupBy. (Default = 500)
    pathToMacs2 = pathToMacs2,
    force = TRUE
)

# Add Peak Matrix
atac_proj <- addPeakMatrix(ArchRProj = atac_proj, force = TRUE)

# Save project
saveArchRProject(atac_proj)

##########################################################################################





# # Add Motif Peak Annotations (Takes some time)
# atac_proj <- addMotifAnnotations(ArchRProj = atac_proj, motifSet = "cisbp", name = "Motif", force = TRUE)

# # Add background peaks
# atac_proj <- addBgdPeaks(atac_proj)

# # (WARNING: kind of slow, 20+ minutes)
# atac_proj <- addDeviationsMatrix(
#   ArchRProj = atac_proj, 
#   peakAnnotation = "Motif",
#   force = TRUE
# )


# # Plot the FRIP/TSS per sample/cluster:
# plotList <- list()
# plotList[[1]] <- plotGroups(
#   ArchRProj = atac_proj, 
#   groupBy = "Sample2", 
#   colorBy = "colData", 
#   name = "FRIP",
#   pal = sample_cmap
# )
# plotList[[2]] <- plotGroups(
#   ArchRProj = atac_proj, 
#   groupBy = "FineClust", 
#   colorBy = "colData", 
#   name = "FRIP",
#   ratioYX = 2
# )
# plotList[[3]] <- plotGroups(
#   ArchRProj = atac_proj, 
#   groupBy = "Sample2", 
#   colorBy = "colData", 
#   name = "TSSEnrichment",
#   pal = sample_cmap
# )
# plotList[[4]] <- plotGroups(
#   ArchRProj = atac_proj, 
#   groupBy = "FineClust", 
#   colorBy = "colData", 
#   name = "TSSEnrichment",
#   ratioYX = 2
# )
# plotPDF(plotList = plotList, name = "FRIP-TSS-Enrichment", width = 4, height = 4,  ArchRProj = atac_proj, addDOC = FALSE)


# #############################################################################