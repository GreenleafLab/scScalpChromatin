#!/usr/bin/env Rscript

#####################################
# Cluster scRNA using iterative LSI
#####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(future)
  library(Matrix)
  library(harmony)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

#### Parameters ####

# Misc
nThreads <- 8
plan("multicore", workers = nThreads)

# Iterative LSI params
nVarGenes <- 4000
nPCs <- 1:25
resolution <- c(0.1, 0.3, 0.6)

# Harmony
harmonize <- c() # Which iterations should be 'harmonized'
covariates <- c()

# UMAP:
umapNeighbors <- 50
umapMinDist <- 0.5
umapDistMetric <- "cosine"

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scRNA_preprocessing/preprocessing_output"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

logfile <- paste0(wd, sprintf("/reclustering_log_%s.txt", format(Sys.time(), "%Y%m%d-%H%M%S")))
con <- file(logfile, open = "wt")
sink(con, type="output")
sink(con, type="message")

# Print all parameters to log file
for ( obj in ls() ) { cat('---',obj,'---\n'); print(get(obj)) }

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin/"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/iterative_LSI.R"))
source(paste0(scriptPath, "/seurat_helpers.R"))

#######################################
# Perform clustering with Iterative LSI
#######################################

message("Reading in previously created Seurat object...")
obj <- readRDS(paste0(wd, "/preclustered.rds"))

# Remove spurious clusters from initial clustering
invalid.clusters <- as.character(c(18, 26, 28, 29))
obj <- subset(obj, idents = c(invalid.clusters), invert=TRUE)

rawCounts <- GetAssayData(object = obj, slot = "counts")

# Identify genes we want to blacklist during clustering

# mitochondrial:
mt.genes <- grep(pattern = "^MT-", x = rownames(rawCounts), value = TRUE)
# Cell cycle (These are loaded by Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# X/Y chromosome genes:
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
geneGR <- GenomicFeatures::genes(txdb)
sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
matchedGeneSymbols <- select(org.Hs.eg.db,
                        keys = sexGenesGR$gene_id,
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "ENTREZID")
sexChr.genes <- matchedGeneSymbols$SYMBOL


# Genes to ignore (just for clustering purposes)
blacklist.genes <- c(
    mt.genes,
    sexChr.genes,
    s.genes,
    g2m.genes
)

#Initialize list for storing iterative LSI output
lsiOut <- list()
clusters <- NULL

# Depth normalize to 10,000, add pseudo count, and then log2 transform
log2CP10k <- sparseLogX(rawCounts, logtype="log2", scale=TRUE, scaleFactor=10^4)
# Store the log2CP10k
obj <- SetAssayData(object = obj, slot = "data", new.data = log2CP10k)

message("Running iterative LSI...")
set.seed(1)
for(i in seq_along(resolution)){
    # If first round, compute variable genes on raw data first
    if(i == 1){
        message(sprintf("Identifying top %s variable genes among all cells...", nVarGenes))
        varGenes <- getVarGenes(log2CP10k, nvar = nVarGenes, blacklist = blacklist.genes)
    }else{
        # For remaining rounds, calculate variable genes using previous clusters
        clusterMat <- edgeR::cpm(groupSums(rawCounts, clusters, sparse = TRUE), log=TRUE, prior.count = 3)
        message(sprintf("Identifying top %s variable genes from round %s LSI...", nVarGenes, i-1))
        varGenes <- getVarGenes(clusterMat, nvar = nVarGenes, blacklist = blacklist.genes)
    }
    # Now run LSI and find clusters
    LSIi <- runLSI(rawCounts[varGenes,], nComponents = max(nPCs), binarize = FALSE)

    # 'Harmonize' SVD PCs, if indicated
    if(i %in% harmonize){
      message(sprintf("Harmonizing LSI SVD PCs for round %s...", i))
      harmonized_pcs <- HarmonyMatrix(
        data_mat  = LSIi$matSVD,
        meta_data = obj@meta.data,
        vars_use  = covariates, # Covariates to 'harmonize'
        do_pca    = FALSE
        )
      LSIi$matSVD <- harmonized_pcs
    }

    reducName <- paste0("LSI_iter",i)
    obj[[reducName]] <- CreateDimReducObject(embeddings = LSIi$matSVD, key = sprintf("LSI%s_", i), assay = "RNA")
    obj <- FindNeighbors(object = obj, reduction = reducName, dims = nPCs, force.recalc = TRUE)
    message(sprintf("Clustering with resolution %s...", resolution[i]))
    obj <- FindClusters(object = obj, resolution = resolution[i])
    clusters <- Idents(obj)
    #Store information
    lsiOut[[reducName]] <- list(
        lsiMat = LSIi$matSVD, 
        varFeatures = varGenes, 
        clusters = clusters
    )
}

# Store cluster information in metadata
obj$Clusters <- Idents(obj)

##################################################
# Run non-linear dimensional reduction (UMAP/tSNE)
##################################################

# Seurat uses the uwot implementation of UMAP by default
message("Calculating UMAP...")
set.seed(1)
obj <- RunUMAP(
  obj,
  reduction = paste0("LSI_iter",length(resolution)), # Use final LSI iteration 
  dims = nPCs,
  n.neighbors = umapNeighbors,
  min.dist = umapMinDist,
  metric = umapDistMetric
)

message("Saving seurat object...")

# Save clustered object here:
saveRDS(obj, file = paste0(wd, "/scalp.rds"))

# Save iterativeLSI info
saveRDS(lsiOut, file = paste0(wd, "/lsiOut_scalp.rds"))

message("Done.")

# Close connections
on.exit({ sink(type = "message"); sink(type = "output"); close(con) })


