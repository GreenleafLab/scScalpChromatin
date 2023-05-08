#!/usr/bin/env Rscript

#####################################
# Cluster scRNA using iterative LSI
#####################################

# Subcluster previously identified groups

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

# Set/Create Working Directory to Folder:
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/harmonized_subclustering"
setwd(wd)

# Subgroups to subcluster:
subgroups <- c("Lymphoid", "Myeloid", "Keratinocytes",  "Fibroblasts", "Endothelial")
#subgroups <- c("HF_Kc")

# Subclustering parameters:
paramDict <- list(
  "Lymphoid" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.1, 0.3),
    "nNeighbors" = 35,
    "minDist" = 0.4,
    "pointSize" = 0.75,
    "harmonize" = c(1,2), # Which iterations should be 'harmonized'
    "covariates" = c("Sample")
    ),
  "Keratinocytes" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.2, 0.3),
    "nNeighbors" = 35,
    "minDist" = 0.4,
    "pointSize" = 0.75,
    "harmonize" = c(1,2),
    "covariates" = c("Sample")
    ),
  "Myeloid" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.1, 0.3),
    "nNeighbors" = 35,
    "minDist" = 0.4,
    "pointSize" = 1.0,
    "harmonize" = c(1,2),
    "covariates" = c("Sample")
    ), 
  "Fibroblasts" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.2, 0.3),
    "nNeighbors" = 35,
    "minDist" = 0.4,
    "pointSize" = 1.00,
    "harmonize" = c(1,2),
    "covariates" = c("Sample")
    ),
  "Endothelial" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.1, 0.3),
    "nNeighbors" = 40,
    "minDist" = 0.35,
    "pointSize" = 1.00,
    "harmonize" = c(1,2),
    "covariates" = c("Sample")
    ),
  "HF_Kc" = list(
    "nVarGenes" = 1500,
    "nPCs" = 1:20,
    "lsiRes" = c(0.2, 0.4),
    "nNeighbors" = 20,
    "minDist" = 0.1,
    "pointSize" = 1.50,
    "harmonize" = c(),
    "covariates" = c() 
    )
  )

# Misc
nThreads <- 8
umapDistMetric <- "cosine"

# change the current plan to access parallelization (for Seurat)
plan("multicore", workers = nThreads)

# Start logging:
logfile <- paste0(wd, sprintf("/subclustering_log_%s.txt", format(Sys.time(), "%Y%m%d-%H%M%S")))
con <- file(logfile, open = "wt")
sink(con, type="output")
sink(con, type="message")

# Print all parameters to log file
for ( obj in ls() ) { cat('---',obj,'---\n'); print(get(obj)) }

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/iterative_LSI.R"))
source(paste0(scriptPath, "/seurat_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# color palettes
sample_cmap <- readRDS("/home/users/boberrey/git_clones/scScalpChromatin/sample_cmap.rds")
disease_cmap <- head(cmaps_BOR$stallion, 3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB") 

# Identify genes we want to blacklist during clustering

# First get list of all genes:
subwd <- paste0(wd, sprintf("/%s", subgroups[1]))
allGenes <- rownames(GetAssayData(object=readRDS(paste0(subwd, sprintf('/%s.rds', subgroups[1]))), slot="counts"))

# Identify genes we want to blacklist during clustering

# mitochondrial:
mt.genes <- grep(pattern="^MT-", x=allGenes, value=TRUE)
# Cell cycle: (These are loaded by Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# Ribosomal:
rp.genes <- grep(pattern="^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", x=allGenes, value=TRUE)

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
    g2m.genes,
    rp.genes
)


# Iteratively cluster subgroups:
for(sg in subgroups){

  subwd <- paste0(wd, sprintf("/%s", sg))

  # Directory for clustering qc plots:
  plotDir <- paste0(subwd, "/clustering_qc")
  dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

  # Read in previously created Seurat subobjects:
  message(sprintf("Reading in data for subgroup %s...", sg))
  obj <- readRDS(paste0(subwd, sprintf('/%s.rds', sg)))

  # Remove any existing DimReductions
  obj <- DietSeurat(obj, features=NULL, assays=NULL, dimreducs=NULL)

  # Subset colors to only those samples present
  samp_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)] %>% unlist()

  #######################################
  # Perform clustering with Iterative LSI
  #######################################

  # Get subgroup parameters:

  # Iterative LSI
  nVarGenes <- paramDict[[sg]]$nVarGenes
  nPCs <- paramDict[[sg]]$nPCs
  resolution <- paramDict[[sg]]$lsiRes

  # UMAP:
  umapNeighbors <- paramDict[[sg]]$nNeighbors
  umapMinDist <- paramDict[[sg]]$minDist

  # Harmony:
  harmonize <- paramDict[[sg]]$harmonize
  covariates <- paramDict[[sg]]$covariates

  rawCounts <- GetAssayData(object = obj, slot = "counts")

  #Initialize list for storing iterative LSI output
  lsiOut <- list()
  clusters <- NULL

  # Depth normalize to 10,000, add pseudo count, and then log2 transform
  log2CP10k <- sparseLogX(rawCounts, logtype="log2", scale=TRUE, scaleFactor=10^4)
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
          svd = LSIi$svd,
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
  saveRDS(obj, file = paste0(subwd, sprintf("/%s.rds", sg)))

  # Save iterativeLSI info
  saveRDS(lsiOut, file = paste0(subwd, sprintf("/lsiOut_%s.rds", sg)))

  # Plot clustering results:
  message("Plotting clustering results...")
  pointSize <- paramDict[[sg]]$pointSize
  plotClusterQC(obj, subgroup=sg, plotDir=plotDir, pointSize=pointSize, sampleCmap=sample_cmap, diseaseCmap=disease_cmap)

  message("Done.")

}



