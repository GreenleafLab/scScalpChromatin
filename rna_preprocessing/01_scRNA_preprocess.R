#!/usr/bin/env Rscript

#####################################
# Cluster scRNA using iterative LSI
#####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(future) # For parallelization
  library(Matrix) # Required to work with sparse matrices
  library(stringr)
})

#### Parameters ####

# Misc
nThreads <- 8
plan("multicore", workers = nThreads)

# Cell quality filter criteria
minFeatures <- 200
maxFeatures <- Inf
minCounts <- 1000
maxCounts <- Inf
maxPctMito <- 20

# Seurat clustering parameters (used only for DoubletDetection and Ambient RNA decontamination)
seuratFeatures <- 2500
seuratDims <- 1:15
seuratRes <- 0.4 # Lower resolution will increase estimated RNA contamination

# DoubletFinder params
useDoubletFinder <- TRUE
estDubRate <- 0.04 # Rough estimate based on targeting ~8k cells

# DecontX params
useDecontX <- TRUE

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/preprocessing_output"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

logfile <- paste0(wd, sprintf("/preprocess_log_%s.txt", format(Sys.time(), "%Y%m%d-%H%M%S")))
con <- file(logfile, open = "wt")
sink(con, type="output")
sink(con, type="message")

# Print all parameters to log file
for ( obj in ls() ) { cat('---',obj,'---\n'); print(get(obj)) }

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin/"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/seurat_helpers.R"))
source(paste0(scriptPath, "/sample_metadata.R"))

plotDir <- paste0(wd,"/raw_qc")
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# Load each of the scalp datasets
message("Reading in data...")
#raw_data_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/scRNA_preprocessing/filtered_feature_bc_matrices/"
raw_data_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scRNA/filtered_feature_bc_matrices"
data_dirs <- list.dirs(path=raw_data_dir, full.names=TRUE, recursive=FALSE)
data_dirs <- data_dirs[grepl("/[A-Z_]+[0-9]+$", data_dirs)]
names(data_dirs) <- str_extract(data_dirs, "(?<=/)[^/]*$")

objs <- list()
for(ix in seq_along(data_dirs)){
  sample <- names(data_dirs)[ix]
  path <- data_dirs[ix]
  message(sprintf("Reading in data from sample %s...", sample))
  data <- Read10X(data.dir=path)
  obj <- CreateSeuratObject(counts=data, project=sample, min.cells=0, min.features=minFeatures) # orig.ident not getting assigned correctly?
  obj$orig.ident <- sample
  objs[[sample]] <- obj
}

objNames <- sapply(objs, function(x) x@project.name)

# Close connections
sink(type = "message")
sink(type = "output")
close(con)


# Perform initial scRNA preprocessing (cell calling, quality filtering, doublet detection, ambient RNA removal):
objs <- lapply(seq_along(objs), function(i){
  objName <- objNames[i]
  obj <- objs[[i]]
  scRNAdataPreProcessing(
    obj, objName, plotDir, # Naming parameters
    minFeatures=minFeatures, maxFeatures=maxFeatures, 
    minCounts=minCounts, maxCounts=maxCounts, maxPctMito=maxPctMito, # Quality filters
    nfeatures=seuratFeatures, dims=seuratDims, res=seuratRes, # Seurat clustering parameters
    runDoubletFinder=useDoubletFinder, runDecontX=useDecontX, # Optional processing steps
    estDubRate=estDubRate, # DoubletDetector parameters
    ncores=1
    )
  })

# Reopen main log
con <- file(logfile, open = "at")
sink(con, type="output", append=TRUE)
sink(con, type="message", append=TRUE)

# Merge individual Seurat objects
obj <- merge(x=objs[[1]], y=objs[2:length(objs)], add.cell.ids=objNames, project="scalp")

# Store sample information in metadata
obj$Sample <- obj$orig.ident

# Store disease status in metadata
obj$diseaseStatus <- NA
obj$diseaseStatus <- ifelse(grepl("C_SD", obj$Sample), "C_SD", obj$diseaseStatus)
obj$diseaseStatus <- ifelse(grepl("C_PB", obj$Sample), "C_PB", obj$diseaseStatus)
obj$diseaseStatus <- ifelse(grepl("AA", obj$Sample), "AA", obj$diseaseStatus)

# Add sample metadata
obj$preservation <- samp.preservation[obj$Sample] %>% unlist() %>% as.factor()
obj$sex <- samp.sex[obj$Sample] %>% unlist() %>% as.factor()
obj$age <- samp.age[obj$Sample] %>% unlist()

# Remove doublets
nPreDub <- dim(obj)[2]
obj <- subset(obj, subset = (DF.classify == "Singlet"))
nPostDub <- dim(obj)[2]
message(sprintf("Removed %s suspected doublets from %s total cells...", nPreDub - nPostDub, nPreDub))

# If contaminated RNA removed, some cells may have very little RNA now
if(useDecontX){
	nPreDecon <- dim(obj)[2]
	obj <- subset(obj, subset = (nFeature_RNA > minFeatures & nCount_RNA > minCounts))
	nPostDecon <- dim(obj)[2]
	message(sprintf("Removed %s cells with too little RNA after decontamination. %s cells remaining.", 
		nPreDecon - nPostDecon, nPostDecon))
}

# Save intermediate Seurat object:
message("Pre-processing complete. Saving intermediate Seurat object...")
saveRDS(obj, file = paste0(wd, "/preprocessed.rds"))

