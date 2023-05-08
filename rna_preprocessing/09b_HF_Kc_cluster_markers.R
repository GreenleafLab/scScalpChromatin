#!/usr/bin/env Rscript

############################################
# Identify marker genes and make some plots
############################################

library(dplyr)
library(tidyr)
library(Seurat)
library(ggrastr)
library(Rmagic)
library(future)
library(data.table)

subgroup <- "HF_Kc"
pointSize <- 2.0
useMagic <- TRUE # Should Rmagic be used for data imputation prior to UMAP plotting?

# change the current plan to access parallelization (for Seurat)
nThreads <- 8
plan("multicore", workers = nThreads)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory to Folder
wd <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/harmonized_subclustering/%s", subgroup)
plotDir <- paste0(wd,"/expression_plots")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

##########################################
# Read in previously created Seurat object
##########################################

obj <- readRDS(paste0(wd, sprintf('/%s.rds', subgroup)))
allGenes <- rownames(obj)

# Now, assign cluster names:
nclust <- length(unique(obj$Clusters))
fineClust <- sapply(1:nclust, function(x) paste0("rHF", x))
names(fineClust) <- 0:(nclust-1)

obj$HFClust <- fineClust[obj$Clusters] %>% unname
Idents(obj) <- "HFClust"

##########################################
# Identify markers per cluster (And GO terms)
##########################################

# find markers for every cluster compared to all remaining cells, report only the positive ones
message("Finding marker genes using Seurat...")
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# Save markers
write.table(obj.markers, file=paste0(wd, sprintf("/marker_genes_%s.tsv", subgroup)), 
    sep="\t", quote=FALSE,row.names=FALSE)

# Use only expressed genes in GO term analyses:
# Define 'expressed genes' as those with at least 2 counts in at least 5 cells
rawCounts <- GetAssayData(object = obj, slot = "counts")
minUMIs <- 2
minCells <- 3
expressedGenes <- rownames(rawCounts[rowSums(rawCounts > minUMIs) > minCells,])

message("Calculating GO terms on each cluster using marker genes...")

clusters <- unique(obj$HFClust)

# Get GO enrichments (expressed genes):
GOresults <- lapply(clusters, function(x){
  message(sprintf("Running GO enrichments on cluster %s...", x))
  df <- obj.markers[obj.markers$cluster == x,]
  # Define custom cutoff for 'interesting genes'
  intGenes <- df[(df$p_val_adj < 0.01) & (df$avg_log2FC > 0.5),"gene"]
  calcTopGo(expressedGenes, interestingGenes=intGenes)
  })
names(GOresults) <- paste0("cluster", clusters)

# Plots of GO term enrichments:
pdf(paste0(plotDir, "/cluster_marker_GO_term_enrichments.pdf"), width=8, height=8)
for(name in names(GOresults)){
    goRes <- GOresults[[name]]
    print(topGObarPlot(goRes, cmap = cmaps_BOR$comet, 
        nterms=10, border_color="black", 
        barwidth=0.85, title=name))
}
dev.off()

# Save a file of the results for each cluster:
goDir <- paste0(wd,"/go_term_enrichments")
dir.create(goDir, showWarnings = FALSE, recursive = TRUE)
for(n in names(GOresults)){
    write.table(GOresults[[n]], file=paste0(goDir, "/", n, "_go_terms.tsv"), quote=FALSE, sep='\t')
}

###########################
# UMAP of high level groups
###########################
message("Plotting selected marker features on UMAP...")

# Set colormaps
qualcmap <- cmaps_BOR$stallion
quantcmap <- cmaps_BOR$sunrise

# Get expression data:
expr <- GetAssayData(obj, slot = 'data') %>% t()
expr <- expr[,Matrix::colSums(expr) > 0] # Remove unexpressed genes

featureSets <- list(
    # https://www.proteinatlas.org/humanproteome/tissue/skin
    "Basal" = c("KRT14", "KRT5", "KRT15", "COL17A1"), # Basal epithelia
    "HFSCs" = c("ITGA6", "ITGB1", "CD200", "LGR5","LHX2", "FRZB", "FZD1", "FZD5", "FZD10",  "IL31RA", "OSMR"), # HFSCs
    "HairGerm" = c("CD34", "CDH3", "LGR5", "CDKN2A", "RUNX1"), # Hair germ markers
    "Matrix" = c("KRT81", "KRT83", "HOXC13", "LEF1"), # Matrix hair keratins/genes
    "Sheath" = c("KRT71", "KRT75"), # IRS/ORS keratins / genes
    "TFs" = c("SOX9", "LHX2", "NFATC1", "TCF3") # Key HFSC TFs
)

selectedGenes <- unlist(featureSets) %>% unname()

flag <- "noMagic"
# Smooth with Rmagic
if(useMagic){
    message("Using MAGIC to impute (smooth) data for plotting...")

    # Run MAGIC directly on the expression matrix
    expr <- magic(expr, genes=selectedGenes, n.jobs = 1, seed = 1)$result
    flag <- "yesMagic"
}

for(name in names(featureSets)){
    features <- featureSets[[name]]
    pdf(paste0(plotDir,"/", name, "_features_UMAP.pdf"))
    for(gene in features){
        if(!gene %in% allGenes){
          message(sprintf("Error: %s is not a valid gene name", gene))
        }else if(!gene %in% colnames(expr)){
          message(sprintf("Error: %s is not expressed in any of these cells", gene))
        }else{
            umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), expr[,gene])        
            colnames(umapDF) <- c("UMAP1", "UMAP2", gene)
            # Clip range of expression:
            upperLim <- quantile(umapDF$gene, probs=c(0.95))
            umapDF[,gene][umapDF[,gene] >= upperLim] <- upperLim
            print(plotUMAP(umapDF, dataType = "quantitative", cmap = quantcmap, covarLabel = gene, point_size = pointSize))
        } 
    }
    dev.off()
}


# Dot plot of marker Genes:
count_mat <- GetAssayData(object = obj, slot = "counts")
avgPctMat <- avgAndPctExpressed(count_mat, obj$Clusters, feature_normalize=TRUE, min_pct=5)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")
wide_df <- prettyOrderMat(wide_df)

grp_order <- colnames(wide_df$mat)
gene_order <- rev(rownames(wide_df$mat))

pdf(paste0(plotDir, "/markers_dot_plot.pdf"), width=6, height=10)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", xorder=grp_order, yorder=gene_order, cmap=cmaps_BOR$sunrise)
dev.off()

# Order labels by frequency:
fineclust_cmap <- cmaps_BOR$stallion[1:length(fineClust)]
names(fineclust_cmap) <- names(getFreqs(obj$HFClust))
# Save color palette for 'NamedClust'
saveRDS(fineclust_cmap, file = paste0(scriptPath, sprintf("/rna_cmap_%s.rds", subgroup)))

### Cluster UMAP ###
umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), obj$HFClust)
# Randomize cells before plotting UMAP
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir, "/HFClust_UMAP.pdf"), width=10, height=10)
print(plotUMAP(umapDF, dataType = "qualitative", cmap = fineclust_cmap, namedColors=TRUE, point_size=pointSize))
dev.off()

# Save object with manual cluster labels
saveRDS(obj, file = paste0(wd, sprintf('/%s.rds', subgroup)))

#################################################################################
