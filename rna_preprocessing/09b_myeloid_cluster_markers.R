#!/usr/bin/env Rscript

############################################################################################################
# Identify marker genes and make some plots
############################################################################################################

library(dplyr)
library(tidyr)
library(Seurat)
library(ggrastr)
library(Rmagic)
library(future)
library(data.table)

subgroup <- "Myeloid"
pointSize <- 1.0
useMagic <- TRUE # Should Rmagic be used for data imputation prior to UMAP plotting?

# change the current plan to access parallelization (for Seurat)
nThreads <- 8
plan("multicore", workers = nThreads)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory to Folder
wd <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/harmonized_subclustering/%s", subgroup)
plotDir <- paste0(wd,"/expression_plots")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

############################################################################################################
# Read in previously created Seurat object
############################################################################################################

obj <- readRDS(paste0(wd, sprintf('/%s.rds', subgroup)))
allGenes <- rownames(obj)

# Now, assign cluster names:
nclust <- length(unique(obj$Clusters))
fineClust <- sapply(1:nclust, function(x) paste0("rMy", x))
names(fineClust) <- 0:(nclust-1)

obj$FineClust <- fineClust[obj$Clusters] %>% unname
Idents(obj) <- "FineClust"

############################################################################################################
# Identify markers per cluster (And GO terms)
############################################################################################################

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
minCells <- 5
expressedGenes <- rownames(rawCounts[rowSums(rawCounts > minUMIs) > minCells,])

message("Calculating GO terms on each cluster using marker genes...")

clusters <- unique(obj.markers$cluster)
Idents(obj) <- "Clusters"

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

############################################################################################################
# UMAP of high level groups
############################################################################################################
message("Plotting selected marker features on UMAP...")

# Set colormaps
qualcmap <- cmaps_BOR$stallion
quantcmap <- cmaps_BOR$sunrise

# Markers for identifying broad classes of cells:
featureSets <- list(
    # APC subtypes:
    #"Mast_cells" = c("KIT", "ENPP3", "FCER1A", "IL1RL1", "TPSB2"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase
    "Macrophages" = c("CD163", "LGMN", "FCGR2A", "C1QB", "C5AR1", "MAFB", "FOLR2"),
    "M1_macs" = c("CCL20", "CXCL3", "IL1B", "IL6", "IL12A", "IFNG", "TNF", "CD163"), # CD163 should be NEGATIVE
    "M2a_macs" = c("CD163", "CD200", "IRF4", "TGFB1", "TGFB2", "CCL2", "STAT6"),
    "TREM2_macs" = c("TREM2", "C3", "FCGBP", "FCGR3A", "OSM", "APOE"),
    "Langerhans_cells" = c("ITGAX", "CD1A", "CLEC1A", "CD207", "EPCAM"), # CD11c = ITGAX, Langerin = CD207, CLEC4K
    "pDC" = c("CCR7", "PTPRC", "CD209", "CLEC4C"), # PTPRC = CD45RA; IFNA1 = interferon alpha
    "moDC" = c("CD14", "CD1A", "CD1C", "ITGAX", "ITGAM", "SIRPA"), # SIRPA = CD172a
    "cDC1" = c("BTLA", "ITGAE", "CD1A", "ITGAM", "CLEC9A", "XCR1", "THBD"), # THBD = CD141 = BDCA-3 (thrombomodulin)
    "cDC2" = c("CD14", "CD163", "CLEC10A", "NOTCH2", "ITGAM", "SIRPA", "CX3CR1", "CD1C", "CD2"), # THBD = CD141 = BDCA-3 (thrombomodulin)
    "TCR_macs" = c("CD3D", "TRAC", "TRBC1", "SPOCK2", "CD14", "CD2")
)

selectedGenes <- unlist(featureSets) %>% unname()

# Get expression data:
expr <- GetAssayData(obj, slot = 'data') %>% t()
expr <- expr[,Matrix::colSums(expr) > 0] # Remove unexpressed genes

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
    pdf(paste0(plotDir,"/", name, "_features_UMAP_", subgroup, ".pdf"))
    for(gene in features){
        if(!gene %in% allGenes){
            message(sprintf("Error: %s is not a valid gene name", gene))
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
avgPctMat <- avgAndPctExpressed(count_mat, obj$FineClust, feature_normalize=TRUE, min_pct=5)

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
names(fineclust_cmap) <- names(getFreqs(obj$FineClust))
# Save color palette for 'NamedClust'
saveRDS(fineclust_cmap, file = paste0(scriptPath, sprintf("/rna_cmap_%s.rds", subgroup)))

### Cluster UMAP ###
umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), obj$FineClust)
# Randomize cells before plotting UMAP
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir, "/FineClust_UMAP.pdf"), width=10, height=10)
print(plotUMAP(umapDF, dataType = "qualitative", cmap = fineclust_cmap, namedColors=TRUE, point_size=pointSize))
dev.off()

### Labeled UMAP ###
source(paste0(scriptPath, "/cluster_labels.R"))
label_cmap <- fineclust_cmap
names(label_cmap) <- unlist(rna.FineClust)[names(label_cmap)]
umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), unlist(rna.FineClust)[obj$FineClust])
# Randomize cells before plotting UMAP
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir, "/LabeledClust_UMAP.pdf"), width=10, height=10)
print(plotUMAP(umapDF, dataType = "qualitative", cmap = label_cmap, namedColors=TRUE, point_size=pointSize))
dev.off()

# Save object with manual cluster labels
saveRDS(obj, file = paste0(wd, sprintf('/%s.rds', subgroup)))


############################################################################################################
# Subset of marker genes for figures
############################################################################################################

source(paste0(scriptPath, "/cluster_labels.R"))

plot_genes <- c(
    "IL15", "CCR7", "CCL19", "CCL17",
    "CD3D", "TRAC",
    "FOLR2", "C1QA", "CD163",
    "CXCL2", "CCL20",
    "TREM2", "OSM", "CD14",
    "FCER1A", "CLEC10A", "CD1C",
    "CLEC9A", "XCR1", 
    "IFNG", "CD207",
    "JCHAIN"
)

clustOrder <- c(
    "rMy5", # "M1.macs", # IL15, IL32, CCR7 (CCL19?, CCL17?)
    "rMy6", # "TCR.macs", # CD3, TCR gene positive macrophages?
    "rMy2", # "M2.macs_1", # C1Qa/b/c, FOLR2, CD14, CD163, (CCL13?)
    "rMy3", # "M2.macs_2", # CXCL2, CXCL3, (CCL20?, S100A8/9?) 
    "rMy7", # "TREM2.macs", # TREM2
    "rMy1", # "cDC2", # CD1c, CLEC10a (conventional DCs - type 2)
    "rMy4", # "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1
    "rMy8" # Plasma cell contamination / doublets
)

# Dot plot of marker Genes:
count_mat <- GetAssayData(object=obj, slot="counts")
avgPctMat <- avgAndPctExpressed(count_mat, obj$FineClust, feature_normalize=TRUE, min_pct=5)

# Subset to genes we care about:
avgPctMat <- avgPctMat[avgPctMat$feature %in% plot_genes,]

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")
wide_df <- prettyOrderMat(wide_df[,clustOrder], clusterCols=FALSE)

# Assign labels
avgPctMat$grp <- unlist(rna.FineClust)[as.character(avgPctMat$grp)]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0

grp_order <- colnames(wide_df$mat)
gene_order <- rev(rownames(wide_df$mat))

pdf(paste0(plotDir, sprintf("/markers_dot_plot_%s.pdf", subgroup)), width=4.5, height=5.5)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", 
    xorder=unlist(rna.FineClust)[grp_order], yorder=rev(gene_order), cmap=cmaps_BOR$sunrise, aspectRatio=1.6)
dev.off()

############################################################################################################

