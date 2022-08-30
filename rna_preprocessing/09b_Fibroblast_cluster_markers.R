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

subgroup <- "Fibroblasts"
pointSize <- 1.5
useMagic <- TRUE # Should Rmagic be used for data imputation prior to UMAP plotting?

# change the current plan to access parallelization (for Seurat)
nThreads <- 8
plan("multiprocess", workers = nThreads)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory to Folder
wd <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scRNA_preprocessing/harmonized_subclustering/%s", subgroup)
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
fineClust <- sapply(1:nclust, function(x) paste0("rFb", x))
names(fineClust) <- 0:(nclust-1)

obj$FineClust <- fineClust[obj$Clusters] %>% unname
Idents(obj) <- "FineClust"

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

# Markers for identifying broad classes of cells:
hla_genes <- grep(pattern = "^HLA-", x = colnames(expr), value = TRUE)

featureSets <- list(
    "Fibroblasts" = c("THY1", "COL1A1", "COL1A2", "COL3A1", "DCN", "MGP", "COL6A2", "CEBPB", "APOD", "CFD"),
    "HF_associated" = c("APCDD1", "VCAN", "CORIN", "PTGDS", "SOX2", "COL11A1"), # Dermal Sheath? Hard to find clearly defined markers...
    "ImmuneRecruiting" = c("CXCL1", "CXCL2", "CXCL14", "CD44"),
    "AAdiff" = c("VEGFA", "CTGF", "FGF7", "OSMR"),
    "TGFBeta" = c("TGFB1", "TGFBR1", "TGFBR2", "TGFBR3", "SMAD1", "SMAD3", "SMAD4", "SMAD7", "ID1", "ID2", "ID3", "ID4"),
    "Papillary_dermis" = c("COL6A5", "APCDD1", "HSPB3", "WIF1", "ENTPD1"), # PMID: 29391249
    "Reticular_dermis" = c("CD36"),
    "Dermal_Papilla" = c("WNT5A", "BMP4", "BMP7", "HHIP", "PTCH1", "SOX18", "RUNX1", "RUNX3", "ALX4")
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
saveRDS(fineclust_cmap, file = sprintf("/home/users/boberrey/git_clones/hairATAC/rna_cmap_%s.rds", subgroup))

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


#####################################
# Subset of marker genes for figures
#####################################

source(paste0(scriptPath, "/cluster_labels.R"))

plot_genes <- c(
    "THY1", "COL1A1", "COL1A2", 
    "CXCL1", "CXCL2", # Fb_1
    "CCL19", "CXCL12",
    "APCDD1", "COL18A1", # rVe3
    "WISP2", "SCN7A", "CCL21", "PDGFD",
    "COL11A1", "EDNRA", "SOX2",
    "HHIP", "PTCH1", "WNT5A"
)

clustOrder <- c(
    "rFb1", # "Fb_1", # CXCL1,2,3
    "rFb3", # "Fb_2", # CCL19, CXCL12
    "rFb4", # "Fb_3", # APCDD1, COL18A1, F13A1
    "rFb5", # "Fb_4", # WISP2, AOX1, ARFGEF3
    "rFb6", # "Fb_5", # NECAB1, SCN7A
    "rFb2", # "D.Sheath", # COL11A1, EDNRA
    "rFb7" # "D.Papilla", # Many markers HHIP, PTCH1, etc.
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

#####################################
