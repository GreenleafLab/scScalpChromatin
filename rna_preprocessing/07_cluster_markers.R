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

subgroup <- "scalp"
pointSize <- 0.5
useMagic <- TRUE # Should Rmagic be used for data imputation prior to UMAP plotting?

# change the current plan to access parallelization (for Seurat)
nThreads <- 8
plan("multicore", workers = nThreads)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin/"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/preprocessing_output"
plotDir <- paste0(wd,"/expression_plots_scalp")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

##########################################
# Read in previously created Seurat object
##########################################
message("Reading in data...")
obj <- readRDS(paste0(wd, sprintf('/%s.rds', subgroup)))

allGenes <- rownames(obj)

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

rm(rawCounts); gc()

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

# Markers for identifying broad classes of cells:
featureSets <- list(
    "Basal_epithelia" = c("KRT15", "KRT5", "COL17A1"),
    "Spinous" = c("KRT1"),
    "HF_keratinocytes" = c("KRT75", "SOX9", "LHX2","ITGB8", "KRT16", "KRT17", "RUNX3"),
    "Glandular" = c("KRT7"),
    "T_cells" = c("CD3D", "CD8A", "CD4", "FOXP3", "IKZF2", "IFNG"),
    "B_cells" = c("MS4A1"), # MS41A = CD20, MZB1 = marginal zone B cell specific protein
    "M1_macs" = c("CCL20", "CD80", "CD86"),
    "M2a_macs" = c("CD163", "TGFB2"),
    "TREM2_macs" = c("TREM2", "OSM"),
    "FOLR2_macs" = c("FOLR2"),
    "CD1a1c_DCs" = c("CD1A", "CD1C", "ITGAX", "ITGAM"), # SIRPA = CD172a
    "CD1a141_DCs" = c("CLEC9A", "XCR1"), # THBD = CD141 = BDCA-3 (thrombomodulin)
    "Mast_cells" = c("KIT", "TPSB2"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase, which are supposed to be common in mast cells
    "Melanocytes" = c("MITF", "SOX10", "MLANA"), # Melanocyte markers
    "Endothlial" = c("VWF", "PECAM1", "SELE"), # Endothlial cells (PECAM1 = CD31), SELE = selectin E (found in cytokine stimulated endothelial cells)
    "Lymphatic" = c("FLT4", "LYVE1", "CCL21"),  # Lymphatic endothelial (FLT4 = VEGFR-3))
    "Angiogenic" = c("SEMA3G"),
    "Muscle" = c("TPM1", "TAGLN"), # TPM = tropomyosin, TAGLN = transgelin (involved crosslinking actin in smooth muscle), MYL = Myosin light chain
    "Fibroblasts" = c("THY1", "COL1A1"),
    "Dermal_sheath" = c("SOX2", "COL11A1"), # Dermal Sheath? Hard to find clearly defined markers...
    "Papillary_dermis" = c("COL6A5", "APCDD1"), # PMID: 29391249
    "Reticular_dermis" = c("CD36"), # CD36 seems more specific for the 'muscle 2' cluster... Myofibroblast?
    "Dermal_Papilla" = c("BMP7", "HHIP", "PTCH1", "SOX18"),
    "cycling" = c("MKI67", "CDK1", "TOP2A")
)

# Get expression data:
expr <- GetAssayData(obj, slot = 'data') %>% t()
expr <- as(expr[,Matrix::colSums(expr) > 0], "sparseMatrix") # Remove unexpressed genes

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
        }else{
            umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), expr[,gene])        
            colnames(umapDF) <- c("UMAP1", "UMAP2", gene)
            # Clip range of expression:
            upperLim <- quantile(umapDF$gene, probs=c(0.95))
            umapDF[,gene][umapDF[,gene] >= upperLim] <- upperLim
            print(plotUMAP(umapDF, dataType = "quantitative", cmap = quantcmap, covarLabel = gene, point_size=pointSize))
        } 
    }
    dev.off()
}

# Markers for identifying broad classes of cells:
featureSets <- list(
    "Keratinocytes" = c("KRT5", "KRT10", "KRT14", "KRT15"), 
    "Fibroblasts" = c("THY1", "COL1A1", "COL11A1"), 
    "T_cells" = c("CD3D", "CD8A", "CD4","IKZF2", "CCL5"), # PTPRC = CD45 
    "B_cells" = c("CD19"), # MS41A = CD20, MZB1 = marginal zone B cell specific protein
    "APCs" = c("CD14", "CD86", "CD74", "CD163"), # Monocyte lineage (FCGR3A = CD16, FCGR1A = CD64, CD74 = HLA-DR antigens-associated invariant chain)
    "Melanocytes" = c("MITF", "SOX10", "MLANA"), # Melanocyte markers
    "Endothlial" = c("VWF", "PECAM1", "SELE"), # Endothlial cells (PECAM1 = CD31), SELE = selectin E (found in cytokine stimulated endothelial cells)
    "Lymphatic" = c("FLT4", "LYVE1"),  # Lymphatic endothelial (FLT4 = VEGFR-3))
    "Muscle" = c("TPM1", "TAGLN", "MYL9"), # TPM = tropomyosin, TAGLN = transgelin (involved crosslinking actin in smooth muscle), MYL = Myosin light chain
    "Mast_cells" = c("KIT", "TPSB2", "HPGD"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase, which are supposed to be common in mast cells
    "HF_surface_markers" = c("ITGB8", "CD200", "SOX9", "LHX2")
)

# Dot plot of cluster markers
count_mat <- GetAssayData(object = obj, slot = "counts")
avgPctMat <- avgAndPctExpressed(count_mat, obj$Clusters, feature_normalize=TRUE, min_pct=5)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")

#wide_df <- prettyOrderMat(wide_df[,rnaOrder], clusterCols=FALSE)
wide_df <- prettyOrderMat(wide_df, clusterCols=TRUE)

grp_order <- colnames(wide_df$mat)
gene_order <- rev(rownames(wide_df$mat))

pdf(paste0(plotDir, "/markers_dot_plot_scalp.pdf"), width=6, height=10)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", xorder=grp_order, yorder=gene_order, cmap=cmaps_BOR$sunrise)
dev.off()
