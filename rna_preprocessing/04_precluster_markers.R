#!/usr/bin/env Rscript

############################################
# Identify marker genes and make some plots
############################################

library(dplyr)
library(tidyr)
library(Seurat)
library(ggrastr)
library(Rmagic)
library(future) # For parallelization
library(data.table)

subgroup <- "preclustered"
pointSize <- 0.2
useMagic <- TRUE # Should Rmagic be used for data imputation prior to UMAP plotting?

# change the current plan to access parallelization (for Seurat)
nThreads <- 8
plan("multicore", workers=nThreads)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin/"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scRNA_preprocessing/preprocessing_output"
plotDir <- paste0(wd,"/expression_plots")
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
# Identify markers per cluster
##########################################

# find markers for every cluster compared to all remaining cells, report only the positive ones
message("Finding marker genes using Seurat...")
obj.markers <- FindAllMarkers(
    obj, 
    only.pos=TRUE, min.pct=0.25, logfc.threshold=1.0,
    max.cells.per.ident = 250 # Don't bother using too many cells for this step
    )

# Save markers
write.table(obj.markers, file=paste0(wd, sprintf("/marker_genes_%s.tsv", subgroup)), 
    sep="\t", quote=FALSE,row.names=FALSE)

##########################################
# UMAP of high level groups
##########################################
message("Plotting selected marker features on UMAP...")

# Set colormaps
qualcmap <- cmaps_BOR$stallion
quantcmap <- cmaps_BOR$sunrise

# Markers for identifying broad classes of cells:
featureSets <- list(
    "Keratinocytes" = c("KRT5", "KRT10", "KRT14", "KRT15", "KRT6A", "KRT6B"), 
    "Fibroblasts" = c("THY1", "COL1A1", "COL1A2"), 
    "T_cells" = c("CD3D", "CD8A", "CD4", "FOXP3", "IKZF2", "CCL5"), # PTPRC = CD45 
    "B_cells" = c("CD19", "MS4A1", "MZB1"), # MS41A = CD20, MZB1 = marginal zone B cell specific protein
    "APCs" = c("CD86", "CD74", "CCR7", "CD163"), # Monocyte lineage (FCGR3A = CD16, FCGR1A = CD64, CD74 = HLA-DR antigens-associated invariant chain)
    "Melanocytes" = c("MITF", "TYR", "SOX10", "MLANA"), # Melanocyte markers
    "Endothlial" = c("VWF", "PECAM1", "SELE"), # Endothlial cells (PECAM1 = CD31), SELE = selectin E (found in cytokine stimulated endothelial cells)
    "Lymphatic" = c("ACKR2", "FLT4", "LYVE1"),  # Lymphatic endothelial (FLT4 = VEGFR-3))
    "Muscle" = c("TPM1", "TAGLN", "MYL9"), # TPM = tropomyosin, TAGLN = transgelin (involved crosslinking actin in smooth muscle), MYL = Myosin light chain
    "Mast_cells" = c("KIT", "FCER1A", "IL1RL1", "TPSB2"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase, which are supposed to be common in mast cells
    "Langerhans_cells" = c("ITGAX", "CD1A", "CLEC1A", "CD207"), # CD11c = ITGAX, Langerin = CD207, CLEC4K
    "HF_surface_markers" = c("GJB6", "ITGB8", "CD200", "FZD7")
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

# Dot plot of cluster markers
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

pdf(paste0(plotDir, "/markers_dot_plot_preclustered.pdf"), width=6, height=10)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", xorder=grp_order, yorder=gene_order, cmap=cmaps_BOR$sunrise)
dev.off()




