#!/usr/bin/env Rscript

########################################
# Prepare figure panels for figure 1
########################################

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

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/fine_clustered"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_proj <- readRDS("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/preprocessing_output/scalp.rds")
plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

# Color Maps
broadClustCmap <- readRDS(paste0(scriptPath, "/scalpClusterColors.rds")) %>% unlist()
atacNamedClustCmap <- readRDS(paste0(scriptPath, "/scATAC_NamedClust_cmap.rds")) %>% unlist()
rnaNamedClustCmap <- readRDS(paste0(scriptPath, "/scRNA_NamedClust_cmap.rds")) %>% unlist()
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
rna_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(rna_proj$Sample)] %>% unlist()
atac_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(atac_proj$Sample2)] %>% unlist()

# Get label cmaps
source(paste0(scriptPath, "/cluster_labels.R"))
broadLabelCmap <- broadClustCmap
names(broadLabelCmap) <- unlist(BroadClust)[names(broadLabelCmap)]
atacLabelClustCmap <- atacNamedClustCmap
names(atacLabelClustCmap) <- unlist(atac.NamedClust)[names(atacNamedClustCmap)]
rnaLabelClustCmap <- rnaNamedClustCmap
names(rnaLabelClustCmap) <- unlist(rna.NamedClust)[names(rnaNamedClustCmap)]

disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

##########################################################################################
# Dot plots of both RNA and ATAC using common set of markers
##########################################################################################

# Dot plot of marker genes 

# NamedClust first:

# Markers for identifying broad classes of cells:
featureSets <- list(
    "Keratinocytes" = c("KRT5", "KRT10", "KRT14", "KRT15"), 
    "Fibroblasts" = c("THY1", "COL1A1", "COL11A1"), 
    "T_cells" = c("CD3D", "CD8A", "CD4","IKZF2", "CCL5"), 
    "B_cells" = c("CD79A"), 
    "APCs" = c("CD14", "CD86", "CD163", "CD1A", "CLEC9A", "XCR1"), 
    "Melanocytes" = c("MITF", "SOX10", "MLANA"), 
    "Endothlial" = c("VWF", "PECAM1", "SELE"), 
    "Lymphatic" = c("FLT4", "LYVE1"),  
    "Muscle" = c("TPM1", "TAGLN", "MYL9"), 
    "Pericyte" = c("TRPC6", "CCL19"), 
    "Mast_cells" = c("KIT", "TPSB2", "HPGD"), 
    "HF_surface_markers" = c("ITGB8", "CD200", "SOX9")
)

NatacOrder <- c(
    "aTc2", 
    "aTc1", 
    "aTc3", 
    "aBc1", 
    "aMy2", 
    "aMy1", 
    "aMy3", 
    "aKc1", 
    "aKc2", 
    "aKc3", 
    "aKc4", 
    "aKc5", 
    "aKc6", 
    "aKc7", 
    "aFb1", 
    "aFb2", 
    "aMu1", 
    "aMu2", 
    "aVe1", 
    "aVe2", 
    "aLe1", 
    "aMe1" 
)

NrnaOrder <- c(
    "rTc3", 
    "rTc2", 
    "rTc1", 
    "rBc1", 
    "rMy4", 
    "rMy2", 
    "rMy1", 
    "rMy3", 
    "rMa1", 
    "rKc5", 
    "rKc4", 
    "rKc1", 
    "rKc2", 
    "rKc3", 
    "rFb1", 
    "rFb2", 
    "rMu1", 
    "rMu2", 
    "rVe1", 
    "rLe1", 
    "rMe1", 
    "rMe2" 
)
LNatacOrder <- unlist(atac.NamedClust)[NatacOrder]
LNrnaOrder <- unlist(rna.NamedClust)[NrnaOrder]

namedClustAspect <- 1.6
fineClustAspect <- 1.6

# Dot plot of RNA cluster markers
count_mat <- GetAssayData(object=rna_proj, slot="counts")
avgPctMat <- avgAndPctExpressed(count_mat, rna_proj$NamedClust, feature_normalize=TRUE, min_pct=0)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]
avgPctMat <- avgPctMat[avgPctMat$grp %in% NrnaOrder,]

# Assign labels
avgPctMat$grp <- unlist(rna.NamedClust)[as.character(avgPctMat$grp)]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")
wide_df <- prettyOrderMat(wide_df[,LNrnaOrder], clusterCols=FALSE)

grp_order <- colnames(wide_df$mat)
gene_order <- rownames(wide_df$mat) %>% rev() # Reverse this if planning on using plot vertically

pdf(paste0(plotDir, "/RNA_NamedClust_markers_dot_plot_scalp.pdf"), width=6, height=10)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", 
  xorder=grp_order, yorder=gene_order, cmap=cmaps_BOR$sunrise, aspectRatio=namedClustAspect)
dev.off()


# Dot plot of GeneScoreMatrix cluster markers
GSM_se <- getMatrixFromProject(atac_proj, useMatrix="GeneScoreMatrix")
GSM_mat <- assays(GSM_se)$GeneScoreMatrix
rownames(GSM_mat) <- rowData(GSM_se)$name

avgPctMat <- avgAndPctExpressed(GSM_mat[,getCellNames(atac_proj)], atac_proj$NamedClust, feature_normalize=TRUE, min_pct=0)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]
avgPctMat <- avgPctMat[avgPctMat$grp %in% NatacOrder,]

# Assign labels
avgPctMat$grp <- unlist(atac.NamedClust)[as.character(avgPctMat$grp)]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 10] <- 0

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")

# Keep same clustering as from RNA:
wide_df <- wide_df[,LNatacOrder]
grp_order <- colnames(wide_df)

pdf(paste0(plotDir, "/GSM_NamedClust_markers_dot_plot_scalp.pdf"), width=6, height=9)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", 
  xorder=grp_order, yorder=gene_order, cmap=cmaps_BOR$horizonExtra, aspectRatio=namedClustAspect)
dev.off()


# Fine Cluster Markers

# Markers for identifying broad classes of cells:
featureSets <- list(
    "Basal_epithelia" = c("KRT15", "KRT5", "COL17A1"),
    "Spinous" = c("KRT1"),
    "HF_keratinocytes" = c("KRT75", "SOX9", "LHX2","ITGB8", "RUNX3", "KRT23"),
    "Glandular" = c("KRT7"),
    "T_cells" = c("CD3D", "CD8A", "CD4", "FOXP3", "IKZF2", "IFNG"),
    "Plasma" = c("MS4A1"), 
    "M1_macs" = c("CCL20", "CD80", "CD86"),
    "M2a_macs" = c("CD163", "TGFB2"),
    "TREM2_macs" = c("TREM2", "OSM"),
    "FOLR2_macs" = c("FOLR2"),
    "CD1a1c_DCs" = c("CD1A", "CD1C", "ITGAX", "ITGAM"), 
    "CD1a141_DCs" = c("CLEC9A", "XCR1"), 
    "Mast_cells" = c("KIT", "TPSB2"), 
    "Melanocytes" = c("MITF", "SOX10", "MLANA"), 
    "Endothlial" = c("VWF", "PECAM1", "SELE"), 
    "Lymphatic" = c("FLT4", "LYVE1", "CCL21"),  
    "Angiogenic" = c("SEMA3G"),
    "Muscle" = c("TPM1", "TAGLN", "CCL19"), 
    "Fibroblasts" = c("THY1", "COL1A1"),
    "Dermal_sheath" = c("SOX2", "COL11A1"), 
    "Papillary_dermis" = c("COL6A5", "APCDD1"), # PMID: 29391249
    "Reticular_dermis" = c("CD36"), 
    "Dermal_Papilla" = c("BMP7", "HHIP", "PTCH1", "SOX18"),
    "cycling" = c("MKI67", "CDK1", "TOP2A")
)

selectedGenes <- unlist(featureSets) %>% unname()


# Specify order of clusters (Fine Clust)
FrnaOrder <- c(
    # Lymphoid / T-cells
    "rTc3", # "Tc",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG 
    "rTc5", # "NK", # Natural Killer cells (XCL1, XCL2, GNLY, NKG7, KLRD1) 
    "rTc1", # "Th_1", # T-helper: NR3C1, RORA, IL7R, CREM, etc. 
    "rTc2", # "Th_2", # T-helper: JUN, FOS, HSP, CD69 etc. 
    "rTc4", # "Treg", # Regulatory T cells: IKZF2, IL2RA, CTLA4, FOXP3 
    "rTc6", # "Cyc.Tc", # Cycling markers (MKI67, TOP2A, etc.)
    # B/Plasma
    "rBc1", # "Plasma",
    "rMy8", # "Plasma", # (MS4A1, IGHM, CD79A, Stray B-cells...)
    # Myeloid
    "rMy5", # "M1.macs", # IL15, IL32, CCR7 (CCL19, CCL17)
    "rMy6", # "TCR.macs", # CD3, TCR gene positive macrophages
    "rMy2", # "M2.macs_1", # C1Qa/b/c, FOLR2, CD14, CD163, (CCL13)
    "rMy3", # "M2.macs_2", # CXCL2, CXCL3, (CCL20, S100A8/9) 
    "rMy7", # "TREM2.macs", # TREM2
    "rMy1", # "cDC2", # CD1c, CLEC10a (conventional DCs - type 2)
    "rMy4", # "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1
    # Keratinocytes
    "rKc1", # "Basal.Kc_1",
    "rKc2", # "Spinous.Kc_1",
    "rKc3", # "Spinous.Kc_2", 
    "rKc4", # "Spinous.Kc_3",
    "rKc6", # "Cyc.Kc_1", # Cycling keratinocytes
    "rKc5", # "HF.Kc_1", # SOX9, DKK3
    "rKc7", # "HF.Kc_2", # Lhx2, LGR5 high
    "rKc8", # "HF.Kc_3", 
    "rKc9", # "HF.Kc_4", # CD200 high
    "rKc10", # "Glandular_1",
    # Fibroblasts
    "rFb1", # "Fb_1", # CXCL1,2,3
    "rFb3", # "Fb_2", # CCL19, CXCL12
    "rFb4", # "Fb_3", # APCDD1, COL18A1, F13A1
    "rFb5", # "Fb_4", # WISP2, AOX1, ARFGEF3
    "rFb6", # "Fb_5", # NECAB1, SCN7A
    "rFb2", # "D.Sheath", # COL11A1, EDNRA
    "rFb7", # "D.Papilla", # HHIP, PTCH1, etc.
    # Endothelial
    "rVe1", # "Vas.Endo_1",
    "rVe2", # "Vas.Endo_2",
    "rVe3", # "Vas.Endo_3",
    "rVe4", # "Vas.Endo_4",
    "rLe1", # "Lymph.Endo_1",
    # Non-subclustered
    "rMa1", # "Mast", # Mast cells
    "rMu1", # "Muscle_1",
    "rMu2", # "Muscle_2",
    "rMe1", # "Melanocytes",
    "rMe2" # "McSC", # Melanocyte stem cells
    #"Other", # "Other"
    #"rVe5" # "Unknown",
)


# Specify order of clusters (Fine Clust)
FatacOrder <- c(
    # Lymphoid / T-cells
    "aTc3", # "Tc",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG  (Also NK cells absorbed here)
    "aTc1", # "Th_1", # T-helper: NR3C1, RORA, IL7R, CREM, etc. 
    "aTc2", # "Th_2", # T-helper: JUN, FOS, HSP, CD69 etc. 
    "aTc5", # "Th_3",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG 
    "aTc4", # "Treg", # Regulatory T cells: IKZF2, IL2RA, CTLA4, FOXP3
    # B/Plasma
    "aBc1", # "Plasma",
    # Myeloid
    "aMy6", # "M1.macs", # IL15, IL32, CCR7 (CCL19, CCL17)
    "aMy1", # "M2.macs_1",
    "aMy3", # "M2.macs_2", # (Kind of between M2 macs and DCs)
    "aMy4", # "M2.macs_3", # CXCL8 
    "aMy2", # "cDC2_1", # CD1c, CLEC10a (conventional DCs - type 2)
    "aMy7", # "cDC2_2", # IL15, IL32, CCR7 (CCL19, CCL17)
    "aMy5", # "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1
    # Keratinocytes
    "aKc1", # "Basal.Kc_1",
    "aKc2", # "Spinous.Kc_2",
    "aKc3", # "Spinous.Kc_1",
    "aKc4", # "Infundibulum", # SOX9, DKK3
    "aKc5", # "Inf.Segment_1", # Lhx2, LGR5 high
    "aKc7", # "Inf.Segment_2", # Lhx2, LGR5 high
    "aKc6", # "Sebaceous", 
    "aKc8", # "Isthmus", # CD200 high
    "aKc9", # "Matrix",  
    "aKc10", # "Eccrine",
    "aKc11", # "Unknown", # likely doublets
    # Fibroblasts
    "aFb1", # "Fb_1", 
    "aFb2", # "Fb_2", 
    "aFb4", # "Fb_3", # NECAB1, SCN7A (rFb5)
    "aFb5", # "Fb_4", 
    "aFb3", # "D.Sheath", # COL11A1
    "aFb6", # "D.Papilla", # HHIP, PTCH1, etc.
    # Endothelial
    "aVe1", # "Vas.Endo_1",
    "aVe2", # "Vas.Endo_2", 
    "aVe3", # "Vas.Endo_3",
    "aLe1", # "Lymph.Endo_1",
    # Non-subclustered
    "aMu1", # "Muscle",
    "aMu2", # "Pericytes",
    "aMe1" # "Melanocytes",
    #"aVe4" # "Unknown",
    #"Other" # "Other"
)
LFatacOrder <- unlist(atac.FineClust)[FatacOrder]
LFrnaOrder <- unlist(rna.FineClust)[FrnaOrder]

# Dot plot of cluster markers
avgPctMat <- avgAndPctExpressed(count_mat, rna_proj$FineClust, feature_normalize=TRUE, min_pct=5)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]
avgPctMat <- avgPctMat[avgPctMat$grp %in% FrnaOrder,]

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")
wide_df <- prettyOrderMat(wide_df[,FrnaOrder], clusterCols=FALSE)

# Assign labels
avgPctMat$grp <- unlist(rna.FineClust)[as.character(avgPctMat$grp)]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0

grp_order <- colnames(wide_df$mat)
gene_order <- rev(rownames(wide_df$mat))

pdf(paste0(plotDir, "/RNA_FineClust_markers_dot_plot_scalp.pdf"), width=12, height=12)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", 
  xorder=unlist(rna.FineClust)[grp_order], yorder=gene_order, cmap=cmaps_BOR$sunrise, aspectRatio=fineClustAspect)
dev.off()


# Dot plot of GeneScoreMatrix cluster markers
avgPctMat <- avgAndPctExpressed(GSM_mat[,getCellNames(atac_proj)], atac_proj$FineClust, feature_normalize=TRUE, min_pct=5)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]

# Determine cluster and gene order (Keep same feature order as from RNA:)
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")
wide_df <- wide_df[,FatacOrder]

# Assign labels
avgPctMat$grp <- unlist(atac.FineClust)[as.character(avgPctMat$grp)]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 10] <- 0

grp_order <- colnames(wide_df)

pdf(paste0(plotDir, "/GSM_FineClust_markers_dot_plot_scalp.pdf"), width=12, height=10)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", 
  xorder=unlist(atac.FineClust)[grp_order], yorder=gene_order, cmap=cmaps_BOR$horizonExtra, aspectRatio=fineClustAspect)
dev.off()


#############################################################################
# Stacked Bar Plots of Sample by Cluster
#############################################################################

broadOrder <- c(
    "Tc", # "Lymphoid",
    "Bc", # "B-cells" # Or plasma
    "My", # "Myeloid",
    "Ma", # "Mast",
    "Kc", # "Keratinocytes",
    "Fb", # "Fibroblasts",
    "Ve", # "Vascular",
    "Le", # "Lymphatic",
    "Mu", # "Muscle", # And pericytes?
    "Me" # "Melanocytes",
)

rnaBroadOrder <- broadOrder[broadOrder %in% unique(rna_proj$BroadClust)]
atacBroadOrder <- broadOrder[broadOrder %in% unique(atac_proj$BroadClust)]
LrnaBroadOrder <- unlist(BroadClust)[rnaBroadOrder]
LatacBroadOrder <- unlist(BroadClust)[atacBroadOrder]

barWidth <- 0.9

# RNA

# BroadClust by Sample barplot
clustBySamp <- fractionXbyY(unlist(BroadClust)[rna_proj$BroadClust], rna_proj$Sample, add_total=TRUE, xname="Cluster", yname="Sample")
clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels=c(LrnaBroadOrder, "total"), ordered=TRUE)
pdf(paste0(plotDir, "/rna_sampleByBroadClustBarPlot.pdf"), height=5, width=6)
stackedBarPlot(clustBySamp, cmap=rna_sample_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# NamedClust by Sample barplot
clustBySamp <- fractionXbyY(unlist(rna.NamedClust)[rna_proj$NamedClust], rna_proj$Sample, add_total=TRUE, xname="Cluster", yname="Sample")
clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels=c(LNrnaOrder, "total"), ordered=TRUE)
pdf(paste0(plotDir, "/rna_sampleByNamedClustBarPlot.pdf"), height=5, width=8)
stackedBarPlot(clustBySamp, cmap=rna_sample_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# FineClust by Sample barplot
clustBySamp <- fractionXbyY(unlist(rna.FineClust)[rna_proj$FineClust], rna_proj$Sample, add_total=TRUE, xname="Cluster", yname="Sample")
clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels=c(LFrnaOrder, "total"), ordered=TRUE)
pdf(paste0(plotDir, "/rna_sampleByFineClustBarPlot.pdf"), height=5, width=12)
stackedBarPlot(clustBySamp, cmap=rna_sample_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# Sample by NamedClust barplot
clustBySamp <- fractionXbyY(rna_proj$Sample, unlist(rna.NamedClust)[rna_proj$NamedClust], add_total=TRUE, xname="Sample", yname="Cluster")
clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels=c(LNrnaOrder, "total"), ordered=TRUE)
pdf(paste0(plotDir, "/rna_NamedClustBySampleBarPlot.pdf"), height=5, width=8)
stackedBarPlot(clustBySamp, cmap=rnaLabelClustCmap, namedColors=TRUE, barwidth=barWidth)
dev.off()


# ATAC

# BroadClust by Sample barplot
clustBySamp <- fractionXbyY(unlist(BroadClust)[atac_proj$BroadClust], atac_proj$Sample2, add_total=TRUE, xname="Cluster", yname="Sample")
clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels=c(LatacBroadOrder, "total"), ordered=TRUE)
pdf(paste0(plotDir, "/atac_sampleByBroadClustBarPlot.pdf"), height=5, width=6)
stackedBarPlot(clustBySamp, cmap=atac_sample_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# NamedClust by Sample barplot
clustBySamp <- fractionXbyY(unlist(atac.NamedClust)[atac_proj$NamedClust], atac_proj$Sample2, add_total=TRUE, xname="Cluster", yname="Sample")
clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels=c(LNatacOrder, "total"), ordered=TRUE)
pdf(paste0(plotDir, "/atac_sampleByNamedClustBarPlot.pdf"), height=5, width=8)
stackedBarPlot(clustBySamp, cmap=atac_sample_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# FineClust by Sample barplot
clustBySamp <- fractionXbyY(unlist(atac.FineClust)[atac_proj$FineClust], atac_proj$Sample2, add_total=TRUE, xname="Cluster", yname="Sample")
clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels=c(LFatacOrder, "total"), ordered=TRUE)
pdf(paste0(plotDir, "/atac_sampleByFineClustBarPlot.pdf"), height=5, width=12)
stackedBarPlot(clustBySamp, cmap=atac_sample_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# Sample by NamedClust barplot
clustBySamp <- fractionXbyY(atac_proj$Sample2, unlist(atac.NamedClust)[atac_proj$NamedClust], add_total=TRUE, xname="Sample", yname="Cluster")
clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels=c(LNatacOrder, "total"), ordered=TRUE)
pdf(paste0(plotDir, "/atac_NamedClustBySampleBarPlot.pdf"), height=5, width=8)
stackedBarPlot(clustBySamp, cmap=atacLabelClustCmap, namedColors=TRUE, barwidth=barWidth)
dev.off()


##########################################################################################
# Identifying Marker Peaks
##########################################################################################

# NamedClust First:
atac_proj$LNamedClust <- unlist(atac.NamedClust)[atac_proj$NamedClust]
use_groups <- unique(atac_proj$LNamedClust)

# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "LNamedClust",
    useGroups = use_groups,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

#Visualize Markers as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks[,LNatacOrder], 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  nLabel = 1, # It still seems like there's not actually a way to NOT plot any labels
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  transpose = FALSE
)
draw(heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapPeaks, name="Peak-Marker-Heatmap-NamedClust", width=10, height=15, ArchRProj=atac_proj, addDOC=FALSE)

##########################################################################################
# Motif Enrichments
##########################################################################################

#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = atac_proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Rename motifs for more aesthetic plotting:
rownames(enrichMotifs) <- lapply(rownames(enrichMotifs), function(x) strsplit(x, "_")[[1]][1]) %>% unlist()

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
    enrichMotifs[,LNatacOrder], 
    n=5, 
    #clusterCols = FALSE, # This is currently bugged and does not work
    transpose=FALSE, 
    cutOff=log10pCut
)

draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapEM, name="Motifs-Enriched-Heatmap-NamedClust", width=8, height=12, ArchRProj=atac_proj, addDOC=FALSE)

# Let's plot a different style heatmap
plot_mat <- plotEnrichHeatmap(enrichMotifs[,LNatacOrder], n=8, transpose=FALSE, 
  cutOff=log10pCut, returnMatrix=TRUE)
plot_mat <- plot_mat[!grepl("^ENSG", rownames(plot_mat)),]

# Get colors for cluster annotation
# (Link FineClusts to BroadClust cmap)
cM <- as.matrix(confusionMatrix(atac_proj$LNamedClust, atac_proj$BroadClust))
map_colors <- apply(cM, 1, function(x)colnames(cM)[which.max(x)])
plotColors <- map_colors[colnames(plot_mat)]

bc_to_nc_map <- invertList(plotColors)

plot_mat <- prettyOrderMat(plot_mat[,LNatacOrder], clusterCols=FALSE, cutOff=1)$mat

pdf(paste0(plotDir, "/Motifs-Enriched-Heatmap-NamedClust.pdf"), width=8, height=12)
fontsize <- 6
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=plotColors[colnames(plot_mat)],col=list(atac_cluster=broadClustCmap), 
  show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))
hm <- BORHeatmap(
  plot_mat, 
  limits=c(0.0,100.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  top_annotation = ta,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = fontsize),
  column_names_gp = gpar(fontsize = fontsize),
  width = ncol(plot_mat)*unit(0.4, "cm"),
  height = nrow(plot_mat)*unit(0.2, "cm"),
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]",
  border_gp = gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()

##########################################################################################

# FineClust Second:
atac_proj$LFineClust <- unlist(atac.FineClust)[atac_proj$FineClust]
use_groups <- unique(atac_proj$LFineClust)

# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "LFineClust",
    useGroups = use_groups,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

#Visualize Markers as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks[,LFatacOrder], 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  nLabel = 1, # It still seems like there's not actually a way to NOT plot any labels
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  transpose = FALSE
)
draw(heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapPeaks, name="Peak-Marker-Heatmap-FineClust", width=10, height=15, ArchRProj=atac_proj, addDOC=FALSE)

##########################################################################################
# Motif Enrichments
##########################################################################################

#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = atac_proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Rename motifs for more aesthetic plotting:
rownames(enrichMotifs) <- lapply(rownames(enrichMotifs), function(x) strsplit(x, "_")[[1]][1]) %>% unlist()

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
    enrichMotifs[,LFatacOrder], 
    n=5, 
    transpose=FALSE, 
    cutOff=log10pCut
)

draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapEM, name="Motifs-Enriched-Heatmap-FineClust", width=8, height=12, ArchRProj=atac_proj, addDOC=FALSE)

# Let's plot a different style heatmap
plot_mat <- plotEnrichHeatmap(enrichMotifs[,LFatacOrder], n=6, transpose=FALSE, 
  cutOff=log10pCut, returnMatrix=TRUE)
plot_mat <- plot_mat[!grepl("^ENSG", rownames(plot_mat)),]

# Get colors for cluster annotation
# (Link FineClusts to BroadClust cmap)
cM <- as.matrix(confusionMatrix(atac_proj$LFineClust, atac_proj$BroadClust))
map_colors <- apply(cM, 1, function(x)colnames(cM)[which.max(x)])
plotColors <- map_colors[colnames(plot_mat)]

bc_to_fc_map <- invertList(plotColors)

plot_mat <- prettyOrderMat(plot_mat[,LFatacOrder], clusterCols=FALSE, cutOff=1)$mat

pdf(paste0(plotDir, "/Motifs-Enriched-Heatmap-FineClust.pdf"), width=8, height=12)
fontsize <- 6
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=plotColors[colnames(plot_mat)],col=list(atac_cluster=broadClustCmap), 
  show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))
hm <- BORHeatmap(
  plot_mat, 
  limits=c(0.0,100.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  top_annotation = ta,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = fontsize),
  column_names_gp = gpar(fontsize = fontsize),
  width = ncol(plot_mat)*unit(0.3, "cm"),
  height = nrow(plot_mat)*unit(0.2, "cm"),
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]",
  border_gp = gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()

#################################################
# Dot plots of cluster colors
#################################################

# RNA cluster colors
df <- data.frame(clusters=factor(LNrnaOrder, ordered=TRUE, levels=LNrnaOrder), 
    y=rep(1, length(NrnaOrder)))

pdf(paste0(plotDir, "/scRNA_cluster_colors.pdf"), width=10, height=2)
p <- (ggplot(df, aes(x=clusters, y=y, color=clusters))
  + geom_point(size=10)
  + scale_color_manual(values=rnaLabelClustCmap)
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
    panel.grid.minor= element_blank(), 
    plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
    legend.position = "none", # Remove legend
    axis.text.x = element_text(angle = 90, hjust = 1)) 
)
p
dev.off()


# ATAC cluster colors
df <- data.frame(clusters=factor(LNatacOrder, ordered=TRUE, levels=LNatacOrder), 
    y=rep(1, length(NrnaOrder)))

pdf(paste0(plotDir, "/scATAC_cluster_colors.pdf"), width=10, height=2)
p <- (ggplot(df, aes(x=clusters, y=y, color=clusters))
  + geom_point(size=10)
  + scale_color_manual(values=atacLabelClustCmap)
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
    panel.grid.minor= element_blank(), 
    plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
    legend.position = "none", # Remove legend
    axis.text.x = element_text(angle = 90, hjust = 1)) 
)
p
dev.off()


#################################################
# Violin plots of QC metrics
#################################################

rna_ccd <- rna_proj@meta.data

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

# scRNA nUMIs / cell
p <- (
    ggplot(rna_ccd, aes(x=Sample, y=nCount_RNA, fill=Sample))
    + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=rna_sample_cmap, limits=names(rna_sample_cmap), name="Sample", na.value="grey")
    + scale_fill_manual(values=rna_sample_cmap)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("Number of UMIs per cell")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,25000), expand = c(0, 0))
)

pdf(paste0(plotDir, "/scRNA_nUMIs_per_cell_violin.pdf"), width=10, height=4)
p
dev.off()


# scRNA nGenes / cell
p <- (
    ggplot(rna_ccd, aes(x=Sample, y=nFeature_RNA, fill=Sample))
    + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=rna_sample_cmap, limits=names(rna_sample_cmap), name="Sample", na.value="grey")
    + scale_fill_manual(values=rna_sample_cmap)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("Number of genes per cell")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,8000), expand = c(0, 0))
)

pdf(paste0(plotDir, "/scRNA_nGenes_per_cell_violin.pdf"), width=10, height=4)
p
dev.off()


# scRNA pctMito / cell
p <- (
    ggplot(rna_ccd, aes(x=Sample, y=percent.mt, fill=Sample))
    + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=rna_sample_cmap, limits=names(rna_sample_cmap), name="Sample", na.value="grey")
    + scale_fill_manual(values=rna_sample_cmap)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("Percent Mitochondrial Reads")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,50), expand = c(0, 0))
)

pdf(paste0(plotDir, "/scRNA_pctMito_per_cell_violin.pdf"), width=10, height=4)
p
dev.off()

atac_ccd <- atac_proj@cellColData %>% as.data.frame()

# scATAC TSS / cell
p <- (
    ggplot(atac_ccd, aes(x=Sample2, y=TSSEnrichment, fill=Sample2))
    + geom_violin(aes(fill=Sample2), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=atac_sample_cmap, limits=names(atac_sample_cmap), name="Sample", na.value="grey")
    + scale_fill_manual(values=atac_sample_cmap)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("TSS Enrichment")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,25), expand = c(0, 0))
)

pdf(paste0(plotDir, "/scATAC_TSS_per_cell_violin.pdf"), width=10, height=4)
p
dev.off()

# scATAC log10 nFrags / cell
p <- (
    ggplot(atac_ccd, aes(x=Sample2, y=log10nFrags, fill=Sample2))
    + geom_violin(aes(fill=Sample2), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=atac_sample_cmap, limits=names(atac_sample_cmap), name="Sample", na.value="grey")
    + scale_fill_manual(values=atac_sample_cmap)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("log10 nFrags")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,5), expand = c(0, 0))
)

pdf(paste0(plotDir, "/scATAC_log10nFrags_per_cell_violin.pdf"), width=10, height=4)
p
dev.off()

#################################################
