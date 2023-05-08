#!/usr/bin/env Rscript

##########################################################################
# Calculate AUROC for cross-validation gkmSVM
##########################################################################

# See:
# https://github.com/Dongwon-Lee/lsgkm/
# https://github.com/kundajelab/lsgkm-svr
# https://github.com/kundajelab/gkmexplain/blob/master/dsQTL/gm12878_sequence_sets/compute_auroc.py

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(PRROC)
  library(ComplexHeatmap)
})

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/plotting_config.R"))

# Read in cross validation results from gkmSVM train
resultDir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/GWAS/gkmSVM/model_predictions"
plotDir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/GWAS/gkmSVM"

res_files <- list.files(
  path=resultDir, 
  pattern="\\.fold0\\.pred\\..*\\.txt$", 
  full.names=TRUE
  )

# Files have the format <path>/<pred_cluster>.pred.<target_type>.[true/null].txt
clusters <- str_replace(basename(res_files), "\\.fold[0-9]+\\..*\\.txt$", "")
targets <- str_match(basename(res_files), "\\.pred\\.([^\\.]*)\\..*")[,2]
truenull <- str_match(basename(res_files), "\\.(true|null)\\.txt$")[,2]

meta_df <- data.frame(cluster=clusters, target=targets, type=truenull, res_file=res_files)

# Calculate AUROC and AUPRC for each prediction:
getPerformancefromFiles <- function(types, files){
  # Read in gkmSVM cv files and calculate AUROC
  file_df <- data.frame(type=types, res_file=files)
  true_dt <- fread(file_df[file_df$type == "true","res_file"])
  null_dt <- fread(file_df[file_df$type == "null","res_file"])
  # If null is longer than true, sample null sequences for calculating AUROC
  if(nrow(true_dt) < nrow(null_dt)){
    set.seed(1)
    null_dt <- null_dt[sample(1:nrow(null_dt), size=nrow(true_dt), replace=FALSE)]
  }
  AUROC <- PRROC::roc.curve(scores.class0 = true_dt[[2]], scores.class1 = null_dt[[2]])$auc
  AUPRC <- PRROC::pr.curve(scores.class0 = true_dt[[2]], scores.class1 = null_dt[[2]], dg.compute=FALSE)$auc.integral
  data.frame("AUROC"=AUROC, "AUPRC"=AUPRC)
}

# Fairly slow...
summary_df <- meta_df %>% group_by(cluster, target) %>% do(getPerformancefromFiles(.$type, .$res_file)) %>% as.data.frame()

auroc_mat <- unmelt(summary_df, row_col="cluster", col_col="target", val_col="AUROC")
auprc_mat <- unmelt(summary_df, row_col="cluster", col_col="target", val_col="AUPRC")
# mat <- prettyOrderMat(mat,clusterCols=TRUE, cutOff=1.0)$mat

source(paste0(scriptPath, "/cluster_labels.R"))

# Specify order of clusters (Fine Clust)
atacOrder <- c(
 # Lymphoid / T-cells
  "aTc3", # "Tc",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG 
  "aTc1", # "Th_1", # T-helper
  "aTc2", # "Th_2", # T-helper
  "aTc5", # "Th_3", # T-helper
  "aTc4", # "Treg", # Regulatory T cells: IKZF2, IL2RA, CTLA4, FOXP3
  # B/Plasma
  "aBc1", # "Plasma",
  # Myeloid
  "aMy6", # "M1.macs", # IL15, IL32, CCR7
  "aMy1", # "M2.macs_1",
  "aMy3", # "M2.macs_2", 
  "aMy4", # "M2.macs_3", 
  "aMy2", # "cDC2_1", # CD1c, CLEC10a (conventional DCs - type 2)
  #"aMy7", # "cDC2_2", 
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
  #"aKc11", # "Unknown",
  # Fibroblasts
  "aFb1", # "Fb_1", 
  "aFb2", # "Fb_2", 
  "aFb4", # "Fb_3", 
  #"aFb5", # "Fb_4", 
  "aFb3", # "D.Sheath", # COL11A1
  "aFb6", # "D.Papilla", # Many markers HHIP, PTCH1, etc.
  # Endothelial
  "aVe1", # "Vas.Endo_1",
  "aVe2", # "Vas.Endo_2", 
  "aVe3", # "Vas.Endo_3",
  "aLe1", # "Lymph.Endo_1",
  # Non-subclustered
  "aMu1", # "Muscle",
  "aMu2", # "Pericytes", 
  "aMe1" # "Melanocytes"
)
LatacOrder <- unlist(atac.FineClust)[atacOrder]

# Relabel
rownames(auroc_mat) <- unlist(atac.FineClust)[rownames(auroc_mat)]
rownames(auprc_mat) <- unlist(atac.FineClust)[rownames(auprc_mat)]
colnames(auroc_mat) <- unlist(atac.FineClust)[colnames(auroc_mat)]
colnames(auprc_mat) <- unlist(atac.FineClust)[colnames(auprc_mat)]

LatacOrder <- LatacOrder[LatacOrder %in% rownames(auroc_mat)]

# AUROC heatmap
pdf(paste0(plotDir, "/crosstype_AUROC_hm_1000bp.pdf"), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  auroc_mat[LatacOrder, LatacOrder], 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  column_title="Target Cell Type", 
  column_title_side="bottom",
  row_title="Model Cell Type",
  dataColors = cmaps_BOR$solar_extra,
  row_names_side = "left",
  width = ncol(auroc_mat)*unit(0.5, "cm"),
  height = nrow(auroc_mat)*unit(0.5, "cm"),
  legendTitle="AUROC",
  border_gp=gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()

# AUPRC heatmap
pdf(paste0(plotDir, "/crosstype_AUPRC_hm_1000bp.pdf"), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  auprc_mat[LatacOrder, LatacOrder], 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  column_title="Target Cell Type", 
  column_title_side="bottom",
  row_title="Model Cell Type",
  dataColors = cmaps_BOR$solar_extra,
  row_names_side = "left",
  width = ncol(auprc_mat)*unit(0.5, "cm"),
  height = nrow(auprc_mat)*unit(0.5, "cm"),
  legendTitle="AUPRC",
  border_gp=gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()
