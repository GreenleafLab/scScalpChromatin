#!/usr/bin/env Rscript

##########################################################################
# Calculate AUROC and AUPRC for cross-validation gkmSVM
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
scriptPath <- "/home/users/boberrey/git_clones/hairATAC"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/plotting_config.R"))

# Read in cross validation results from gkmSVM train
resultDir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/model_predictions"
plotDir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM"

# Use capture group to get only files that predict their own cell type
res_files <- list.files(
  path=resultDir, 
  pattern="^(.*\\.)fold.*\\.pred\\.\\1.*\\.txt$", 
  full.names=TRUE
  )

# Files have the format <path>/<cluster>.<fold#>.pred.<target_cell_type>.txt
clusters <- str_replace(basename(res_files), "\\.fold[0-9]+\\..*\\.txt$", "")
fold <- str_extract(basename(res_files), "\\.fold[0-9]+\\.") %>% str_extract(.,"[0-9]+") %>% as.numeric()
truenull <- str_match(basename(res_files), "\\.(true|null)\\.txt$")[,2]
meta_df <- data.frame(cluster=clusters, fold=fold, type=truenull, res_file=res_files)


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

summary_df <- meta_df %>% group_by(cluster, fold) %>% do(getPerformancefromFiles(as.character(.$type), as.character(.$res_file))) %>% as.data.frame()

# Translate to labeled clusters
source(paste0(scriptPath, "/cluster_labels.R"))
summary_df$Lcluster <- unlist(atac.NamedClust)[summary_df$cluster]

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
  "aMy5", # "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1 https://www.frontiersin.org/articles/10.3389/fimmu.2014.00239/full
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

summary_df$Lcluster <- factor(summary_df$Lcluster, levels=LatacOrder, ordered=TRUE)
summary_df <- summary_df[summary_df$Lcluster %in% LatacOrder,]

# Map FineClust colors to NamedClust
library(ArchR)
atac_proj <- loadArchRProject("/scratch/groups/wjg/boberrey/hairATAC/analyses/scATAC_preprocessing/fine_clustered", force=TRUE)
cM <- as.matrix(confusionMatrix(atac_proj$FineClust, atac_proj$NamedClust))
allFineClust <- unique(atac_proj$FineClust)
map_colors <- apply(cM, 1, function(x)colnames(cM)[which.max(x)])
atacNamedClustCmap <- readRDS(paste0(scriptPath, "/scATAC_NamedClust_cmap.rds")) %>% unlist()
cmap <- atacNamedClustCmap[map_colors[allFineClust]]
names(cmap) <- unlist(atac.FineClust)[allFineClust]

# Plot a violin / box plot
dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

# Plot AUROC
p <- (
  ggplot(summary_df, aes(x=Lcluster, y=AUROC, color=Lcluster, fill=Lcluster))
  + geom_boxplot(alpha=0.5, outlier.shape=NA) # Hide fliers (we show them with geom_jitter)
  + geom_jitter(aes(group=cluster), size=0.75, color="black",
     position=position_jitterdodge(seed=1, jitter.width=7.0, jitter.height=0.0, dodge.width=dodge_width))
  + scale_y_continuous(limits=c(0.5,1.0), expand=c(0,0))
  + scale_color_manual(values=cmap)
  + scale_fill_manual(values=cmap)
  + xlab("")
  + ylab("AUROC")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          legend.position="none", # Remove legend
          axis.text.x = element_text(angle=90, hjust=1)) 
)

pdf(paste0(plotDir, "/cv_AUROC_boxplot_1000bpNC_fullLim.pdf"), width=8, height=6)
p
dev.off()

# Plot AUPRC
p <- (
  ggplot(summary_df, aes(x=Lcluster, y=AUPRC, color=Lcluster, fill=Lcluster))
  + geom_boxplot(alpha=0.5, outlier.shape=NA) # Hide fliers (we show them with geom_jitter)
  + geom_jitter(aes(group=cluster), size=0.75, color="black",
     position=position_jitterdodge(seed=1, jitter.width=7.0, jitter.height=0.0, dodge.width=dodge_width))
  + scale_y_continuous(limits=c(0.5,1.0), expand=c(0,0))
  + scale_color_manual(values = cmap)
  + scale_fill_manual(values = cmap)
  + xlab("")
  + ylab("AUPRC")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          legend.position="none", # Remove legend
          axis.text.x = element_text(angle=90, hjust=1)) 
)

pdf(paste0(plotDir, "/cv_AUPRC_boxplot_1000bpNC_fullLim.pdf"), width=8, height=6)
p
dev.off()

###################################################################################