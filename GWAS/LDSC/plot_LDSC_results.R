#!/usr/bin/env Rscript

##########################################################################
# Make plots of LDSC results
##########################################################################


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(reshape2)
  library(ComplexHeatmap)
})

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin/"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# Load colormaps
broadClustCmap <- readRDS(paste0(scriptPath, "scalpClusterColors.rds")) %>% unlist()
atacNamedClustCmap <- readRDS(paste0(scriptPath, "scATAC_NamedClust_cmap.rds")) %>% unlist()

# Identify all result files
resultDir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/GWAS/ldsc/FineClust_h2_results"
setwd(resultDir)
resultFiles <- list.files(path=resultDir, pattern=".results$", full.names=TRUE)

readH2results <- function(filename){
  # Read in h2 *.results file
  splits <- strsplit(basename(filename), split="\\.")[[1]]
  ct <- splits[1]
  gwas <- paste0(splits[2:(length(splits)-1)], collapse=".")
  # The first row corresponds to the category of interest, with the remaining being components
  # of the baseline model
  res <- fread(filename)[1,] %>% unlist()
  res[1] <- ct
  res <- c(res, "gwas"=gwas)
  res
}

# Do not use 'allPeaks' or clusters that had <10000 specific peaks
exclude <- c("allPeaks", "aMy7", "aKc11")

results <- lapply(resultFiles, readH2results) %>% do.call(rbind,.) %>% as.data.frame()
results <- results[results$Category %ni% exclude,]
results[,2:7] <- sapply(results[,2:7], as.numeric)
results$Enrichment_FDR <- p.adjust(results$Enrichment_p, method="fdr")
results <- results[order(results$Enrichment_FDR, decreasing=FALSE),]

# Heatmap of enrichment with *** indicating FDR thresholds
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

# Add cluster labels
source(paste0(scriptPath, "/cluster_labels.R"))
results$LCategory <- unlist(atac.FineClust)[results$Category]

# Save copy of results
write.table(results, file=paste0(resultDir, "/FineClust_h2_full_results.tsv"), 
  quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE) 

atacOrder <- atacOrder[atacOrder %in% unique(results$Category)]
LatacOrder <- unlist(atac.FineClust)[atacOrder]

mat <- unmelt(results, row_col="gwas", col_col="LCategory", val_col="Enrichment")
pmat <- unmelt(results, row_col="gwas", col_col="LCategory", val_col="Enrichment_FDR")

# Separate negative traits (messes with clustering)
neg_mat <- mat[c("PASS_Schizophrenia", "UKB_460K.body_BMIz"), LatacOrder]
pos_mat <- mat[rownames(mat) %ni% rownames(neg_mat),]

pos_mat <- pos_mat[,LatacOrder]
pos_mat <- prettyOrderMat(pos_mat, clusterCols=FALSE, cutOff=1.0)$mat
mat <- rbind(pos_mat, neg_mat)
pmat <- pmat[rownames(mat), colnames(mat)]

# Get colors for each cluster annotation
inv.FineClust <- invertList(atac.FineClust)
broadClust <- unlist(inv.FineClust)[LatacOrder] %>% gsub('[0-9]+', '', .) %>% sub('.', '', .)
colors <- broadClustCmap[broadClust]
names(colors) <- LatacOrder

pdf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/GWAS/ldsc/h2_results_hm_FineClust.pdf", width=15, height=8)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=colnames(mat),col=list(atac_cluster=colors), 
  show_legend=c("atac_cluster"=FALSE))
hm <- BORHeatmap(
  mat, 
  limits=c(-75,75), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$brewer_yes,
  top_annotation = ta,
  row_names_side = "left",
  width = ncol(mat)*unit(0.5, "cm"),
  height = nrow(mat)*unit(0.5, "cm"),
  legendTitle="Enrichment",
  border_gp=gpar(col="black"), # Add a black border to entire heatmap
  cell_fun=function(j,i,x,y,w,h,col){
    if(pmat[i,j] < 0.0005){
      grid.text("***",x,y)
    }else if(pmat[i,j] < 0.005){
      grid.text("**",x,y)
    }else if(pmat[i,j] < 0.05){
      grid.text("*",x,y)
    }else{
      grid.text("",x,y)
    }
  }
  )
draw(hm)
dev.off()

##########################################################################

