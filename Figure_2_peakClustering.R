#!/usr/bin/env Rscript

##########################################################
# Analyses of full project peak to gene linkages
##########################################################

#Load ArchR (and associated libraries)
library(ArchR)
library(igraph)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(ggrastr)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered"
plotDir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/p2gLink_plots"

#Set/Create Working Directory to Folder
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 0.25
barwidth <- 0.9

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)

# Color Maps
broadClustCmap <- readRDS(paste0(scriptPath, "/scalpClusterColors.rds")) %>% unlist()
atacNamedClustCmap <- readRDS(paste0(scriptPath, "/scATAC_NamedClust_cmap.rds")) %>% unlist()
rnaNamedClustCmap <- readRDS(paste0(scriptPath, "/scRNA_NamedClust_cmap.rds")) %>% unlist()
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
atac_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(atac_proj$Sample2)] %>% unlist()

# Get label cmaps
source(paste0(scriptPath, "/cluster_labels.R"))
atacLabelClustCmap <- atacNamedClustCmap
names(atacLabelClustCmap) <- unlist(atac.NamedClust)[names(atacNamedClustCmap)]
rnaLabelClustCmap <- rnaNamedClustCmap
names(rnaLabelClustCmap) <- unlist(rna.NamedClust)[names(rnaNamedClustCmap)]

# Add labels to project
source(paste0(scriptPath, "/cluster_labels.R"))
atac_proj$LNamedClust <- unlist(atac.NamedClust)[atac_proj$NamedClust]

disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

# P2G definition cutoffs
corrCutoff <- 0.5       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Coaccessibility cutoffs
coAccCorrCutoff <- 0.0  # Default in getCoAccessibility is 0.5

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

# Get full peak matrix
pSE <- getMatrixFromProject(atac_proj, useMatrix="PeakMatrix", binarize=FALSE)
pmat <- assays(pSE)$PeakMatrix
rownames(pmat) <- rowRanges(pSE) %>% {paste0(seqnames(.), "_", start(.), "_", end(.))}
pmat <- pmat[,getCellNames(atac_proj)] # Need to get cells in order of ArchR project...

# Get GeneScore matrix
gsmSE <- getMatrixFromProject(atac_proj, useMatrix="GeneScoreMatrix", binarize=FALSE)
gsmat <- assays(gsmSE)$GeneScoreMatrix
rownames(gsmat) <- rowData(gsmSE)$name
gsmat <- gsmat[,getCellNames(atac_proj)] # Need to get cells in order of ArchR project...

# Get GeneIntegration matrix
gimSE <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix", binarize=FALSE)
gimat <- assays(gimSE)$GeneIntegrationMatrix
rownames(gimat) <- rowData(gimSE)$name
gimat <- gimat[,getCellNames(atac_proj)] # Need to get cells in order of ArchR project...

##########################################################################################
# Prepare full-project peak to gene linkages, loops, and coaccessibility (full and subproject links)
##########################################################################################

# Load lists of p2g objects, etc.
full_p2gGR <- readRDS(file=paste0(wd, "/multilevel_p2gGR.rds")) # NOT merged or correlation filtered
full_coaccessibility <- readRDS(file=paste0(wd, "/multilevel_coaccessibility.rds"))
plot_loop_list <- readRDS(file=paste0(wd, "/multilevel_plot_loops.rds"))

##########################################################################################
# Filter redundant peak to gene links
##########################################################################################

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(atac_proj@peakSet)$Peak2GeneLinks # (Save original P2G links just in case)
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$symbol, "-", full_p2gGR$peakName))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF

##########################################################################################
# Identify 'highly regulated' genes and 'highly-regulating' peaks
##########################################################################################

# Get distribution of peaks to gene linkages and identify 'highly-regulated' genes
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)
p2gFreqs <- getFreqs(p2gGR$symbol)

x <- 1:length(p2gFreqs)
rank_df <- data.frame(npeaks=p2gFreqs, rank=x)

##########################################################################################
# Get pseudo-bulked data for downstream analyses
##########################################################################################

ccd <- atac_proj@cellColData

# Get pseudo-bulks
knn_groups <- getLowOverlapAggregates(atac_proj, target.agg=250, k=250, 
  overlapCutoff=0.8, dimReduc="IterativeLSI", seed=123)

kgrps <- unique(knn_groups$group)

# Get pseudo-bulked peak matrix
pmatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowSums(pmat[,use_cells])
  }) %>% do.call(cbind,.)
colnames(pmatPsB) <- kgrps

# Get pseudo-bulked GeneScore matrix
gsmatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(gsmat[,use_cells]) # GImatrix is already scaled
  }) %>% do.call(cbind,.)
colnames(gsmatPsB) <- kgrps

# Get pseudo-bulked GeneIntegration matrix
gimatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(gimat[,use_cells]) # GImatrix is already scaled
  }) %>% do.call(cbind,.)
colnames(gimatPsB) <- kgrps

# Get scaled peak matrix (used for summing accessibility across linked peaks)
scaled_pmatPsB <- t(t(pmatPsB)/colSums(pmatPsB)) * 10000 #(The peak matrix is not normalized prior to this)

# Get log transformed GeneIntegration matrix
L2gsmatPsB <- sparseLogX(gsmatPsB, logtype="log2", scale=FALSE) #(The GeneScore matrix is already scaled)
L2gimatPsB <- sparseLogX(gimatPsB, logtype="log2", scale=FALSE) #(The GeneIntegration matrix is already scaled)

# Identify likely cluster based on majority label
psb_labels <- sapply(kgrps, function(kg){
  clustLabels <- ccd[knn_groups$cell_name[knn_groups$group == kg],]$NamedClust
  names(sort(table(clustLabels),decreasing=TRUE))[1]
  })

# Normalize PseudoBulked peaks (used only for plotting heatmaps)
norm_pmatPsB <- sparseLogX(pmatPsB, logtype="log2", scale=TRUE, scaleFactor=10^6) %>% preprocessCore::normalize.quantiles()
colnames(norm_pmatPsB) <- colnames(pmatPsB)
rownames(norm_pmatPsB) <- rownames(pmatPsB)

# Row-scale peak matrix
centered_pmatPsB <- t(scale(t(norm_pmatPsB), center=TRUE, scale=TRUE))


#################################################################################
# H-cluster peaks and explore relationship between number of peaks and resulting expression
#################################################################################

names(p2gGR) <- p2gGR$peakName
hm_outdir <- paste0(plotDir, "/clustering_linkedPeaks")
dir.create(hm_outdir, showWarnings=FALSE, recursive=TRUE)

getLinkedPeakMat <- function(gene, p2gGR, pmat){
  # Return matrix of clustered linked peaks for a given gene
  ################################
  # gene = name of gene to build graph for
  # p2gGR = full peak to gene GR (with peakNames)
  # coaccGR = coaccessibility GR
  # pmat = pseudo-bulked matrix of peak accessibility

  # Get peaks linked to gene
  peaks <- p2gGR[p2gGR$symbol == gene]$peakName
  peak_mat <- pmat[peaks,]

  # hierarchical cluster peaks and identify regulatory 'modules'
  hc <- hclust(dist(peak_mat), method="complete", members=NULL)
  clusters <- cutree(hc, h=23) # Cutoff of 25 yields 1-5 clusters per gene

  return(list(mat=peak_mat[hc$order,], clusters=clusters[hc$order], hc=hc))
}


# Get heatmap labeled by cluster membership
gene_set <- c("RUNX3", "RUNX1", "AQP3", "HLA-DRB1")

for(gene in gene_set){
  gene_pmat_list <- getLinkedPeakMat(gene, p2gGR, centered_pmatPsB)
  hmat <- gene_pmat_list$mat
  clusters <- gene_pmat_list$clusters
  cmap <- getColorMap(cmaps_BOR$stallion, n=length(unique(clusters)))
  names(cmap) <- names(getFreqs(clusters))

  # Get integrated gene expression
  expr_df <- data.frame(
    geneExpr=L2gimatPsB[gene,],
    geneActivity=L2gsmatPsB[gene,],
    peakSum=log2(colSums(scaled_pmatPsB[rownames(hmat),])+1), # log2 of sum of scaled peaks
    cellType=psb_labels
  )
  peak_rsq <- cor.test(expr_df$peakSum, expr_df$geneExpr)$estimate**2
  gs_rsq <- cor.test(expr_df$geneActivity, expr_df$geneExpr)$estimate**2

  yrange <- max(expr_df$geneExpr) - min(expr_df$geneExpr)
  plot_ylims <- c(min(expr_df$geneExpr) - 0.05*yrange, max(expr_df$geneExpr) + 0.1*yrange)
  xrange <- max(expr_df$peakSum) - min(expr_df$peakSum)
  plot_xlims <- c(min(expr_df$peakSum) - 0.05*xrange, max(expr_df$peakSum) + 0.1*xrange)

  pointSize <- 2.5
  p1 <- (
    ggplot(expr_df, aes(x=peakSum, y=geneExpr, color=cellType))
    + geom_point(size=pointSize)
    + geom_smooth(method='lm', color="red", se=FALSE)
    + scale_color_manual(values=atacNamedClustCmap, limits=names(atacNamedClustCmap))
    + ylab("Log2(Integrated Gene Expression)")
    + xlab("Log2(Linked Peak Accessibility)")
    + theme_BOR(border=FALSE)
    + ggtitle(sprintf("%s, Peak Sum Rsq = %s", gene, round(peak_rsq, 3)))
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 1,
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(limits=plot_ylims, expand = c(0, 0))
    + scale_x_continuous(limits=plot_xlims, expand = c(0, 0))
  )
  xrange <- max(expr_df$geneActivity) - min(expr_df$geneActivity)
  plot_xlims <- c(min(expr_df$geneActivity) - 0.1*xrange, max(expr_df$geneActivity) + 0.1*xrange)

  p2 <- (
    ggplot(expr_df, aes(x=geneActivity, y=geneExpr, color=cellType))
    + geom_point(size=pointSize)
    + geom_smooth(method='lm', color="red", se=FALSE)
    + scale_color_manual(values=atacNamedClustCmap, limits=names(atacNamedClustCmap))
    + ylab("Log2(Integrated Gene Expression)")
    + xlab("Log2(GeneActivityScore)")
    + theme_BOR(border=FALSE)
    + ggtitle(sprintf("%s, GeneActivity Rsq = %s", gene, round(gs_rsq, 3)))
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 1,
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(limits=plot_ylims, expand = c(0, 0))
    + scale_x_continuous(limits=plot_xlims, expand = c(0, 0))
  )

  # Determine cluster order from accessibility (already set)
  peak_order <- rownames(hmat)
  exprRange <- quantile(expr_df$geneExpr, c(0, 0.95))
  use_cmap <- cmaps_BOR$sunrise
  col_fun <- circlize::colorRamp2(seq(exprRange[1], exprRange[2], length=length(use_cmap)), use_cmap)

  pdf(paste0(hm_outdir, sprintf("/hclust_linked_peaks_%s.pdf", gene)), width=7, height=6)
  ht_opt$simple_anno_size <- unit(0.25, "cm")
  ta <- HeatmapAnnotation(cellType=expr_df$cellType, geneExpr=expr_df$geneExpr,
    col=list(cellType=atacNamedClustCmap, geneExpr=col_fun))
  la <- HeatmapAnnotation(cluster=clusters[rownames(hmat)],col=list(cluster=cmap), 
    which="row", show_legend=c("cluster"=TRUE))
  hm <- BORHeatmap(
    hmat[peak_order,],
    limits=c(-2.,2.),
    clusterCols=TRUE, clusterRows=FALSE,
    labelCols=FALSE, labelRows=FALSE,
    dataColors = cmaps_BOR$solar_extra,
    top_annotation = ta,
    left_annotation = la,
    row_names_side = "left",
    width = ncol(hmat)*unit(0.05, "cm"),
    height = nrow(hmat)*unit(0.1, "cm"),
    border_gp=gpar(col="black"), # Add a black border to entire heatmap
    legendTitle="Row Z-score"
    )
  draw(hm)
  print(p1)
  print(p2)
  dev.off()
}


# Compare linked peak rsq vs gene score rsq for all highly linked genes
top_genes <- rownames(rank_df[rank_df$npeaks > 20,])
cor_res_list <- list()
for(gene in top_genes){

  gene_pmat_list <- getLinkedPeakMat(gene, p2gGR, centered_pmatPsB)
  hmat <- gene_pmat_list$mat

  # Get integrated gene expression
  expr_df <- data.frame(
    geneExpr=L2gimatPsB[gene,],
    geneActivity=L2gsmatPsB[gene,],
    peakSum=log2(colSums(scaled_pmatPsB[rownames(hmat),])+1)
  )

  peakSumFit <- lm(formula = expr_df$geneExpr ~ expr_df$peakSum)
  slope <- summary(peakSumFit)$coefficients[2,1]
  intercept <- summary(peakSumFit)$coefficients[1,1]

  peak_rsq <- cor.test(expr_df$peakSum, expr_df$geneExpr)$estimate**2
  gs_rsq <- cor.test(expr_df$geneActivity, expr_df$geneExpr)$estimate**2
  cor_res_list[[gene]] <- c(peakSumSlope=slope, peakSumIntercept=intercept, peak_rsq=peak_rsq, gs_rsq=gs_rsq)
}

cor_df <- as.data.frame(do.call(rbind, cor_res_list))
cor_df$npeaks <- rank_df[rownames(cor_df),]$npeaks


melt_df <- reshape2::melt(cor_df[,c(3,4)])

p <- (
  ggplot(melt_df, aes(x=variable, y=value, color=variable, fill=variable))
  + geom_violin(alpha=0.5)
  + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
  #+ geom_jitter(aes(group=variable), size=1.0, color="black",
  #  position=position_jitterdodge(seed=1, jitter.width=0.25, jitter.height=0.0, dodge.width=0.75))
  + scale_fill_manual(values=cmaps_BOR$stallion)
  + scale_color_manual(values=cmaps_BOR$stallion)
  #+ xlab(colnames(df)[1])
  + ylab("R-squared to integrated expression")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          #aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
          axis.text.x = element_text(angle = 90, hjust = 1)) 
  + scale_y_continuous(limits=c(0,1.0), expand = c(0, 0))
)

pdf(paste0(hm_outdir, "/peakScore_vs_geneScore_topGenes.pdf"), width=6, height=6)
p
dev.off()

#################################################################################