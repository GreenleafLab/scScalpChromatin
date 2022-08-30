#!/usr/bin/env Rscript

##########################################################
# Analyses of full project peak to gene linkages
##########################################################

#Load ArchR (and associated libraries)
library(ArchR)
library(Seurat)
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
rna_proj <- readRDS("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scRNA_preprocessing/preprocessing_output/scalp.rds")

# Color Maps
broadClustCmap <- readRDS(paste0(scriptPath, "/scalpClusterColors.rds")) %>% unlist()
atacNamedClustCmap <- readRDS(paste0(scriptPath, "/scATAC_NamedClust_cmap.rds")) %>% unlist()
rnaNamedClustCmap <- readRDS(paste0(scriptPath, "/scRNA_NamedClust_cmap.rds")) %>% unlist()
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
rna_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(rna_proj$Sample)] %>% unlist()
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
corrCutoff <- 0.5      # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})

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
p2gMeta <- metadata(atac_proj@peakSet)$Peak2GeneLinks %>% metadata()

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$symbol, "-", full_p2gGR$peakName))]

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF


atacOrder <- c(
    "aTc2", # "CD8.Tc"
    "aTc1", # "CD4.Tc"
    "aTc3", # "Tregs"
    "aBc1", # "B.cells"
    "aMy2", # "Macs_1"
    "aMy1", # "DCs_1"
    "aMy3", # "CLEC9a.DC"
    "aKc1", # "Basal.Kc_1"
    "aKc2", # "Spinous.Kc_1"
    "aKc3", # "Spinous.Kc_2"
    "aKc4", # "HF.Kc_1"
    "aKc5", # "HF.Kc_2"
    "aKc6", # "HF.Kc_3",
    "aKc7", # "HF.Kc_4",
    "aFb1", # "D.Fib" # Papillary/Reticular dermal fibroblasts
    "aFb2", # "D.Sheath" # Dermal sheath
    "aMu1", # "Muscle_1"
    "aMu2", # "Muscle_2" # Myofibroblasts?
    "aVe1", # "Vas.Endo_1"
    "aVe2", # "Vas.Endo_2"
    "aLe1", # "Lymph.Endo"
    "aMe1" # "Melanocytes"
)

label_genes <- c(
  "ICOS", "RUNX3", "TWIST2", "CD84", "CTLA4", "KRT14", "IKZF1", "COL1A1",
  "HLA-DRB1", "CD28", "EGFR", "CD3D", "ITGAX", "CXCR6",
  "TNF", "RUNX1", "MITF", "FOSL2", "FZD7", "MALAT1", "POU2F3"
  )

# Tracks of genes:
# (Define plot region based on bracketing linked peaks)
promoterGR <- promoters(getGenes(atac_proj))

# markerGenes <- c("IL21", "RUNX3")
mPromoterGR <- promoterGR[promoterGR$symbol %in% label_genes]
mP2G_GR <- p2gGR[p2gGR$symbol %in% label_genes]

# Restrict to only loops linking genes of interest (full project loops)
plotLoops <- getPeak2GeneLinks(atac_proj, corCutOff=corrCutoff, resolution = 100)[[1]]
sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
plotLoops <- plotLoops[width(plotLoops) > 100]

# Create copy of individual project plot loops
sub_plot_loop_list <- list()
for(pn in names(plot_loop_list)){
  subPlotLoops <- plot_loop_list[[pn]]
  sol <- findOverlaps(resize(subPlotLoops, width=1, fix="start"), mPromoterGR)
  eol <- findOverlaps(resize(subPlotLoops, width=1, fix="end"), mPromoterGR)
  subPlotLoops <- c(subPlotLoops[from(sol)], subPlotLoops[from(eol)])
  subPlotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
  sub_plot_loop_list[[pn]] <- subPlotLoops[width(subPlotLoops) > 100]
}

# Bracket plot regions around loops
plotRegions <- lapply(label_genes, function(x){
  gr <- range(plotLoops[plotLoops$symbol == x])
  lims <- grLims(gr)
  gr <- GRanges(
      seqnames = seqnames(gr)[1],
      ranges = IRanges(start=lims[1], end=lims[2])
    )
  gr
  }) %>% as(., "GRangesList") %>% unlist()
plotRegions <- resize(plotRegions, 
  width=width(plotRegions) + 0.05*width(plotRegions), 
  fix="center")


# Tracks of genes (scalp):
p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "LNamedClust", 
    useGroups = unlist(atac.NamedClust)[atacOrder],
    pal = atacLabelClustCmap,
    plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
    sizes = c(7, 0.2, 1.25, 2.5),
    geneSymbol = label_genes, 
    region = plotRegions, 
    loops = sub_plot_loop_list$scalp,
    tileSize=500,
    minCells=200
)

plotPDF(plotList = p, 
    name = "Super-Enhancer-Tracks-scalpOnly-p2gLinks.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 6, height = 7)


# Tracks of genes (Keratinocytes):
p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "LNamedClust", 
    useGroups = unlist(atac.NamedClust)[atacOrder],
    pal = atacLabelClustCmap,
    plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
    sizes = c(7, 0.2, 1.25, 2.5),
    geneSymbol = label_genes, 
    region = plotRegions, 
    loops = sub_plot_loop_list$Keratinocytes,
    tileSize=500,
    minCells=200
)

plotPDF(plotList = p, 
    name = "Super-Enhancer-Tracks-KeratinocytesOnly-p2gLinks.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 6, height = 7)


# Broad Cluster Tracks:
broadOrder <- c(
    "Tc", # "Lymphoid",
    "My", # "Myeloid",
    "Bc", # "B-cells" # Or plasma?
    "Ma", # "Mast",
    "Ve", # "Vascular",
    "Le", # "Lymphatic",
    "Fb", # "Fibroblasts",
    "Mu", # "Muscle", # And pericytes?
    "Me", # "Melanocytes",
    "Kc" # "Keratinocytes",
)
atacBroadOrder <- broadOrder[broadOrder %in% unique(atac_proj$BroadClust)]
LatacBroadOrder <- unlist(BroadClust)[atacBroadOrder]

# Tracks of genes (scalp):
p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "BroadClust", 
    useGroups = broadOrder,
    pal = broadClustCmap,
    plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
    sizes = c(7, 0.2, 1.25, 2.5),
    geneSymbol = label_genes, 
    region = plotRegions, 
    loops = sub_plot_loop_list$scalp,
    tileSize=500,
    minCells=200
)

plotPDF(plotList = p, 
    name = "Super-Enhancer-Tracks-BroadClustScalpOnly_p2gLinks.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 6, height = 7)


# Subclustered Keratinocytes:
subgroup <- "Keratinocytes"
sub_dir <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/subclustered_%s", subgroup)
sub_proj <- loadArchRProject(sub_dir, force=TRUE)
kc_sub_cmap <- readRDS(paste0(scriptPath, sprintf("/atac_cmap_%s.rds", subgroup)))

kc_sub_order <- c(
  "aKc1", # "Basal.Kc_1",
  "aKc2", # "Spinous.Kc_2",
  "aKc3", # "Spinous.Kc_1",
  "aKc4", # "Infundibulum", # SOX9, DKK3
  "aKc5", # "Inf.Segment_1", # Lhx2, LGR5 high
  "aKc7", # "Inf.Segment_2", # Lhx2, LGR5 high
  "aKc6", # "Sebaceous", 
  "aKc8", # "Isthmus", # CD200 high
  "aKc9", # "Matrix", 
  "aKc10" # "Eccrine",
  #"aKc11" # "Unknown", (doublet?)
)

# Tracks of genes (scalp):
p <- plotBrowserTrack(
    ArchRProj = sub_proj, 
    groupBy = "FineClust", 
    useGroups = kc_sub_order,
    pal = kc_sub_cmap,
    plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
    sizes = c(7, 0.2, 1.25, 2.5),
    geneSymbol = label_genes, 
    region = plotRegions, 
    loops = sub_plot_loop_list$Keratinocytes,
    tileSize=500,
    minCells=200
)

plotPDF(plotList = p, 
    name = "Super-Enhancer-Tracks-FineClustKeratinocytesOnly_p2gLinks.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 6, height = 7)

##########################################################################################
# Violin plots of (integrated) RNA expression for select genes
##########################################################################################

# WARNING: this seems to fail if you have re-assigned the P2G links above. 
GImat <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix")
data_mat <- assays(GImat)[[1]]
rownames(data_mat) <- rowData(GImat)$name
sub_mat <- data_mat[label_genes,]

# These DO NOT match the order of the above matrix by default
grouping_data <- data.frame(cluster=factor(atac_proj$BroadClust, 
  ordered=TRUE, levels=atacBroadOrder))
rownames(grouping_data) <- getCellNames(atac_proj)
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

pList <- list()
for(gn in label_genes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$cluster %in% atacBroadOrder,]

  covarLabel <- "cluster"  

  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=gene, fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
    #+ geom_jitter(aes(group=Sample), size=0.025, 
    #  position=position_jitterdodge(seed=1, jitter.width=0.05, jitter.height=0.0, dodge.width=dodge_width))
    #+ stat_summary(fun="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
    # width=0.75, position=dodge,show.legend = FALSE)
    + scale_color_manual(values=broadClustCmap, limits=names(broadClustCmap), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=broadClustCmap)
    + guides(fill=guide_legend(title=covarLabel), 
      colour=guide_legend(override.aes = list(size=5)))
    + ggtitle(gn)
    + xlab("")
    + ylab("Integrated RNA Expression")
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Expression_Violin_byBroadClust.pdf"), width=10, height=4)
pList
dev.off()

# Keratinocytes only:

# These DO NOT match the order of the above matrix by default
kc_freqs <- getFreqs(sub_proj$FineClust)[kc_sub_order]
use_kc <- names(kc_freqs[kc_freqs > 200])
grouping_data <- data.frame(cluster=factor(sub_proj$FineClust, 
  ordered=TRUE, levels=kc_sub_order))
rownames(grouping_data) <- getCellNames(sub_proj)
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

pList <- list()
for(gn in label_genes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$cluster %in% use_kc,]

  covarLabel <- "cluster"  

  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=gene, fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
    #+ geom_jitter(aes(group=Sample), size=0.025, 
    #  position=position_jitterdodge(seed=1, jitter.width=0.05, jitter.height=0.0, dodge.width=dodge_width))
    #+ stat_summary(fun="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
    # width=0.75, position=dodge,show.legend = FALSE)
    + scale_color_manual(values=kc_sub_cmap, limits=names(kc_sub_cmap), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=kc_sub_cmap)
    + guides(fill=guide_legend(title=covarLabel), 
      colour=guide_legend(override.aes = list(size=5)))
    + ggtitle(gn)
    + xlab("")
    + ylab("Integrated RNA Expression")
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Expression_Violin_byKeratinocyteFineClust.pdf"), width=10, height=4)
pList
dev.off()

