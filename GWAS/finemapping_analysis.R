#!/usr/bin/env Rscript

##########################################################################
# Analysis using finemapped SNPs
##########################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(plyranges)
  library(data.table)
  library(stringr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(parallel)
  library(ggrepel)
  library(ComplexHeatmap)
})

# Set Threads to be used
ncores <- 8
addArchRThreads(threads = ncores)


# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/fine_clustered"
fm_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gwas/PICS2"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

topN <- 80 # Number of genes to plot in heatmap

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_proj <- readRDS("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/preprocessing_output/scalp.rds")
plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

raw_finemapped_gr <- readRDS(paste0(fm_dir, "/unfiltered_finemapping_genomic_range.rds"))

# Some of the fine-mapped SNPs are duplicated (i.e. the Finacune SNPs sometimes have both FINEMAP and SuSiE finemapping results)
# Deduplicate trait-SNP pairs prior to proceeding with enrichment analyses:
raw_finemapped_gr <- raw_finemapped_gr[order(raw_finemapped_gr$fm_prob, decreasing=TRUE)]
raw_finemapped_gr$trait_snp <- paste0(raw_finemapped_gr$disease_trait, "_", raw_finemapped_gr$linked_SNP)
raw_finemapped_gr <- raw_finemapped_gr[!duplicated(raw_finemapped_gr$trait_snp)] %>% sort()

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

allGenesGR <- getGenes(atac_proj)

# P2G definition cutoffs
corrCutoff <- 0.5       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

##########################################################################################
# Plot finemapping probability by peak overlaps
##########################################################################################

prob_bins <- c(0, 0.01, 0.05, 0.1, 0.25, 1)

getPeakOLpct <- function(full_gr, peaks_gr, breaks=c(0, 0.05, 0.1, 0.25, 1)){
  # Calculate percent of SNPs overlapping peaks by bins of finemapped probability
  bin_ids <- makeBins(full_gr$fm_prob, breaks=breaks)$ids
  peak_ol_pct <- sapply(1:(length(breaks)-1), function(b){
    gr <- full_gr[bin_ids == b]
    nol <- overlapsAny(gr, peaks_gr) %>% sum()
    nol/length(gr)
  })
  bin_names <- sapply(2:length(breaks), function(i) paste0(breaks[i-1], "-", breaks[i]))
  names(peak_ol_pct) <- bin_names
  # Get number of SNPs in each bin
  freqs <- getFreqs(bin_ids)
  freqs <- freqs[order(as.integer(names(freqs)))]
  names(freqs) <- bin_names
  list(ol_pct=peak_ol_pct, ol_nSNPs=freqs)
}

# (We will put these in order for matching prior plots)
disease_traits <- list(
  skin_hair=c(
    "Alopecia areata",  # PICS AA SNPs (from GWC via PICS database)
    "Psoriasis",
    "Eczema",
    "Cutaneous squamous cell carcinoma",
    "Cutaneous malignant melanoma",
    "Vitiligo",
    "Male-pattern baldness", # PICS via GWC
    "Balding_Type4", # Finacune finemapped
    "Sunburns",
    "Hair color"
  ),
  other_autoimmune=c(
    "Asthma",
    "Crohn's disease",
    "Systemic lupus erythematosus",
    "Celiac disease",
    "Ulcerative colitis"
  ),
  brain_traits=c(
    "Parkinson's disease",
    "Schizophrenia",
    "Major depressive disorder",
    "Educational attainment (years of education)",
    "Neuroticism"
  ),
  other_traits=c(
    "Red blood cell count",
    "BMI", # Finacune finemapped
    "Height",
    "SBP" # Finacune finemapped
  )
)
all_disease_traits <- do.call(c,disease_traits) %>% unname()

pList <- lapply(all_disease_traits, function(dt){
  message(sprintf("Testing trait %s...", dt))
  ol <- getPeakOLpct(raw_finemapped_gr[raw_finemapped_gr$disease_trait == dt], peaks_gr=allPeaksGR, breaks=prob_bins)
  pct_ol <- ol$ol_pct
  nSNPs <- ol$ol_nSNPs
  xlabs <- paste0(names(nSNPs), sprintf("\n(%s)", nSNPs))
  df <- data.frame(fm_probs=names(pct_ol), percent_snps_overlapping=pct_ol)
  df$fm_probs <- factor(df$fm_probs, levels=names(pct_ol), ordered=TRUE)
  qcBarPlot(df, cmap="grey80", barwidth=0.9, border_color="black") + 
    scale_y_continuous(limits=c(0, 0.4), expand = c(0, 0)) + 
    ggtitle(dt) + scale_x_discrete(labels=xlabs)
  })

pdf(paste0(plotDir, "/pctPeaksOL_by_fm_probs_eachTrait.pdf"), width=5, height=5)
pList
dev.off()

# Plot box plot showing increase in fraction of snps overlapping peaks in large groups of SNP categories
use_groups <- c("skin_hair", "other_autoimmune", "brain_traits", "other_traits")
sub_traits_list <- disease_traits[use_groups]

ol_df <- lapply(names(sub_traits_list), function(x){
  message(sprintf("Calculating overlaps for group %s...", x))
  dts <- sub_traits_list[[x]]
  df <- lapply(dts, function(dt){
      ol <- getPeakOLpct(raw_finemapped_gr[raw_finemapped_gr$disease_trait == dt], peaks=allPeaksGR, breaks=prob_bins)
      pct_ol <- ol$ol_pct
      nSNPs <- ol$ol_nSNPs
      df <- data.frame(fm_probs=names(pct_ol), nSNPs=nSNPs, pct_snps_ol=pct_ol)
      df$fm_probs <- factor(df$fm_probs, levels=names(pct_ol), ordered=TRUE)
      df$trait <- dt
      df
    }) %>% do.call(rbind,.)
  df$trait_set <- x
  df
  }) %>% do.call(rbind,.)

ol_df$trait_set <- factor(ol_df$trait_set, ordered=TRUE, levels=use_groups)

# Plot multiple category box plot
cmap <- cmaps_BOR$stallion
dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

p <- (
  ggplot(ol_df, aes(x=fm_probs, y=pct_snps_ol, color=trait_set, fill=trait_set))
  + geom_boxplot(alpha=0.5, outlier.shape = NA) # Hide fliers (we show them with geom_jitter)
  + geom_jitter(aes(group=trait_set), size=1.0, color="black",
     position=position_jitterdodge(seed=1, jitter.width=0.25, jitter.height=0.0, dodge.width=dodge_width))
  + scale_fill_manual(values=cmap)
  + scale_color_manual(values=cmap)
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          #aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
          axis.text.x = element_text(angle = 90, hjust = 1)) 
  + scale_y_continuous(limits=c(0, 0.39), expand = c(0, 0))
)
pdf(paste0(plotDir, "/pctPeaksOL_by_fm_probs.pdf"), width=10, height=4)
p
dev.off()


# Read in cell type specific peaks
cts_peak_files <- list.files(
  path="/oak/stanford/groups/wjg/boberrey/hairATAC/results/GWAS/ldsc/FineClust_specific_peaks",
  pattern="*.bed$",
  full.names=TRUE
)
names(cts_peak_files) <- basename(cts_peak_files) %>% sapply(., function(x) gsub("_specific_peaks.bed", "", x)) %>% unname()

cts_peaks <- lapply(cts_peak_files, function(pf){
  dt <- fread(pf)
  GRanges(seqnames=dt$V1, IRanges(dt$V2+1, end=dt$V3+1))
  })
names(cts_peaks) <- names(cts_peak_files)

fisherTestSNPs <- function(peaks_gr, snps_gr, disease_trait){
  # Calculate fisher enrichment for a specific trait in a given peak set
  ######################################################################
  # peaks_gr = peakset to use for testing enrichment
  # snps_gr = full fine-mapping GR
  # disease_trait = which trait to test enrichment for

  # SNPs overlapping peakset
  nondis_gr <- snps_gr[snps_gr$disease_trait != disease_trait]
  nondis_gr <- nondis_gr[!duplicated(nondis_gr$linked_SNP)] # remove duplicated SNPs
  nol <- overlapsAny(nondis_gr, peaks_gr) %>% sum() # n non-disease SNPs overlapping
  nnol <- length(nondis_gr) - nol # n non-disease SNPs not overlapping

  # Trait SNPs overlapping peakset
  dis_gr <- snps_gr[snps_gr$disease_trait == disease_trait]
  dis_gr <- dis_gr[!duplicated(dis_gr$linked_SNP)] # remove duplicated SNPs (shouldn't be any)
  dnol <- overlapsAny(dis_gr, peaks_gr) %>% sum() # n disease SNPs overlapping
  dnnol <- length(dis_gr) - dnol # n disease SNPs not overlapping

  OR <- (dnol/dnnol)/(nol/nnol)
  pval <- fisher.test(matrix(c(dnol, dnnol, nol, nnol),2,2), alternative="greater")$p.value

  list(trait=disease_trait, ol_dis_snps=dnol, nol_dis_snps=dnnol, ol_snps=nol, nol_snps=nnol, OR=OR, fisher_pval=pval)
}

# (Both PICS and finacune finemapped SNPs have ~20% with prob > 0.01)
filt_snp_gr <- raw_finemapped_gr[raw_finemapped_gr$fm_prob >= 0.01]

res_df <- lapply(names(cts_peaks), function(ct){
    message(sprintf("Testing %s peaks...", ct))
    peaks_gr <- cts_peaks[[ct]]
    dis_res <- lapply(all_disease_traits, function(dt){
      fisherTestSNPs(peaks_gr, filt_snp_gr, dt)
    })
    dis_res <- do.call(rbind, lapply(dis_res, data.frame))
    dis_res$cluster <- ct
    dis_res
  }) %>% do.call(rbind,.)

res_df$padj <- p.adjust(res_df$fisher_pval, method="fdr")
res_df <- res_df[order(res_df$padj, decreasing=FALSE),]
res_df$mlog10padj <- -log10(res_df$padj)

# Plot Dot Plot of fischer enrichment
# Specify order of clusters (Fine Clust)
FatacOrder <- c(
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
LFatacOrder <- unlist(atac.FineClust)[FatacOrder]
FatacOrder <- FatacOrder[FatacOrder %in% unique(res_df$cluster)]

# Plot separate dot plot for each group:
max_mlogpval <- 15
colorLims <- c(0, max_mlogpval)
sizeLims <- c(min(res_df$OR), max(res_df$OR))

pList <- list()
for(grp in names(disease_traits)){

  sub_res_df <- res_df[res_df$trait %in% disease_traits[[grp]],]
  # Determine cluster and gene order:
  wide_df <- unmelt(sub_res_df, row_col="trait", col_col="cluster", val_col="OR")
  wide_df <- prettyOrderMat(wide_df[,FatacOrder], clusterCols=FALSE)$mat

  # Prepare for plotting
  plot_df <- sub_res_df[sub_res_df$cluster %in% FatacOrder,]
  plot_df$grp <- unlist(atac.FineClust)[as.character(plot_df$cluster)]
  plot_df$mlog10padj <- ifelse(plot_df$mlog10padj > max_mlogpval, max_mlogpval, plot_df$mlog10padj)

  grp_order <- colnames(wide_df)
  #trait_order <- rownames(wide_df) %>% rev()
  trait_order <- disease_traits[[grp]] %>% rev()

  pList[[grp]] <- dotPlot(plot_df, xcol="grp", ycol="trait", color_col="mlog10padj", size_col="OR", 
    xorder=unlist(atac.FineClust)[grp_order], yorder=trait_order, 
    cmap=cmaps_BOR$wolfgang_extra, aspectRatio=nrow(wide_df)/ncol(wide_df), 
    sizeLims=sizeLims, colorLims=colorLims)
}

pdf(paste0(plotDir, "/fmSNP_enrichment_fisher_clusters_dotPlots.pdf"), width=12, height=8)
pList
dev.off()

# Save results table
save_df <- res_df
save_df$Lcluster <- unlist(atac.FineClust)[as.character(save_df$cluster)]
save_df <- save_df[!is.na(save_df$Lcluster),]
write.table(save_df, file=paste0(wd, "/fmSNP_enrichment_fisher_clusters_results.tsv"), 
  quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

rm(raw_finemapped_gr); gc() 

##########################################################################################
# Links Fine-mapped SNPs to candidate genes using Peak-to-Gene links
##########################################################################################

finemapped_GR <- readRDS(paste0(fm_dir, "/filtered_finemapping_genomic_range.rds"))
# Some of the fine-mapped SNPs are duplicated (i.e. the Finacune SNPs sometimes have both FINEMAP and SuSiE finemapping results)
# Deduplicate trait-SNP pairs prior to proceeding with enrichment analyses:
finemapped_GR <- finemapped_GR[order(finemapped_GR$fm_prob, decreasing=TRUE)]
finemapped_GR$trait_snp <- paste0(finemapped_GR$disease_trait, "_", finemapped_GR$linked_SNP)
finemapped_GR <- finemapped_GR[!duplicated(finemapped_GR$trait_snp)] %>% sort()

# Load full project p2g links, plot loops, etc.
full_p2gGR <- readRDS(file=paste0(wd, "/multilevel_p2gGR.rds")) # NOT merged or correlation filtered
full_coaccessibility <- readRDS(file=paste0(wd, "/multilevel_coaccessibility.rds"))
plot_loop_list <- readRDS(file=paste0(wd, "/multilevel_plot_loops.rds"))

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(atac_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$symbol, "-", full_p2gGR$peakName))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF

##########################################################################################
# Identify Finemapped SNPs linked to genes
##########################################################################################

p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff, 
  varCutoffATAC=varCutoffATAC, varCutoffRNA=varCutoffRNA, filtNA=TRUE)

pol <- findOverlaps(p2gGR, finemapped_GR, type="any", maxgap=-1L, ignore.strand=TRUE)
expandFmGR <- finemapped_GR[to(pol)]
expandFmGR$linkedGene <- p2gGR[from(pol)]$symbol
expandFmGR$linkedPeak <- p2gGR[from(pol)]$peakName
expandFmGR$p2gCorr <- p2gGR[from(pol)]$Correlation
expandFmGR$SNP_to_gene <- paste(expandFmGR$linked_SNP, expandFmGR$linkedGene, sep="_")

# Heatmaps of linked genes and expression by cell type
source(paste0(scriptPath, "/cluster_labels.R"))
rna_proj$LFineClust <- unlist(rna.FineClust)[rna_proj$FineClust]
exclude_clust <- c("Unknown", "Cyc.Tc", "Plasma_contam", "TCR.macs", "McSC")
rna_proj <- rna_proj[,rna_proj$LFineClust %ni% exclude_clust]

countMat <- GetAssayData(object=rna_proj, slot="counts")
groupedCountMat <- averageExpr(countMat, rna_proj$FineClust) # Calculates average log2cp10k, so must only subset afterwards

# Specify order of clusters (Fine Clust)
rnaOrder <- c(
    # Lymphoid / T-cells
    "rTc3", # "Tc",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG 
    "rTc5", # "NK", # Natural Killer cells (XCL1, XCL2, GNLY, NKG7, KLRD1) 
    "rTc1", # "Th_1", # T-helper
    "rTc2", # "Th_2", # T-helper
    "rTc4", # "Treg", # Regulatory T cells: IKZF2, IL2RA, CTLA4, FOXP3 
    "rTc6", # "Cyc.Tc", # Cycling markers (MKI67, TOP2A, etc.)
    # B/Plasma
    "rBc1", # "Plasma",
    # Myeloid
    "rMy5", # "M1.macs", # IL15, IL32, CCR7 (CCL19, CCL17)
    "rMy6", # "TCR.macs", # CD3, TCR gene positive macrophages?
    "rMy2", # "M2.macs_1", # C1Qa/b/c, FOLR2, CD14, CD163, (CCL13?)
    "rMy3", # "M2.macs_2", # CXCL2, CXCL3, (CCL20?, S100A8/9?) 
    "rMy7", # "TREM2.macs", # TREM2
    "rMy1", # "cDC2", # CD1c, CLEC10a (conventional DCs - type 2)
    "rMy4", # "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1 
    # Keratinocytes
    "rKc1", # "Basal.Kc_1",
    "rKc2", # "Spinous.Kc_1",
    "rKc3", # "Spinous.Kc_2", 
    "rKc4", # "Spinous.Kc_3",
    "rKc6", # "Cyc.Kc_1", # Cycling keratinocytes
    "rKc5", # "Infundibulum", # SOX9, DKK3
    "rKc7", # "Inf.Segment", # Lhx2, LGR5 high
    "rKc8", # "Sebaceous", 
    "rKc9", # "Isthmus", # CD200 high
    "rKc10", # "Eccrine",
    # Fibroblasts
    "rFb1", # "Fb_1", # CXCL1,2,3
    "rFb3", # "Fb_2", # CCL19, CXCL12
    "rFb4", # "Fb_3", # APCDD1, COL18A1, F13A1
    "rFb5", # "Fb_4", # WISP2, AOX1, ARFGEF3
    "rFb6", # "Fb_5", # NECAB1, SCN7A
    "rFb2", # "D.Sheath", # COL11A1, EDNRA
    "rFb7", # "D.Papilla", # Many markers HHIP, PTCH1, etc.
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
    "rMe1" # "Melanocytes"
)

rnaOrder <- rnaOrder[rnaOrder %in% colnames(groupedCountMat)]

# Get colors for cluster annotation
# (Link FineClusts to BroadClust cmap)
cM <- as.matrix(confusionMatrix(rna_proj$LFineClust, rna_proj$BroadClust))
map_colors <- apply(cM, 1, function(x)colnames(cM)[which.max(x)])
plotColors <- map_colors[unlist(rna.FineClust)[rnaOrder]]


plotFMsnpsToGenes <- function(expandFmGR, groupedCountMat, trait, cmap, plotDir, backgroundGenes, clusterOrder=NULL, topN=80){
  # Plot finemapped SNPs to genes through peak to gene linkages
  #############################################################
  # expandFmGR = 'expanded' GR of finemapped SNPs linked to genes
  # groupedCountMat = matrix of expression grouped by cell types
  # trait = the GWAS trait/disease to investigate
  # cmap = color map of clusters in groupedCountMat
  # plotDir = plot directory path
  # backgroundGenes = list of genes to use as background for GO term enrichment
  # clusterOrder = order of clusters for heatmap (if null, will bicluster heatmap)
  # topN = number of top genes (by cumulative finemapping probability)
  trait_name <- strsplit(trait, split=" ")[[1]] %>% paste(., collapse="_")
  trait_gr <- expandFmGR[expandFmGR$disease_trait %in% c(trait)]
  trait_gr <- trait_gr[order(trait_gr$fm_prob, decreasing=TRUE)]
  # Further filter by finemapping posterior probability
  trait_gr <- trait_gr[trait_gr$fm_prob >= 0.01]
  trait_gr <- trait_gr[!is.na(trait_gr$linkedGene)]
  trait_gr <- trait_gr[!duplicated(trait_gr$SNP_to_gene)]
  # Save fmGWAS-linked genes
  message(sprintf("Saving fmGWAS-genes for traint %s...", trait))
  saveRDS(trait_gr, paste0(plotDir, sprintf("/%s_fmGWAS_genes.rds", trait_name)))
  # Cumulative fine-mapping posterior probability per gene
  totalFMprobs <- trait_gr %>% as.data.frame() %>% group_by(linkedGene) %>% summarize(total_FM_prob=sum(fm_prob)) %>% as.data.frame()
  namedFMprobs <- totalFMprobs$total_FM_prob
  names(namedFMprobs) <- totalFMprobs$linkedGene
  message(sprintf("Trait %s has a total of %s unique fmGWAS-linked genes", trait, length(namedFMprobs)))
  # Keep only top N genes for plotting:
  plotGenes <- totalFMprobs[order(totalFMprobs$total_FM_prob, decreasing=TRUE),]$linkedGene %>% head(topN)
  # GO term enrichments
  traitGO <- rbind(
    calcTopGo(backgroundGenes, 
      interestingGenes=plotGenes,
      ontology="MF"),
    calcTopGo(backgroundGenes, 
      interestingGenes=plotGenes,
      ontology="BP")
    )
  traitGO <- traitGO[order(as.numeric(traitGO$pvalue), decreasing=FALSE),]
  pdf(paste0(plotDir, sprintf("/%s_GO_term_enrichments_MFBP_top%s.pdf", trait_name, topN)), width=8, height=6)
  print(topGObarPlot(traitGO, cmap = cmaps_BOR$comet, 
          nterms=6, border_color="black", 
          barwidth=0.85, title=sprintf("%s GWAS linked Genes", trait)))
  dev.off()
  avgMat <- groupedCountMat[plotGenes,] 
  avgMat <- t(scale(t(avgMat)))
  # Cluster for heatmap
  if(is.null(clusterOrder)){
    plotMat <- prettyOrderMat(avgMat[plotGenes,], clusterCols=TRUE, cutOff=1)$mat
  }else{
    plotMat <- prettyOrderMat(avgMat[plotGenes,clusterOrder], clusterCols=FALSE, cutOff=1)$mat
  }
  plotMat[plotMat > 3] <- 3
  plotMat[plotMat < -3] <- -3
  colnames(plotMat) <- unlist(rna.FineClust)[colnames(plotMat)]
  # Barplot for number of linked peaks per gene
  linkedGeneFreq <- getFreqs(trait_gr$linkedGene)
  pdf(paste0(plotDir, sprintf("/%s_GWASlinkedGenesHeatmap.pdf", trait_name)), width=15, height=12)
  fontsize <- 6
  ht_opt$simple_anno_size <- unit(0.25, "cm")
  ta <- HeatmapAnnotation(rna_cluster=plotColors,col=list(rna_cluster=broadClustCmap), 
    show_legend=c(rna_cluster=FALSE), show_annotation_name=c(rna_cluster=FALSE))
  ra <- HeatmapAnnotation(
    nPeaks=anno_barplot(
      linkedGeneFreq[rownames(plotMat)], 
      bar_width=1, height=unit(3, "cm"),
      gp=gpar(fill="grey", fontsize=fontsize),
      labels_gp=gpar(fontsize=fontsize)
    ), 
    FMP=anno_barplot(
      namedFMprobs[rownames(plotMat)],
      bar_width=1, height=unit(3, "cm"),
      gp=gpar(fill="blue", fontsize=fontsize), ylim=c(0, max(1, max(namedFMprobs)))
    ),
    which="row")
  hm <- BORHeatmap(
    plotMat, 
    limits=c(-2.0,2.0), 
    clusterCols=FALSE, clusterRows=FALSE,
    labelCols=TRUE, labelRows=TRUE,
    dataColors = cmaps_BOR$sunrise,
    top_annotation = ta,
    right_annotation = ra,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = fontsize),
    column_names_gp = gpar(fontsize = fontsize),
    width = ncol(plotMat)*unit(0.4, "cm"),
    height = nrow(plotMat)*unit(0.25, "cm"),
    legendTitle="row Z-score",
    border_gp = gpar(col="black") # Add a black border to entire heatmap
    )
  draw(hm)
  dev.off()
}


# Plot selected traits
plotFMsnpsToGenes(expandFmGR, groupedCountMat, "Alopecia areata", cmap=plotColors, plotDir=plotDir, 
  backgroundGenes=unique(p2gGR$symbol), clusterOrder=rnaOrder, topN=topN)

plotFMsnpsToGenes(expandFmGR, groupedCountMat, "Balding_Type4", cmap=plotColors, plotDir=plotDir, 
  backgroundGenes=unique(p2gGR$symbol), clusterOrder=rnaOrder, topN=topN)

plotFMsnpsToGenes(expandFmGR, groupedCountMat, "Eczema", cmap=plotColors, plotDir=plotDir, 
  backgroundGenes=unique(p2gGR$symbol), clusterOrder=rnaOrder, topN=topN)

plotFMsnpsToGenes(expandFmGR, groupedCountMat, "Hair color", cmap=plotColors, plotDir=plotDir, 
  backgroundGenes=unique(p2gGR$symbol), clusterOrder=rnaOrder, topN=topN)


############################################################################################################

# Enrichment of TFs for AGA

trait_gr <- expandFmGR[expandFmGR$disease_trait %in% c("Balding_Type4")]
trait_gr <- trait_gr[order(trait_gr$fm_prob, decreasing=TRUE)]
# Further filter by fine-mapping posterior probability
trait_gr <- trait_gr[trait_gr$fm_prob >= 0.01]
trait_gr <- trait_gr[!is.na(trait_gr$linkedGene)]
trait_gr <- trait_gr[!duplicated(trait_gr$SNP_to_gene)]

# Cumulative fine-mapping posterior probability per gene
totalFMprobs <- trait_gr %>% as.data.frame() %>% group_by(linkedGene) %>% summarize(total_FM_prob=sum(fm_prob)) %>% as.data.frame()
namedFMprobs <- totalFMprobs$total_FM_prob
names(namedFMprobs) <- totalFMprobs$linkedGene

# Keep only top N genes for plotting:
plotGenes <- totalFMprobs[order(totalFMprobs$total_FM_prob, decreasing=TRUE),]$linkedGene %>% head(topN)

trait_genes <- trait_gr$linkedGene %>% unique()
bg_genes <- unique(p2gGR$symbol)
TFdb <- fread("/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gene_sets/lambert_2018_TF_master_list.csv")
TFs <- TFdb %>% filter(Is_TF == "Yes") %>% pull(Symbol)
valid_TFs <- TFs[TFs %in% bg_genes]

# Hypergeometric enrichment
q <- sum(trait_genes %in% valid_TFs)  # q = number of white balls drawn without replacement
m <- length(valid_TFs)                # m = number of white balls in urn
n <- length(bg_genes) -  m            # n = number of black balls in urn
k <- length(trait_genes)              # k = number of balls drawn from urn
TFphypPval <- phyper(q, m, n, k, lower.tail=FALSE, log.p=FALSE)

enrichment <- (q/length(trait_genes))/(m/length(bg_genes))
OR <- (q/(k-q))/(m/n)
