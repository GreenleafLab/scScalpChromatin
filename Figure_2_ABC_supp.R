#!/usr/bin/env Rscript

##########################################################
# Analyses of full project peak to gene linkages
##########################################################

#Load ArchR (and associated libraries)
library(ArchR)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrastr)
library(InteractionSet)
library(GenomicInteractions)

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
plotDir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/p2gLink_plots"

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

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

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
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$peakName, "_", full_p2gGR$symbol))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF

##########################################################################################
# Load ABC model data to test for enrichment in peak-to-gene links
##########################################################################################

p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)

# Nasser et al 2021 ABC model predictions
abc_dt <- fread("/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/nasser2021_ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
keep_cols <- c(
  "chr", "start", "end", "name", "class", 
  "activity_base", "TargetGene", 
  #"TargetGeneTSS", "TargetGeneExpression", "TargetGenePromoterActivityQuantile", "TargetGeneIsExpressed", 
  "distance", 
  #"isSelfPromoter", "hic_contact", "powerlaw_contact", "powerlaw_contact_reference", "hic_contact_pl_scaled", 
  #"hic_pseudocount", "hic_contact_pl_scaled_adj", "ABC.Score.Numerator", "powerlaw.Score.Numerator", "powerlaw.Score", 
  "ABC.Score", 
  "CellType")
abc_dt <- abc_dt[,..keep_cols]

abc_gr <- makeGRangesFromDataFrame(abc_dt, keep.extra.columns=TRUE, ignore.strand=TRUE, 
    seqnames.field="chr", start.field="start", end.field="end")

chain <- rtracklayer::import.chain("/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/liftover/hg19ToHg38.over.chain")
abc_gr <- rtracklayer::liftOver(abc_gr, chain) %>% unlist()
abc_gr <- abc_gr[width(abc_gr) > 100]

# Rename gene column to use GInteraction constructor
colnames(mcols(abc_gr))[4] <- "symbol"

# Get GR of all promoter regions (to use for gene overlapping)
promoterGR <- promoters(getGenes(atac_proj))


convertP2GtoGInt <- function(p2gGR, promoterGR){
  # Create a genomic interaction object from a peak-to-gene genomic range
  # Uses a promoterGR as the non-peak anchor
  promoterGR <- promoterGR[!is.na(promoterGR$symbol)]
  names(promoterGR) <- promoterGR$symbol

  # Make sure we are only using genes for which we have promoter information
  p2gGR <- p2gGR[p2gGR$symbol %in% promoterGR$symbol]

  gr1 <- granges(p2gGR)
  gr2 <- granges(promoterGR[p2gGR$symbol])
  strand(gr2) <- NA # Need to remove strand information for swapAnchors to work correctly
  GInt <- GenomicInteractions(gr1, gr2)
  GInt <- swapAnchors(GInt, mode="order") %>% sort()
  mcols(GInt) <- mcols(p2gGR)
  GInt
}

sampleMatchDist <- function(target, background, size, 
  replace=FALSE, nbins=20, bin_type="size"){
  # Sample indices from target vector to match distribution of background vector
  # target and background must be numeric
  # target = numeric vector
  # background = numeric vector
  # size = how many indicees to return
  # replace = should indices be sampled with replacement
  # nbins = how many bins for selecting matching targets
  # bin_type = c("size" or "width"): should bins have equal numbers of elements (size),
  #   or equal linear numeric width (width)
  if(!replace){
    stopifnot(length(target) > size)
  }

  # Set target indices
  names(target) <- 1:length(target)

  # Make bins
  breaks <- makeBins(background, bins=nbins, bin.type=bin_type)$breaks
  breaks <- c(breaks, Inf)
  bin_names <- 1:length(breaks)
  names(breaks) <- bin_names

  # Approximate distribution by selecting desired frequency from each bin
  bg_bins <-sapply(background, function(x) min(which(breaks > x)))
  bg_bins <- factor(bg_bins, levels=bin_names)
  bg_bin_freqs <- getFreqs(bg_bins)
  bg_bin_freqs <- bg_bin_freqs[bin_names]
  bg_bin_probs <- bg_bin_freqs/sum(bg_bin_freqs)

  # Identify bins for target values
  tg_bins <- sapply(target, function(x) min(which(breaks > x)))
  tg_bins <- factor(tg_bins, levels=bin_names)

  # Take the appropriate number of indices from each bin to match distribution
  bin_target <- sapply(size*bg_bin_probs, ceiling)
  samp_idx <- lapply(names(bin_target), function(x){
    target_bin <- target[tg_bins == x]
    names(target_bin)[sample(1:length(target_bin), 
      size=min(bin_target[x], length(target_bin)), replace=replace)]
    }) %>% do.call(c,.)

  as.numeric(samp_idx)
}

###################################################################################################
# Plot overlap of P2G links in ABC links 
###################################################################################################

# All valid p2g links GInteraction
p2gGInt <- convertP2GtoGInt(getP2G_GR(atac_proj, corrCutoff=corrCutoff), promoterGR)

# All possible p2g links GInteraction
all_p2gGR <- getP2G_GR(atac_proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
all_p2gGInt <- convertP2GtoGInt(all_p2gGR, promoterGR)

# ABC model GInteraction
abc_GInt <- convertP2GtoGInt(abc_gr, promoterGR)

p2gGInt$pdist <- pairdist(p2gGInt, type="mid")
all_p2gGInt$pdist <- pairdist(all_p2gGInt, type="mid")

# Overlap for all valid p2g links
valid_p2g_olap <- length(p2gGInt[overlapsAny(p2gGInt, abc_GInt, type="any", ignore.strand=TRUE)])/length(p2gGInt)

# Overlap for all possible p2g links
all_p2g_olap <- length(all_p2gGInt[overlapsAny(all_p2gGInt, abc_GInt, type="any", ignore.strand=TRUE)])/length(all_p2gGInt)

########################################################################
# Permute distance-matched background and get overlap 
nboot <- 100
dm_p2g_olap <- mclapply(1:nboot, function(x){
  message(sprintf("Iteration %s...", x))
  set.seed(x)
  samp_p2gGInt <- all_p2gGInt[sampleMatchDist(all_p2gGInt$pdist, p2gGInt$pdist, size=length(p2gGInt))] %>% sort()
  length(samp_p2gGInt[overlapsAny(samp_p2gGInt, abc_GInt, type="any", ignore.strand=TRUE)])/length(samp_p2gGInt)
  }, mc.cores=5) %>% unlist()

# (Non-parallel)
# dm_p2g_olap <- lapply(1:nboot, function(x){
#   message(sprintf("Iteration %s...", x))
#   set.seed(x)
#   samp_p2gGInt <- all_p2gGInt[sampleMatchDist(all_p2gGInt$pdist, p2gGInt$pdist, size=length(p2gGInt))] %>% sort()
#   length(samp_p2gGInt[overlapsAny(samp_p2gGInt, abc_GInt, type="any", ignore.strand=TRUE)])/length(samp_p2gGInt)
#   }) %>% unlist()
########################################################################

dm_p2g_olap_mean <- mean(dm_p2g_olap)

# Overlap for subgroup p2g links
filt_full_p2gGR <- full_p2gGR[full_p2gGR$Correlation > corrCutoff & 
    full_p2gGR$VarQATAC > varCutoffATAC & 
    full_p2gGR$VarQRNA > varCutoffRNA]
filt_full_p2gGR$p2gID <- paste0(filt_full_p2gGR$peakName, "_", filt_full_p2gGR$symbol)

subsets <- unique(filt_full_p2gGR$source)

sub_p2g_olap <- lapply(subsets, function(ss){
  subP2G <- filt_full_p2gGR[filt_full_p2gGR$source == ss]
  subP2GInt <- convertP2GtoGInt(subP2G, promoterGR)
  length(subP2GInt[overlapsAny(subP2GInt, abc_GInt, type="any", ignore.strand=TRUE)])/length(subP2GInt)
  }) %>% unlist()
names(sub_p2g_olap) <- subsets

ol_res <- c(all_bg=all_p2g_olap, dist_matched=dm_p2g_olap_mean, all_p2g=valid_p2g_olap, sub_p2g_olap) %>% sort()
ol_df <- data.frame(group=names(ol_res), overlap=ol_res)
ol_df$group <- factor(ol_df$group, levels=names(ol_res), ordered=TRUE)

pdf(paste0(plotDir, "/fracOL_p2g_ABC_model.pdf"), width=8, height=6)
qcBarPlot(ol_df, cmap="royalblue1", barwidth=0.9, border_color=NA)
dev.off()

###################################################################################################
# Calculate p-value for enrichment of ABC model overlaps using one-sided Fisher's exact tests
###################################################################################################

# Construct peak-to-gene interaction obects to test for pair-wise overlaps and calculate distances, etc.
# (First, get all possible peak-to-gene comparisons to use as background)
all_p2gGR <- getP2G_GR(atac_proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
all_p2gGR$p2gID <- paste0(all_p2gGR$peakName, "_", all_p2gGR$symbol)
all_p2gGInt <- convertP2GtoGInt(all_p2gGR, promoterGR)
all_p2gGInt$pdist <- pairdist(all_p2gGInt, type="mid")

# Remove valid peak-to-gene linkages from background set (For Odds ratio / Fisher's exact test)
all_p2gGR <- all_p2gGR[all_p2gGR$p2gID %ni% filt_full_p2gGR$p2gID]

exp_subsets <- c("full_dataset", subsets)

olap_pvals <- sapply(exp_subsets, function(ss){
  message(sprintf("Calculating enrichment for subset %s...", ss))
  set.seed(1)
  if(ss == "full_dataset"){
    sub_P2G <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)
  }else{
    sub_P2G <- filt_full_p2gGR[filt_full_p2gGR$source == ss]
  }
  sub_P2GInt <- convertP2GtoGInt(sub_P2G, promoterGR)

  # Add pair-distances so we can select distance-matched background
  sub_P2GInt$pdist <- pairdist(sub_P2GInt, type="mid")

  # Sample distance-matched background peaks
  samp_p2gGInt <- all_p2gGInt[sampleMatchDist(all_p2gGInt$pdist, sub_P2GInt$pdist, size=length(sub_P2GInt))] %>% sort()

  olap <- sub_P2GInt[overlapsAny(sub_P2GInt, abc_GInt, 
    type="any", ignore.strand=TRUE)]

  unfilt_olap <- samp_p2gGInt[overlapsAny(samp_p2gGInt, abc_GInt, 
    type="any", ignore.strand=TRUE)]

  length(olap)/length(p2gGInt)
  length(unfilt_olap)/length(samp_p2gGInt)

  dnol <- length(olap) # Number of peak-to-gene links that overlap an ABC-model prediction
  dnnol <- length(sub_P2GInt) - dnol # Number of peak-to-gene links that do not overlap an ABC-model prediction
  nol <- length(unfilt_olap) # Number of non peak-to-gene links that overlap an ABC-model prediction
  nnol <- length(samp_p2gGInt) - nol # Number of non peak-to-gene links that do not overlap an ABC model prediction
  OR <- (dnol/dnnol)/(nol/nnol) # Odds ratio
  pval <- fisher.test(matrix(c(dnol, dnnol, nol, nnol),2,2), alternative="greater")$p.value
  c(OR=OR, pval=pval)
}) %>% t()

# full_dataset OR = 1.86, p-value ~ 0 (too low to be calculated)
# Scalp (non-subclustered) OR = 1.97, p-value ~ 0 (too low to be calculated)
# Keratinocytes OR = 1.59, p-value = 1.02e-105
# Fibroblasts OR = 1.81, p-value = 3.71e-92
# Endothelial OR = 2.18, p-value ~ 0 (too low to be calculated)
# T-lymphocytes OR = 2.30, p-value = 1.54e-126
# Myeloid OR = 1.92, p-value ~ 0 (too low to be calculated)

###################################################################################################