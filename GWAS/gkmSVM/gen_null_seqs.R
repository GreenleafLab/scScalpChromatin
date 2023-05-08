#!/usr/bin/env Rscript

##########################################################################
# Generate fasta files for training gkmSVM models
##########################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  library(parallel)
  library(ComplexHeatmap)
})

# Set Threads to be used
ncores <- 8
addArchRThreads(threads = ncores)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/plotting_config.R"))

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/fine_clustered"
outdir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/GWAS/gkmSVM"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Set/Create Working Directory to Folder
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

# Set genome
genome <- BSgenome.Hsapiens.UCSC.hg38.masked

# Load ArchR project
atac_proj <- loadArchRProject(wd, force=TRUE)

# Peak width to use
peak_width <- 1001

# Min cells per cell type to use
min_cells <- 200

# Number of peaks to use per cell type
use_peaks <- 75000
nmarker_peaks <- 7500

##########################################################################################
# Preparing Data
##########################################################################################

##### Functions #######

getRandomPos <- function(n, genome, use_chr=NULL, width=1, 
  blacklist_gr=NULL, non_overlapping=TRUE){
  # Generate random genomic positions given provided BSgenome object
  ###################################################################
  # n = number of sequences to return
  # genome = BSgenome object to get sequences from
  # use_chr = vector of valid chromosomes to use (will use all if not provided)
  # width = sequence width
  # blacklist_gr = Genomic Range of blacklisted regions (e.g. Repeat Masks, etc.)
  # non_overlapping = return only non-overlapping regions

  # Get chromosome lengths
  if(is.null(use_chr)){
    use_chr <- seqnames(genome)
  }
  chr_lengths <- seqlengths(genome)[use_chr]
  nseqs <- 0

  full_gr <- GRanges()
  while(nseqs < n){
    nneed <- n - nseqs
    nsample <- min(nneed*5, n)
    message(sprintf("Still need %s more sequences...", nneed))
    # Sample random chromosomes (with respect to chr length)
    random_chr <- sample(x=use_chr, size=nsample, prob=chr_lengths, replace=TRUE)

    # Sample random positions
    message(sprintf("Sampling %s sequences...", nsample))
    random_pos <- sapply(random_chr, function(chr){
      sample(chr_lengths[chr]-width, size=1)
      })
    gr <- GRanges(random_chr, IRanges(random_pos, random_pos+width), seqinfo=seqinfo(genome)[use_chr])
    if(!is.null(blacklist_gr)){
      message("Filtering sequences by provided blacklist...")
      gr <- gr[!overlapsAny(gr, blacklist_gr, ignore.strand=TRUE)]
    }
    full_gr <- c(full_gr, gr) %>% trim_oob() %>% trim_N_seqs(.,genome)
    if(non_overlapping){
      full_gr <- resize(reduce(full_gr), width=width, fix="center")
    }
    nseqs <- length(full_gr)
  }
  message("Sorting GR...")
  full_gr <- sort(full_gr[sample(1:length(full_gr), size=n, replace=FALSE)])
  message("Done.")
  return(full_gr)
}

selectTargetSeqs <- function(peaks_gr, targets_gr, blacklist_gr, nseqs,
  nbins=20, bin_type="size", seed=123){
  # Identify random genomic sequences to be used for model training
  # peaks_gr = true sequences for cell type
  # targets_gr = target null sequences to sample from
  # nseqs = vector of number of sequences to take from each list
 
  # Add peaks to blacklist
  message("Preparing targets...")
  blacklist_gr <- sort(reduce(c(blacklist_gr, granges(peaks_gr))))

  # Filter targets by blacklist
  targets_gr <- targets_gr[!targets_gr %over% blacklist_gr]

  # Prepare a 'bin grid' of parameters to match distribution on
  message("Preparing bins by GC content for sampling random sequences...")
  GC_breaks <- makeBins(peaks_gr$GC, bins=nbins, bin.type=bin_type)$breaks
  GC_breaks <- c(GC_breaks, Inf)
  names(GC_breaks) <- 1:length(GC_breaks)

  # Identify all possible bins, and then estimate distribution to sample from
  bin_names <- names(GC_breaks)

  # Approximate true distribution by binning
  peak_bins <- sapply(peaks_gr$GC, function(x) min(which(GC_breaks > x)))
  peak_bins <- factor(peak_bins, levels=bin_names)
  peak_bin_freqs <- getFreqs(peak_bins)
  peak_bin_freqs <- peak_bin_freqs[bin_names]
  peak_bin_probs <- peak_bin_freqs/sum(peak_bin_freqs)

  # Identify bins on target regions
  GC_target_bins <- sapply(targets_gr$GC, function(x) min(which(GC_breaks > x)))
  GC_target_bins <- factor(GC_target_bins, levels=bin_names)
  
  # Take the appropriate number of sequences from each bin to match original distribution
  set.seed(seed)
  bin_target <- sapply(nseqs*peak_bin_probs, ceiling)
  sub_target_gr <- lapply(names(bin_target), function(x){
    bin_gr <- targets_gr[GC_target_bins==x]
    samp_gr <- bin_gr[sample(1:length(bin_gr), size=min(bin_target[x], length(bin_gr)), replace=FALSE)]
    samp_gr
  }) %>% as(., "GRangesList") %>% unlist() %>% sort()

  return(sub_target_gr)
}

##########################################################################################
# Prepare 'True' cluster peaks
##########################################################################################

# Get all peaks from ArchR project
peaks_gr <- getPeakSet(atac_proj)
peaks_gr$peakName <- (peaks_gr %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
seqinfo(peaks_gr) <- seqinfo(genome)[seqlevels(peaks_gr)]

# Identify blacklist regions in genome
use_chr <- unique(seqnames(peaks_gr))
use_masks <- c("AGAPS", "AMB")

# Obtain full mask GR
masked_gr <- lapply(use_chr, function(chr){
  message(sprintf("Getting masks from %s...", chr))
  masks <- Biostrings::masks(genome[[as.character(chr)]])
  valid_masks <- masks@NAMES[sapply(masks@nir_list, length)>0]
  valid_masks <- valid_masks[valid_masks %in% use_masks]
  full_range <- sapply(valid_masks, function(msk){GRanges(chr, as(masks[[msk]], "IRanges"))}) %>%
    as(., "GRangesList") %>% unlist()
  full_range$mask <- names(full_range)
  full_range
  }) %>% as(., "GRangesList") %>% unlist() %>% sort()

# Now, identify peaks per cell type and write fasta files (separated by chr for cross validation
# purposes.)
cell_types <- getFreqs(atac_proj$FineClust)
cell_types <- names(cell_types)[cell_types > min_cells]

# Identify marker peaks for each retained cluster
# We will use the 'top' N peaks (by score) for each cluster as well as all
# 'marker peaks' for that cluster to enhance cell-type specificity of model training data

# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "FineClust",
    useGroups = cell_types,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markerPeaks, cutOff="FDR <= 0.1 & Log2FC >= 0.5")
markerList <- lapply(markerList, function(mdf){
  mdf$peakName <- paste0(mdf$seqnames, "_", mdf$start, "_", mdf$end)
  head(mdf, nmarker_peaks) # Subset to most significant N marker peaks
  })

# Identify 'top' peaks for cell type and combine with appropriate set of marker peaks
cts_peaks <- lapply(cell_types, function(ct){
  message(sprintf("Gettting peaks from cluster %s...", ct))
  peaks <- getClusterPeaks(atac_proj, clusterNames=c(ct), peakGR=peaks_gr, originalScore=TRUE)
  # Sort peaks by score
  peaks <- peaks[order(peaks$score, decreasing=TRUE)]
  # Prioritize 'marker' peaks for each cell type
  marker_peak_names <- markerList[[ct]]$peakName
  peaks <- c(peaks[peaks$peakName %in% marker_peak_names], peaks)
  # Deduplicate and clean up
  peaks <- peaks[!duplicated(peaks$peakName)]
  peaks <- resize(peaks,width=peak_width, fix="center") %>% trim_oob() %>% trim_N_seqs(.,genome)
  peaks <- peaks %>% head(use_peaks) %>% sort()
  peaks$GC <- gcContent(peaks, genome)
  peaks
  })
names(cts_peaks) <- cell_types

# Only train models for cell types that have sufficient peaks for model training
cts_peaks <- cts_peaks[unlist(lapply(cts_peaks, length)) > use_peaks*0.75]

##########################################################################################
# Show similarity between input training data (Jaccard index)
##########################################################################################

cell_types <- names(cts_peaks)

# Calculate jaccard of peaks used in each model
jaccard <- function(v1, v2){
  # Calculate jaccard index between two vectors v1 and v2
  length(intersect(v1,v2))/length(unique(c(v1,v2)))
}

jacc_mat <- matrix(0, length(cell_types), length(cell_types))
for(i in 1:nrow(jacc_mat)){
  for(j in 1:ncol(jacc_mat)){
    v1 <- cts_peaks[[i]]$peakName
    v2 <- cts_peaks[[j]]$peakName
    jacc_mat[i,j] <- jaccard(v1, v2)
  }
}
rownames(jacc_mat) <- cell_types
colnames(jacc_mat) <- cell_types

jacc_mat[jacc_mat == 1.0] <- NA

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

pmat <- jacc_mat

# Add cluster labels
source(paste0(scriptPath, "/cluster_labels.R"))
LatacOrder <- unlist(atac.FineClust)[atacOrder]
rownames(pmat) <- unlist(atac.FineClust)[rownames(pmat)]
colnames(pmat) <- unlist(atac.FineClust)[colnames(pmat)]
LatacOrder <- LatacOrder[LatacOrder %in% colnames(pmat)]
pmat <- pmat[LatacOrder, LatacOrder]

# Jaccard heatmap
pdf(paste0(outdir, "/TrainingPeaks_jaccard_index_HM_1000bp.pdf"), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  pmat,
  limits=c(0,0.75),
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  column_title="FineClust", 
  column_title_side="bottom",
  row_title="FineClust",
  dataColors = cmaps_BOR$solar_extra,
  row_names_side = "left",
  width = ncol(pmat)*unit(0.5, "cm"),
  height = nrow(pmat)*unit(0.5, "cm"),
  legendTitle="Jaccard Index",
  border_gp=gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()

##########################################################################################
# Create 'null' training data
##########################################################################################

# Create full ATAC peak-set from multiple tissue types (can be overlapping, but non-identical)
all_peaks <- peaks_gr %>% resize(peak_width, fix="center") %>% trim_oob() %>% trim_N_seqs(.,genome)
#all_peaks <- c(all_peaks, brain_peaks, mpal_peaks, tcga_peaks) %>% granges() %>% sort() %>% unique()

# Create random peak set
set.seed(123)
random_peaks <- getRandomPos(4e6, genome=genome, 
    use_chr=seqlevels(all_peaks), width=peak_width, 
    blacklist_gr=c(masked_gr, all_peaks), non_overlapping=FALSE)

# Add GC content to peak sets
#all_peaks$GC <- gcContent(all_peaks, genome)
random_peaks$GC <- gcContent(random_peaks, genome)

# Pre-filter random_peaks to better match target GC distribution
all_cts_peaks <- as(cts_peaks, "GRangesList") %>% unlist()
random_peaks <- selectTargetSeqs(all_cts_peaks, targets_gr=random_peaks, 
  blacklist_gr=c(masked_gr, all_peaks), nseqs=2e6)

#####################################################################################
# Find null seqs for each cell type in parallel

#peak_group <- invertList(labelHierarchy)

null_regions <- mclapply(names(cts_peaks), function(x){

  set.seed(123)
  nseqs <- length(cts_peaks[[x]])
  full_null <- selectTargetSeqs(cts_peaks[[x]], targets_gr=random_peaks, blacklist_gr=masked_gr, nseqs=nseqs)
  message(sprintf("Successfully found %s seqs for %s", length(full_null), x))
  print(t.test(cts_peaks[[x]]$GC, full_null$GC)) # usually within 1% GC content
  full_null
  }, mc.cores=ncores)

names(null_regions) <- names(cts_peaks)

#####################################################################################

# Write true and null sequences to fasta files (split by chromosome)
fasta_dir <- paste0(outdir, "/fastas_1000bp_randOnly")
dir.create(fasta_dir, showWarnings = FALSE, recursive = TRUE)

for(ct in names(cts_peaks)){
  # write true fastas
  message(sprintf("Writing true sequences for %s...", ct))
  peaks <- cts_peaks[[ct]]
  peaks <- split(peaks, seqnames(peaks))
  chr_names <- names(peaks)
  for(chr in chr_names){
    true_seqs <- BSgenome::getSeq(genome, names=peaks[[chr]])
    #names(true_seqs) <- peaks[[chr]]$peakName
    names(true_seqs) <- (peaks[[chr]] %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
    writeXStringSet(true_seqs, file=paste0(fasta_dir, sprintf("/%s_%s_true_seqs.fasta", ct, chr)))
  }
  
  # write null fastas
  message(sprintf("Writing null sequences for %s...", ct))
  peaks <- null_regions[[ct]]
  peaks <- split(peaks, seqnames(peaks))
  chr_names <- names(peaks)
  for(chr in chr_names){
    null_seqs <- BSgenome::getSeq(genome, names=peaks[[chr]])
    names(null_seqs) <- (peaks[[chr]] %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
    writeXStringSet(null_seqs, file=paste0(fasta_dir, sprintf("/%s_%s_null_seqs.fasta", ct, chr)))
  }
}

#####################################################################################