#!/usr/bin/env Rscript

##########################################################################
# Generate an region bed files to make annot files for LD score regression
##########################################################################

# See here: https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
# Annot file will be made from peaks identified in scATAC data
# Annot file has the following columns:
# CHR, BP, SNP, and CM (centimorgan), followed by one column per annotation
# The value of the annotation for each SNP can be 0/1 for binary categories, 
# or arbitrary numbers for continuous annotations.
# The file can have many categories, or just a single category

# It MUST have the same SNPs in the same order as the .bim file used for the 
# computation of LD scores

# This script will just prepare the 'gene_coord_file' to be used in 'make_annot.py' from ldsc
# Usage: Rscript get_cell_type_peaks.R <ident> <ArchR_dir> <out_dir>
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 3){
  stop("Must supply ident for getting peaks from ArchR project! \n \
    Usage: Rscript get_cell_type_peaks.R <ident> <ArchR_dir> <out_dir>")
}

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(rtracklayer)
})

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/plotting_config.R"))

message("Generating cell type-specific peak sets...")

# Set Threads to be used
addArchRThreads(threads = 8)

ident <- args[1]
archrDir <- args[2]
outdir <- args[3]
minPeaks <- 5000
minCells <- 40
resizePeakWidth <- 500
liftToHg19 <- FALSE # Should peak coordinates be lifted over to hg19?

# set output directory
outdir <- paste0(outdir, sprintf("/%s_specific_peaks", ident))

#Set/Create Working Directory to Folder
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load ArchR project
proj <- loadArchRProject(archrDir, force=TRUE)

# Load full peakset
peaksGR <- getPeakSet(proj)
names(peaksGR) <- (peaksGR %>% {paste(seqnames(.), start(.), end(.), sep="_")})

#########################################################
# Peak selection using all peaks filtered for overlap
#########################################################

ccd <- proj@cellColData
clust_freqs <- getFreqs(ccd[,ident])
use_clust <- names(clust_freqs)[clust_freqs > minCells]

# Set the max number of overlapping clusters per peak to be 1/4 the total number of clusters
cutoff <- floor(length(use_clust)/4)

cts_peaks <- lapply(use_clust, function(x) getClusterPeaks(proj, clusterNames=x, peakGR=peaksGR))
names(cts_peaks) <- use_clust

# Construct matrix of peak usage by cluster
peak_mat <- matrix(0, nrow=length(peaksGR), ncol=length(use_clust))
rownames(peak_mat) <- names(peaksGR)
colnames(peak_mat) <- use_clust

# Fill in matrix
for(clust in use_clust){
  peak_mat[names(cts_peaks[[clust]]), clust] <- 1
}

# Determine overlap for each peak
peak_usage <- Matrix::rowSums(peak_mat)
names(peak_usage) <- rownames(peak_mat)

# Plot peak-sharing information
nclust_per_peak <- getFreqs(peak_usage)
df <- data.frame(nClust=as.integer(names(nclust_per_peak)), nPeaks=nclust_per_peak)
df <- df[order(df$nClust),]

pdf(paste0(outdir, "/nClusters_per_peak.pdf"), width=8, height=6)
qcBarPlot(df, cmap="royalblue1", barwidth=0.9, border_color=NA) + geom_vline(xintercept=cutoff + 0.5, linetype="dashed")
dev.off()

# Percent of total peaks
df <- rbind(c(0,0), df)
df$cumulative <- cumsum(df$nPeaks)
df$cumulative <- df$cumulative / max(df$cumulative)
pct_peaks_passing <- df[df$nClust == cutoff,]$cumulative %>% round(.,2)

pcol <- "royalblue1"
p <- (
  ggplot(df, aes(x=df[,1], y=df[,3]))
  + geom_line(size=2, color=pcol)
  + scale_fill_manual(values = pcol)
  + xlab(colnames(df)[1])
  + ylab(colnames(df)[3])
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          axis.text.x = element_text(angle = 90, hjust = 1))
  + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
  )
p <- p + geom_vline(xintercept=cutoff, linetype="dashed") + geom_hline(yintercept=pct_peaks_passing, linetype="dashed")
p <- p + annotate("text",  x=-Inf, y = Inf, label = paste0(pct_peaks_passing*100, "%"), vjust=1, hjust=-1)

pdf(paste0(outdir, "/cumulative_pct_peaks.pdf"), width=8, height=6)
p
dev.off()

# Identify which peaks to use for each cluster
filt_peak_mat <- peak_mat[peak_usage <= cutoff,]

markerGRList <- lapply(use_clust, function(clust){
  clust_vec <- filt_peak_mat[,clust]
  peaksGR[names(clust_vec[clust_vec>0])]
  })
names(markerGRList) <- use_clust

# Filter clusters that had too few peaks:
markerGRList <- markerGRList[lapply(markerGRList, length) > minPeaks]

# Add a GR of all peaks to use as background
markerGRList[["allPeaks"]] <- getPeakSet(proj)


if(liftToHg19){
  message("Lifting over peak coordinates from hg38 to hg19...")
  chain <- import.chain("/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/liftover/hg38ToHg19.over.chain")
  nms <- names(markerGRList)
  markerGRList <- lapply(names(markerGRList), function(x){
    liftOver(markerGRList[[x]], chain) %>% unlist()
    })
  names(markerGRList) <- nms
}

# Save coordinate bed files to use for generating annot files
for(ctype in names(markerGRList)){
  outfile <- paste0(outdir, sprintf("/%s_specific_peaks.bed", ctype))
  gr <- markerGRList[[ctype]] %>% sort() %>% resize(., width=resizePeakWidth, fix="center") %>% reduce(.)
  df <- DataFrame(
    chrom=seqnames(gr),
    start=start(gr)-1,
    end=end(gr)
  )
  write.table(df, outfile, quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)
}
message("Done")

