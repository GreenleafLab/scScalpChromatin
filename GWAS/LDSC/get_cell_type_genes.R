#!/usr/bin/env Rscript

##########################################################
# Generate an region bed files for FineMapped SNP analyses
##########################################################

# This script will just identify 'cell type specific genes' to be used for FineMapped SNP analyses
# cell type specific genes genes will be compared to cell type specific peaks to contrast
# genes vs chromatin information for SNP enrichment/assignment analyses


#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(Seurat)
  library(ArchR)
  library(future)
  library(dplyr)
  library(tidyr)
  library(rtracklayer)
})

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/plotting_config.R"))

# change the current plan to access parallelization (for Seurat)
nThreads <- 8
plan("multicore", workers = nThreads)

message("Generating cell type-specific gene sets...")

ident <- "FineClust"
atac_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/fine_clustered"
rna_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/preprocessing_output"
outdir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/GWAS/rna_ldsc/"
minCells <- 40
minGenes <- 200
windowSize <- 100000  # The size of the window to expand around cell type specific genes
# (In the 2018 NG Finacune paper, they added 100kb around each gene body prior to performing LDSC.
# This makes very little difference in the actual results compared to not expanding the window.)

# set output directory
outdir <- paste0(outdir, sprintf("/%s_specific_genes", ident))

#Set/Create Working Directory to Folder
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# Load Seurat project
message("Reading in data...")
obj <- readRDS(paste0(rna_dir, '/scalp.rds'))

# Load ArchR project to get genes GR (to ensure compatibility)
atac_proj <- loadArchRProject(atac_dir, force=TRUE)
allGenesGR <- getGenes(atac_proj)

# Drop genes with 'NA' symbol
allGenesGR <- allGenesGR[!is.na(allGenesGR$symbol)]
names(allGenesGR) <- allGenesGR$symbol

#########################################################
# Gene selection using all genes filtered for overlap
#########################################################

allGenes <- allGenesGR$symbol
allGenes <- unname(allGenes)

# We will first identify marker genes using Seurat's marker gene identification,
# then we will filter for cts-specific genes in a manner analagous to how we filtered cts-peaks

ccd <- obj@meta.data
clust_freqs <- getFreqs(ccd[,ident])
use_clust <- names(clust_freqs)[clust_freqs > minCells]

# Set the max number of overlapping clusters per gene to be 1/4 the total number of clusters
cutoff <- floor(length(use_clust)/4)

# Use Seurat's FindMarkers to get marker genes for each cluster
message("Finding marker genes using Seurat...")
Idents(obj) <- ident
# (Default is min.pct=0.1 and logfc.threshold=0.25)
obj.markers <- FindAllMarkers(obj, only.pos=TRUE, 
  min.pct=0.1, logfc.threshold=0.25, genes.use=allGenes)

# Remove clusters that do not have any marker genes
use_clust <- use_clust[use_clust %in% unique(obj.markers$cluster)]

# Remove marker genes not in all genes
obj.markers <- obj.markers[obj.markers$gene %in% allGenes,]

# Construct matrix of gene usage by cluster
gene_mat <- matrix(0, nrow=length(allGenes), ncol=length(use_clust))
rownames(gene_mat) <- allGenes
colnames(gene_mat) <- use_clust

# Fill in matrix
for(clust in use_clust){
  gene_mat[obj.markers[obj.markers$cluster == clust,7], clust] <- 1
}

# Determine overlap for each gene
gene_usage <- Matrix::rowSums(gene_mat)
names(gene_usage) <- rownames(gene_mat)

# Plot peak-sharing information
nclust_per_gene <- getFreqs(gene_usage)
df <- data.frame(nClust=as.integer(names(nclust_per_gene)), nGenes=nclust_per_gene)
df <- df[order(df$nClust),]

pdf(paste0(outdir, "/nClusters_per_gene.pdf"), width=8, height=6)
qcBarPlot(df, cmap="royalblue1", barwidth=0.9, border_color=NA) + geom_vline(xintercept=cutoff + 0.5, linetype="dashed")
dev.off()

# Percent of total peaks
df$cumulative <- cumsum(df$nGenes)
df$cumulative <- df$cumulative / max(df$cumulative)

# Need to subtract the genes associated with 0 clusters
df$cumulative <- df$cumulative - df[1,3]

pct_genes_passing <- df[df$nClust == cutoff,]$cumulative %>% round(.,2)

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
  + scale_y_continuous(expand=c(0, 0), limits=c(0,1)) # Make bars start at the axis
  )
p <- p + geom_vline(xintercept=cutoff, linetype="dashed") + geom_hline(yintercept=pct_genes_passing, linetype="dashed")
p <- p + annotate("text",  x=-Inf, y = Inf, label = paste0(pct_genes_passing*100, "%"), vjust=1, hjust=-1)

pdf(paste0(outdir, "/cumulative_pct_genes.pdf"), width=8, height=6)
p
dev.off()

# Identify which genes to use for each cluster
filt_gene_mat <- gene_mat[gene_usage <= cutoff,]

markerGRList <- lapply(use_clust, function(clust){
  clust_vec <- filt_gene_mat[,clust]
  allGenesGR[names(clust_vec[clust_vec>0])]
  })
names(markerGRList) <- use_clust

# Filter clusters that had too few genes:
markerGRList <- markerGRList[lapply(markerGRList, length) > minGenes]

# Add a GR of all marker genes to use as background
markerGRList[["allGenes"]] <- allGenesGR[unique(obj.markers$gene)] %>% sort(., ignore.strand=TRUE)

# Save coordinate bed files to use for generating annot files
for(ctype in names(markerGRList)){

  # First, save un-expanded GR
  gr <- markerGRList[[ctype]]
  saveRDS(gr, file=paste0(outdir, sprintf("/%s_specific_genes.rds", ctype)))

  # Now, save expanded range bed file
  outbed <- paste0(outdir, sprintf("/%s_specific_genes.bed", ctype))
  gr <- sort(gr, ignore.strand=TRUE) %>% 
    suppressWarnings(resize(., width=width(.)+windowSize, fix="center", ignore.strand=TRUE)) %>% 
    reduce(., ignore.strand=TRUE) %>% trim()
  df <- DataFrame(
    chrom=seqnames(gr),
    start=start(gr)-1,
    end=end(gr)
  )
  write.table(df, outbed, quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)
}
message("Done")


##########################################################
