#!/usr/bin/env Rscript

#############################
# Processing scalp data
#############################

#Load ArchR (and associated libraries)
library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory
subgroup <- "Keratinocytes"
wd <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/subclustered_%s", subgroup)
full_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/fine_clustered"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 1.0

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_proj <- readRDS(sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/harmonized_subclustering/%s/%s.rds", subgroup, subgroup))

plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

# Colormaps
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
atac_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(atac_proj$Sample2)] %>% unlist()
disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

rna_sub_cmap <- readRDS(paste0(scriptPath, sprintf("/rna_cmap_%s.rds", subgroup)))
atac_sub_cmap <- readRDS(paste0(scriptPath, sprintf("/atac_cmap_%s.rds", subgroup)))

# Load labels from file
source(paste0(scriptPath, "/cluster_labels.R"))

rna_label_cmap <- rna_sub_cmap
names(rna_label_cmap) <- unlist(rna.FineClust)[names(rna_label_cmap)]
atac_label_cmap <- atac_sub_cmap
names(atac_label_cmap) <- unlist(atac.FineClust)[names(atac_label_cmap)]

atac_proj$LFineClust <- unlist(atac.FineClust)[atac_proj$FineClust]
rna_proj$LFineClust <- unlist(rna.FineClust)[rna_proj$FineClust]

###########################################################################################
# Do not proceed prior to calling peaks
###########################################################################################

# Compute group coverages
atac_proj <- addGroupCoverages(
  ArchRProj=atac_proj, 
  groupBy="FineClust", 
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
  )

# Get peaks that were called on this subproject's subclusters from full ArchR project
full_proj <- loadArchRProject(full_dir, force=TRUE)
full_peaks <- getPeakSet(full_proj)
peaks <- getClusterPeaks(full_proj, clusterNames=unique(atac_proj$FineClust), peakGR=full_peaks)
rm(full_proj); gc()

# Now add these peaks to the subproject and generate peak matrix
atac_proj <- addPeakSet(atac_proj, peakSet=peaks, force=TRUE)
atac_proj <- addPeakMatrix(atac_proj, force=TRUE)
atac_proj <- addMotifAnnotations(atac_proj, motifSet="cisbp", name="Motif", force=TRUE)

# Calculate coaccessibility
atac_proj <- addCoAccessibility(
    ArchRProj = atac_proj,
    reducedDims = "Harmony"
)

# Calculate peak-to-gene links
atac_proj <- addPeak2GeneLinks(
    ArchRProj = atac_proj,
    reducedDims = "Harmony"
)

# Add background peaks
atac_proj <- addBgdPeaks(atac_proj, force = TRUE)

# Add deviations matrix
atac_proj <- addDeviationsMatrix(
  ArchRProj = atac_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

###########################################################################################
# Analyses using P2G links and integrated RNA
###########################################################################################

################################
# Plot Peak2Gene heatmap
corcutoff <- 0.45

# Add labels to project
source(paste0(scriptPath, "/cluster_labels.R"))
atac_proj$LabelClust <- unlist(atac.FineClust)[atac_proj$FineClust]

nclust <- 12
p <- plotPeak2GeneHeatmap(
  atac_proj, 
  corCutOff = corcutoff, 
  groupBy="LabelClust", 
  nPlot = 1000000, returnMatrices=FALSE, 
  palGroup=atac_label_cmap, 
  k=nclust, seed=1
  )
pdf(paste0(plotDir, "/peakToGeneHeatmap_LabelClust.pdf"), width=16, height=10)
p
dev.off()
################################

# You need to force it to plot all peaks if you want to match the labeling when you 'returnMatrices'.
p2gMat <- plotPeak2GeneHeatmap(
  atac_proj, 
  corCutOff = corcutoff, 
  groupBy="LabelClust",
  nPlot = 1000000, returnMatrices=TRUE, 
  k=nclust, seed=1) # Need to set the nPlot number correctly

# Get association of peaks to clusters
kclust_df <- data.frame(
  kclust=p2gMat$ATAC$kmeansId,
  peakName=p2gMat$Peak2GeneLinks$peak,
  gene=p2gMat$Peak2GeneLinks$gene
  )

# Fix peakname
kclust_df$peakName <- sapply(kclust_df$peakName, function(x) strsplit(x, ":|-")[[1]] %>% paste(.,collapse="_"))

# Get full peakset
peaksGR <- getPeakSet(atac_proj)
names(peaksGR) <- (peaksGR %>% {paste(seqnames(.), start(.), end(.), sep="_")})

# Get motif matches
matches <- getMatches(atac_proj, "Motif")
r1 <- SummarizedExperiment::rowRanges(matches)
rownames(matches) <- paste(seqnames(r1),start(r1),end(r1),sep="_")
matches <- matches[names(peaksGR)]

clusters <- unique(kclust_df$kclust) %>% sort()

enrichList <- lapply(clusters, function(x){
  cPeaks <- kclust_df[kclust_df$kclust == x,]$peakName %>% unique()
  ArchR:::.computeEnrichment(matches, which(names(peaksGR) %in% cPeaks), seq_len(nrow(matches)))
  }) %>% SimpleList
names(enrichList) <- clusters

# Format output to match ArchR's enrichment output
assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
    d <- lapply(seq_along(enrichList), function(y){
        enrichList[[y]][colnames(matches),x,drop=FALSE]
      }) %>% Reduce("cbind",.)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
names(assays) <- colnames(enrichList[[1]])
assays <- rev(assays)
res <- SummarizedExperiment::SummarizedExperiment(assays=assays)

formatEnrichMat <- function(mat, topN, minSig, clustCols=TRUE){
  plotFactors <- lapply(colnames(mat), function(x){
    ord <- mat[order(mat[,x], decreasing=TRUE),]
    ord <- ord[ord[,x]>minSig,]
    rownames(head(ord, n=topN))
  }) %>% do.call(c,.) %>% unique()
  pMat <- mat[plotFactors,]
  prettyOrderMat(pMat, clusterCols=clustCols)$mat
}

pMat <- formatEnrichMat(assays(res)$mlog10Padj, 5, 10, clustCols=TRUE)

pdf(paste0(plotDir, "/enrichedMotifs_kclust_p2gHM.pdf"), width=6, height=10)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  pMat, 
  limits=c(0,100), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  row_names_side = "left",
  width = ncol(pMat)*unit(0.5, "cm"),
  height = nrow(pMat)*unit(0.33, "cm"),
  border_gp=gpar(col="black"), # Add a black border to entire heatmap
  legendTitle="-log10 FDR"
  )
draw(hm)
dev.off()

# And finally, GO enrichments of genes per cluster
source(paste0(scriptPath, "/GO_wrappers.R"))

kclust <- unique(kclust_df$kclust) %>% sort()
all_genes <- kclust_df$gene %>% unique() %>% sort()

GOresults <- lapply(kclust, function(k){
  message(sprintf("Running GO enrichments on k cluster %s...", k))
  clust_genes <- kclust_df[kclust_df$kclust == k,]$gene %>% unique()
  upGO <- rbind(
    calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="BP"), 
    calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="MF")
    #calcTopGo(all_genes, interestingGenes=upGenes, nodeSize=5, ontology="CC")
    )
  upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
  })

names(GOresults) <- paste0("cluster_", kclust)

# Plots of GO term enrichments:
pdf(paste0(plotDir, "/kclust_GO_terms.pdf"), width=8, height=5)
for(name in names(GOresults)){
    goRes <- GOresults[[name]]
    if(nrow(goRes)>1){
      print(topGObarPlot(goRes, cmap = cmaps_BOR$comet, 
        nterms=10, border_color="black", 
        barwidth=0.85, title=name))
    }
}
dev.off()

##########################################################################################

