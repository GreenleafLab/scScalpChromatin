#!/usr/bin/env Rscript

#######################################################################################
# Analyses on keratinocytes
# (This assumes you have already performed peak assignments, RNA integration, etc.)
#######################################################################################

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

###########################################################################################
# Do not proceed prior to calling peaks on whole 're-merged' project
# Also do not proceed prior to assigning peaks to keratinocyte sub-project
###########################################################################################

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

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

####################################################################################
# Plot UMAPs of each disease and each sample independently
####################################################################################

## ATAC:

na.color <- "grey80"

ccd <- atac_proj@cellColData

pList <- list()
for(dis in c("AA", "C_SD", "C_PB")){

    # Disease alone (single color)
    umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="diseaseStatus")
    colnames(umapDF) <- c("UMAP1", "UMAP2", "value")
    use_vals <- ccd[rownames(umapDF), "diseaseStatus"]
    umapDF$value <- ifelse(umapDF$value == dis, use_vals, NA)

    # Randomize cells before plotting UMAP
    set.seed(1)
    umapDF <- umapDF[sample(nrow(umapDF)),] %>% arrange(desc(is.na(value)))
    pList[[paste0(dis, "_only")]] <- plotUMAP(umapDF, dataType="qualitative", cmap=disease_cmap, 
      namedColors=TRUE, point_size=pointSize, covarLabel="diseaseStatus", na.value=na.color)

    # Disease by sample (multi-color)
    umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="diseaseStatus")
    colnames(umapDF) <- c("UMAP1", "UMAP2", "value")
    use_vals <- ccd[rownames(umapDF), "Sample2"]
    umapDF$value <- ifelse(umapDF$value == dis, use_vals, NA)
    # Randomize cells before plotting UMAP
    set.seed(1)
    umapDF <- umapDF[sample(nrow(umapDF)),] %>% arrange(desc(is.na(value)))
    pList[[paste0(dis, "_bySamp")]] <- plotUMAP(umapDF, dataType="qualitative", cmap=atac_sample_cmap, 
      namedColors=TRUE, point_size=pointSize, covarLabel="Sample", na.value=na.color)

    # Disease by Cluster (multi-color)
    umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="diseaseStatus")
    colnames(umapDF) <- c("UMAP1", "UMAP2", "value")
    use_vals <- ccd[rownames(umapDF), "LFineClust"]
    umapDF$value <- ifelse(umapDF$value == dis, use_vals, NA)
    # Randomize cells before plotting UMAP
    set.seed(1)
    umapDF <- umapDF[sample(nrow(umapDF)),] %>% arrange(desc(is.na(value)))
    pList[[paste0(dis, "_byClust")]] <- plotUMAP(umapDF, dataType="qualitative", cmap=atac_label_cmap, 
      namedColors=TRUE, point_size=pointSize, covarLabel="Cluster", na.value=na.color)
}

pdf(paste0(plotDir, "/disease_UMAPs_atac.pdf"), width=7, height=7)
pList
dev.off()

## RNA:

pList <- list()
for(dis in c("AA", "C_SD", "C_PB")){

    # Disease alone (single color)
    umapDF <- data.frame(Embeddings(object=rna_proj, reduction="umap"), 
        ifelse(rna_proj$diseaseStatus == dis, rna_proj$diseaseStatus, NA))
    colnames(umapDF) <- c("UMAP1", "UMAP2", "value")
    # Randomize cells before plotting UMAP
    set.seed(1)
    umapDF <- umapDF[sample(nrow(umapDF)),] %>% arrange(desc(is.na(value)))
    pList[[paste0(dis, "_only")]] <- plotUMAP(umapDF, dataType="qualitative", cmap=disease_cmap, 
        point_size=pointSize, namedColors=TRUE, na.value=na.color)

    # Disease by sample (multi-color)
    umapDF <- data.frame(Embeddings(object=rna_proj, reduction="umap"), 
        ifelse(rna_proj$diseaseStatus == dis, rna_proj$Sample, NA))
    colnames(umapDF) <- c("UMAP1", "UMAP2", "value")
    # Randomize cells before plotting UMAP
    set.seed(1)
    umapDF <- umapDF[sample(nrow(umapDF)),] %>% arrange(desc(is.na(value)))
    pList[[paste0(dis, "_bySamp")]] <- plotUMAP(umapDF, dataType="qualitative", cmap=sample_cmap, 
        point_size=pointSize, namedColors=TRUE, na.value=na.color)

    # Disease by Cluster (multi-color)
    umapDF <- data.frame(Embeddings(object=rna_proj, reduction="umap"), 
        ifelse(rna_proj$diseaseStatus == dis, rna_proj$LFineClust, NA))
    colnames(umapDF) <- c("UMAP1", "UMAP2", "value")
    # Randomize cells before plotting UMAP
    set.seed(1)
    umapDF <- umapDF[sample(nrow(umapDF)),] %>% arrange(desc(is.na(value)))
    pList[[paste0(dis, "_byClust")]] <- plotUMAP(umapDF, dataType="qualitative", cmap=rna_label_cmap, 
        point_size=pointSize, namedColors=TRUE, na.value=na.color)
}

pdf(paste0(plotDir, "/disease_UMAPs_rna.pdf"), width=7, height=7)
pList
dev.off()


####################################################################################
# Milo on scATAC
####################################################################################

library(miloR)
library(SingleCellExperiment)

# Convert ArchR project into a singleCellExperiment

# Subset Milo object to remove samples that are too lowly represented:
ccd <- atac_proj@cellColData
samp_freqs <- getFreqs(ccd$Sample2)
valid_samps <- samp_freqs[samp_freqs > 50] %>% names()
valid_cells <- rownames(ccd[ccd$Sample2 %in% valid_samps,])

# First get just the peak matrix (what we will use as our 'counts' matrix)
counts_mat <- getMatrixFromProject(atac_proj, useMatrix="PeakMatrix")
counts_mat <- counts_mat[,valid_cells]

sce_atac <- SingleCellExperiment(counts_mat)
names(assays(sce_atac)) <- names(assays(counts_mat))
dim_reduc <- atac_proj@reducedDims$Harmony$matDR
umap <- atac_proj@embeddings$UMAP$df
reducedDims(sce_atac)[["LSI"]] <- dim_reduc[valid_cells,]
reducedDims(sce_atac)[["UMAP"]] <- umap[valid_cells,]
colData(sce_atac) <- colData(counts_mat)

# Compare AA to C (all controls)
sce_atac$diseaseStatus <- factor(ifelse(sce_atac$diseaseStatus == "AA", "AA", "C"), levels=c("C", "AA"), ordered=TRUE)

# Create Milo object
set.seed(123)
milo_proj <- Milo(sce_atac)

# build KNN graph using LSI reduced dimensions
d <- ncol(reducedDims(sce_atac)$LSI)
k <- 30 # 30 is default
reduced.dims <- "LSI" # Milo capitalizes the dim reduc names
milo_proj <- buildGraph(milo_proj, k=k, d=d, reduced.dim=reduced.dims)

# Defining representative neighborhoods on the KNN graph
prop <- 0.1 # "We suggest using prop=0.1 for datasets of less than 30k cells. For bigger datasets using prop=0.05 should be sufficient (and makes computation faster)"
milo_proj <- makeNhoods(milo_proj, prop=prop, k=k, d=d, refined=TRUE, reduced_dims=reduced.dims)

pdf(paste0(plotDir, "/miloR_nhood_size_hist.pdf"), width=8, height=6)
plotNhoodSizeHist(milo_proj)
dev.off()

# Counting cells in neighborhoods
milo_proj <- countCells(milo_proj, meta.data=as.data.frame(colData(milo_proj)), sample="Sample2")

# Defining experimental design
milo_design <- data.frame(colData(milo_proj))[,c("Sample2", "diseaseStatus")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$Sample

# Computing neighborhood connectivity
milo_proj <- calcNhoodDistance(milo_proj, d=d, reduced.dim=reduced.dims)

# Testing differential abundance
da_results <- testNhoods(milo_proj, design = ~ diseaseStatus, design.df=milo_design, reduced.dim=reduced.dims)

# Visualize results:
milo_proj <- buildNhoodGraph(milo_proj)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo_proj, da_results, layout="UMAP", alpha=0.1) 

pdf(paste0(plotDir, sprintf("/miloR_atac_DA_UMAP_%s.pdf", subgroup)), width=7, height=6)
nh_graph_pl
dev.off()

# Save differential abundance results

# First need to map Milo neighborhoods to majority cluster IDs
nhoodmat <- milo_proj@nhoods
cell_to_clust <- atac_proj$FineClust
names(cell_to_clust) <- atac_proj$cellNames

# For each nhood, identify the most common cluster
nhood_to_clust <- apply(nhoodmat, 2, function(x){
  cellnames <- names(x[x>0]) # Which cells are part of this nhood
  names(sort(table(cell_to_clust[cellnames]), decreasing=TRUE))[1] # Most frequent cluster
  })

# Assign majority cluster to nhoods
da_results$majority_cluster <- nhood_to_clust[da_results$Nhood]
da_results$majority_lclust <- unlist(atac.FineClust)[da_results$majority_cluster]
da_results <- da_results[order(da_results$PValue, decreasing=FALSE),]

# Save table
table_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/supplemental_tables"
write.table(da_results, file=paste0(table_dir, "/keratinocyte_subclustered_Milo_results.tsv"), quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

####################################################################################
# Marker Genes
####################################################################################

markerGenes  <- c(
  # https://www.proteinatlas.org/humanproteome/tissue/skin
  "KRT14", "KRT5", "KRT15", "COL17A1", # Basal epithelia
  "KRT1", "KRT10", # Spinous epithelia
  "DSC1", "KRT2", "IVL", # Granular epithelia 
  "KRT7", # Glandular
  "ITGA6", "ITGB1", "LGR5","LHX2", "FRZB", "IL31RA", "OSMR", # HFSCs?
  "CD34", "CDH3", "LGR5", "CDKN2A", # Hair germ markers (some say HG is CD34 neg?)
  "KRT81", "KRT83", "HOXC13", "LEF1", # Matrix hair keratins/genes
  "KRT71", # IRS keratins / genes
  "KRT75", # ORS keratins / genes
  "AR", "CD200", "LGR6", "LRIG1", # Isthmus markers
  "ELOVL3", "PRDM1", "IHH", # Sebaceous (PRDM1 = Blimp1) (IHH = Indian hedgehog)
  "KRT23", "KRT18", "KRT19", "EDAR", # Mystery cluster markers
  "SOX9", "LHX2", "NFATC1", "TCF3", # Key HFSC TFs
  "HLA-A", "HLA-B", "HLA-C", # MHC1
  "WNT3", "WNT5A", "WNT10A", "WNT10B", # WNTs
  "FOXC1", "POU2F3", "TP63", "KLF4", "KLF5",  # Other interesting skin TFs
  "RUNX1", "RUNX2", "RUNX3",
  "NFIC", "FOSL1", "FOSL2", "FOXN3",
  "EGR1", "EGR2", "EGR3", "EGR4"
) %>% unique()

# GeneScores first:

# Identify Marker Gene through Pairwise Test vs Bias-Matched Background:
markersGS <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "LFineClust",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersGS)$name

# Lists of 'marker genes'
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.50")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.0", 
  labelMarkers = markerGenes,
  transpose = FALSE # Can't transpose if we want pretty ordered heatmaps...
)

draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 6, height = 8, ArchRProj = atac_proj, addDOC = FALSE)


p <- plotEmbedding(
    ArchRProj = atac_proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(atac_proj), 
    size = pointSize, plotAs = "points"
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)

# Now Imputed RNA:

markersPRNA <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "GeneIntegrationMatrix", 
    groupBy = "LFineClust",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersPRNA)$name %>% unique()

# Lists of pseudo-RNA 'marker genes'
markerList <- getMarkers(markersPRNA, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")


##############################################
# ("Top" genes are defined as having the most peak to gene links)
source(paste0(scriptPath, "/GO_wrappers.R"))

markers <- markerList$HF.Kc_3$name %>% head(100)

upGO <- rbind(
    calcTopGo(geneNames, interestingGenes=markers, nodeSize=5, ontology="BP"),
    calcTopGo(geneNames, interestingGenes=markers, nodeSize=5, ontology="MF")
    #calcTopGo(all_genes, interestingGenes=upGenes, nodeSize=5, ontology="CC")
)
upGO <- upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]

##############################################

heatmapPRNA <- plotMarkerHeatmap(
  seMarker = markersPRNA, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.5", 
  labelMarkers = markerGenes,
  transpose = FALSE # Can't transpose if we want pretty ordered heatmaps...
)

draw(heatmapPRNA, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPRNA, name = "pseudoRNA-Marker-Heatmap", width = 6, height = 8, ArchRProj = atac_proj, addDOC = FALSE)


#Plot the UMAP Embedding with pseudoRNA Marker Genes Overlayed
p <- plotEmbedding(
    ArchRProj = atac_proj, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes,
    pal = cmaps_BOR$sunrise,
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL, 
    plotAs="points", size = pointSize
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-pseudoRNA_Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)


# Tracks of genes:

atacOrder <- c(
    "aKc1", # "Basal.Kc_1",
    "aKc2", # "Spinous.Kc_2",
    "aKc3", # "Spinous.Kc_1",
    "aKc4", # "Infundibulum", # SOX9, DKK3
    "aKc5", # "Inf.Segment_1", # Lhx2, LGR5 high
    "aKc7", # "Inf.Segment_2", # Lhx2, LGR5 high
    "aKc6", # "Sebaceous", 
    "aKc8", # "Isthmus", # CD200 high
    "aKc9", # "Matrix", 
    "aKc10" # "Eccrine"
    #"aKc11", # "Unknown",
)

p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "LFineClust", 
    useGroups = unlist(atac.FineClust)[atacOrder],
    pal = atac_label_cmap,
    plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
    sizes = c(6, 0.2, 1.25, 0.5),
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(atac_proj),
    tileSize=250,
    minCells=100
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 6, height = 7)

##########################################################################################
# Identifying Marker Peaks
##########################################################################################

# Exclude certain groups
all_LFineClust <- unique(atac_proj$LFineClust)
exclude <- c("Unknown_1")
use_groups <- all_LFineClust[all_LFineClust %ni% exclude]

# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "LFineClust",
    useGroups = use_groups,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

#Visualize Markers as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  nLabel = 1, # It still seems like there's not actually a way to NOT plot any labels
  binaryClusterRows = TRUE,
  transpose = FALSE
)
draw(heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapPeaks, name="Peak-Marker-Heatmap", width=10, height=15, ArchRProj=atac_proj, addDOC=FALSE)

##########################################################################################
# Motif Enrichments
##########################################################################################

#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = atac_proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Rename motifs for more aesthetic plotting:
rownames(enrichMotifs) <- lapply(rownames(enrichMotifs), function(x) strsplit(x, "_")[[1]][1]) %>% unlist()

# Subset to clusters that have at least some enrichment
log10pCut <- 5
keepClust <- colnames(enrichMotifs)[apply(assays(enrichMotifs)$mlog10Padj, 2, max) > log10pCut]

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(enrichMotifs[,keepClust], n=5, transpose=TRUE, cutOff=log10pCut)
#heatmapEM <- plotEnrichHeatmap(enrichMotifs, n=5, transpose=TRUE, cutOff=log10pCut)
draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapEM, name="Motifs-Enriched-Marker-Heatmap", width=12, height=8, ArchRProj=atac_proj, addDOC=FALSE)

#############################
# Motif Footprinting
#############################

# Just a few that might be particularly interesting 
motifPositions <- getPositions(atac_proj)

motifs <- c(
  "AR_", "TP63", "SOX4", "SOX9", "JUND", "FOS", "BATF", "FOXC1", "LHX2",
  "PRRX1", "PRRX2", "POU2F3", "RUNX3", "CEBPB", "HOXA1", "HOXC13", # Highly Regulated TFs
  "NFIC", "KLF4", "POU2F3", "EGR1", "EGR2", "EGR3", "EGR4"
  )
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = atac_proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "FineClust"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = atac_proj, 
  pal = atac_sub_cmap,
  normMethod = "Divide",
  plotName = "FineClust-Footprints",
  addDOC = FALSE,
  smoothWindow = 5
)
###########################

plotVarDev <- getVarDeviations(atac_proj, plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = atac_proj, addDOC = FALSE)

markerMotifs <- getFeatures(atac_proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

# Plot over UMAP:

p <- plotEmbedding(
    ArchRProj = atac_proj, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(atac_proj),
    plotAs="points", size = pointSize
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-ChromVar.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, width = 5, height = 5)

####################################################################################
# Now further subcluster HF cells
####################################################################################

subClusterArchR <- function(proj, subCells, outdir){
  # Subset an ArchR project for focused analysis

  message(sprintf("Subgroup has %s cells.", length(subCells)))
  sub_proj <- subsetArchRProject(
      ArchRProj = proj,
      cells = subCells,
      outputDirectory = outdir,
      dropCells = TRUE,
      force = TRUE
  )
  saveArchRProject(sub_proj)

  # Delete stuff we don't want to copy...
  unlink(paste0(outdir, "/Plots/*"), recursive = TRUE)
  unlink(paste0(outdir, "/IterativeLSI/*"), recursive = TRUE)
  unlink(paste0(outdir, "/Embeddings/*"), recursive = TRUE)
  unlink(paste0(outdir, "/ImputeWeights/*"), recursive = TRUE)
}

# Hair follicle keratinocytes
sg <- "HF_Kc"
sg_clust <- c("HF.Kc_2", "HF.Kc_4", "HF.Kc_6")
sg_cells <- getCellNames(atac_proj)[as.character(atac_proj@cellColData[["LFineClust"]]) %in% sg_clust]

outdir <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/subclustered_%s", sg)
subClusterArchR(atac_proj, subCells=sg_cells, outdir=outdir)

####################################################################################
# IFE keratinocyte trajectory
####################################################################################

atac_proj <- addSlingShotTrajectories(
  ArchRProj = atac_proj,
  name = "IFE",
  useGroups = c("Basal.Kc_1", "Spinous.Kc_1", "Spinous.Kc_2"),
  principalGroup = "Basal.Kc_1",
  groupBy = "LFineClust",
  embedding = "UMAP",
  force = TRUE,
  seed = 1
)


# Plot trajectory
p1 <- plotTrajectory(
  atac_proj, 
  trajectory = "IFE.Curve1", 
  embedding = "UMAP",
  colorBy = "cellColData", 
  name = "IFE.Curve1", 
  continuousSet = "fireworks2",
  plotAs = "points",
  smoothWindow = 5,
  size = pointSize
  )

plotPDF(p1, name="IFE_slingshot_trajectory.pdf", ArchRProj=atac_proj, addDOC=FALSE, width=5, height=5)

####################################################################################
# More trajectory plots
####################################################################################

# Get list of genes we want to highlight (e.g. genes involved in HF development)
library(org.Hs.eg.db)
library(GO.db)
go_id = GOID(GOTERM[Term(GOTERM) == "cornification"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
cornification_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

go_id = GOID(GOTERM[Term(GOTERM) == "hemidesmosome assembly"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
hemidesmosome_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

go_id = GOID(GOTERM[Term(GOTERM) == "hair follicle development"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
hfdev_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

label_genes <- c(cornification_genes, hemidesmosome_genes, markerGenes) %>% unique() %>% sort()

# ChromVAR motifs:
IFE.trajMM  <- getTrajectory(ArchRProj=atac_proj, name="IFE.Curve1", 
    useMatrix="MotifMatrix", log2Norm=FALSE, groupEvery=1.5)

# Dump the deviations rows
IFE.trajMM <- IFE.trajMM[grepl("^z:", rownames(IFE.trajMM))]

# For some reason, it adds the z to the TF
rownames(IFE.trajMM) <- rownames(IFE.trajMM) %>% strsplit(":") %>% sapply('[',2) %>% strsplit("_") %>% sapply('[',1)

topVarTFs <- assays(IFE.trajMM)$smoothMat
topVarTFs <- rownames(topVarTFs[order(rowVars(topVarTFs), decreasing=TRUE),]) %>% head(10)
labelTFs <- c("TP63", "KLF3", "KLF4", "KLF5", "RUNX1", "LHX2", "POU2F3", "NFIC", "RORA", "REL", topVarTFs) %>% unique()

var.cutoff <- 0.9 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)

p1 <- plotTrajectoryHeatmap(IFE.trajMM, pal=paletteContinuous(set="solarExtra"), 
    varCutOff=var.cutoff, limits=lims, labelMarkers=labelTFs, labelTop=0)

plotPDF(p1, name="IFE_MM_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=7, height=10)

# GeneScoreMatrix :
IFE.trajGSM  <- getTrajectory(ArchRProj=atac_proj, name="IFE.Curve1", 
    useMatrix="GeneScoreMatrix", log2Norm=TRUE, groupEvery=1.5) 

# For some reason, it adds the chromosome to the gene name
rownames(IFE.trajGSM) <- rownames(IFE.trajGSM) %>% strsplit(":") %>% sapply('[',2)

var.cutoff <- 0.9 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)
labeltop <- 10 # Default is 50

p1 <- plotTrajectoryHeatmap(IFE.trajGSM, pal=paletteContinuous(set="horizonExtra"), 
    varCutOff=var.cutoff, limits=lims, labelMarkers=label_genes, labelTop=labeltop)

plotPDF(p1, name="IFE_GSM_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=7, height=10)

# GeneIntegrationMatrix :
IFE.trajGIM  <- getTrajectory(ArchRProj=atac_proj, name="IFE.Curve1", 
    useMatrix="GeneIntegrationMatrix", log2Norm=TRUE, groupEvery=1.5)

rownames(IFE.trajGIM) <- rownames(IFE.trajGIM) %>% strsplit(":") %>% sapply('[',2)

var.cutoff <- 0.9 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)

p1 <- plotTrajectoryHeatmap(IFE.trajGIM, pal=cmaps_BOR$sunrise, 
    varCutOff=var.cutoff, limits=lims, labelMarkers=label_genes, 
    labelTop=labeltop, returnMatrix=FALSE) # set returnMatrix to false for plotting

plotPDF(p1, name="IFE_GIM_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=7, height=10)

# Peaks :
IFE.trajPeaks  <- getTrajectory(ArchRProj=atac_proj, name="IFE.Curve1", 
    useMatrix="PeakMatrix", log2Norm=TRUE, groupEvery=1.5) 


var.cutoff <- 0.9 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)

p1 <- plotTrajectoryHeatmap(IFE.trajPeaks, pal=paletteContinuous(set="solarExtra"), 
    varCutOff=var.cutoff, limits=lims, labelTop=5, 
    maxFeatures = 100000, returnMatrix=FALSE) # set returnMatrix to false for plotting

plotPDF(p1, name="IFE_PeakMatrix_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=7, height=10)

####################################################################################
# Integrative pseudo-time analysis
####################################################################################

# GeneScoreMatrix:
corGSM_MM <- correlateTrajectories(IFE.trajGSM, IFE.trajMM)

trajGSM2 <- IFE.trajGSM[corGSM_MM[[1]]$matchname1, ]
trajMM2 <- IFE.trajMM[corGSM_MM[[1]]$name2, ]

trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat=TRUE, varCutOff=0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal=paletteContinuous(set="horizonExtra"),  varCutOff=0, rowOrder=rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal=paletteContinuous(set="solarExtra"), varCutOff=0, rowOrder=rowOrder)
plotPDF(ht1 + ht2, name="IFE_GSM_MM_corr_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=12, height=7)

# GeneIntegrationMatrix:
corGIM_MM <- correlateTrajectories(IFE.trajGIM, IFE.trajMM)

trajGIM2 <- IFE.trajGIM[corGIM_MM[[1]]$matchname1, ]
trajMM2 <- IFE.trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat=TRUE, varCutOff=0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))

ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal=cmaps_BOR$sunrise,  varCutOff=0, rowOrder=rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal=paletteContinuous(set="solarExtra"), varCutOff=0, rowOrder=rowOrder)
plotPDF(ht1 + ht2, name="IFE_GIM_MM_corr_trajectories.pdf", ArchRProj=atac_proj, 
    addDOC=FALSE, width=12, height=7)

####################################################################################
# Identify most highly-regulated keratins
####################################################################################

# identify 'highly-regulated' genes
p2gGR <- getP2G_GR(atac_proj, corrCutoff=0.45)
p2gFreqs <- getFreqs(p2gGR$symbol)
x <- 1:length(p2gFreqs)
p2g_freq_df <- data.frame(npeaks=p2gFreqs, rank=x)

allGenes <- getGenes(atac_proj)$symbol %>% unname() %>% sort()

# identify groups of cells long trajectory

getQuantileTrajGroups <- function(proj, trajName, groupEvery=10){
    # Return list of groups of cell names along trajectory
    ###############################################################
    # proj = ArchR project
    # trajName = trajectory name (in cellColData)
    # groupEvery = how frequently to bin cells (percentiles)
    trajectory <- getCellColData(proj, trajName)
    trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]
    breaks <- seq(0, 100, groupEvery)
    groupList <- lapply(seq_along(breaks), function(x){
          if(x == 1){
              NULL
          }else{
              rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] & trajectory[,1] <= breaks[x])]
          }
      })[-1]
    names(groupList) <- paste0("Tbin_", breaks[-length(breaks)], "_", breaks[-1])
    groupList
}

# Assign trajectory groups
trajBins <- getQuantileTrajGroups(atac_proj, "IFE.Curve1", groupEvery=20)
binLabels <- invertList(trajBins) %>% unlist()
atac_proj$Tgroup <- ifelse(getCellNames(atac_proj) %in% names(binLabels), binLabels[getCellNames(atac_proj)], "None")

# Get motifs to plot as marker regions
motifPositions <- getPositions(atac_proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

show_motifs <- c("TP63", "KLF4")
motif_names <- unique(motifGR$motifName)
motif_names <- sapply(show_motifs, function(m) motif_names[grepl(m, motif_names)]) %>% as.character()
flist <- lapply(motif_names, function(m){
    motifGR[motifGR$motifName == m] %>% resize(250, fix="center")
    })
names(flist) <- motif_names
flist[["peaks"]] <- allPeaksGR


# (Define plot region based on bracketing linked peaks)
promoterGR <- promoters(getGenes(atac_proj))
plotGenes <- c("ITGA6", "ITGB1", "KRT14", "KRT1", "KRT10", "DMKN") 

mPromoterGR <- promoterGR[promoterGR$symbol %in% plotGenes]
mP2G_GR <- p2gGR[p2gGR$symbol %in% plotGenes]

# Restrict to only loops linking genes of interest
plotLoops <- getPeak2GeneLinks(atac_proj, corCutOff=0.45, resolution = 100)[[1]]
sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
plotLoops <- plotLoops[width(plotLoops) > 100]

# Bracket plot regions around SNPs
plotRegions <- lapply(plotGenes, function(x){
  gr <- range(plotLoops[plotLoops$symbol == x])
  lims <- grLims(gr)
  gr <- GRanges(
      seqnames = seqnames(gr)[1],
      ranges = IRanges(start=lims[1], end=lims[2])
    )
  gr
  }) %>% as(., "GRangesList") %>% unlist()

widths <- width(plotRegions) + 0.05*width(plotRegions)
widths <- ifelse(widths > 100000, widths, 100000)
plotRegions <- resize(plotRegions, width=widths, fix="center")

# Tracks of genes:
p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "Tgroup", 
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    features = flist,
    sizes = c(6, 0.33, 1, 1),
    pal = getColorMap(ArchRPalettes$fireworks2, length(names(trajBins)), type="quantitative"),
    geneSymbol = plotGenes, 
    useGroups = names(trajBins),
    region = plotRegions, 
    loops = plotLoops,
    tileSize=300,
    minCells=100
)

plotPDF(plotList = p, 
    name = "IFE_traj_track_plots_keratinocyteP2G_g20.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 7, height = 6)



##########################################################################################
# Violin plots of (integrated) RNA expression for select genes
##########################################################################################

GImat <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix")
data_mat <- assays(GImat)[[1]]
rownames(data_mat) <- rowData(GImat)$name
sub_mat <- data_mat[plotGenes,]

# These DO NOT match the order of the above matrix by default
grouping_data <- data.frame(tgroup=factor(atac_proj$Tgroup, 
  ordered=TRUE, levels=names(trajBins)))
rownames(grouping_data) <- getCellNames(atac_proj)
grouping_data <- grouping_data[!is.na(grouping_data$tgroup),,drop=FALSE]
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

tgroup_cmap <- getColorMap(ArchRPalettes$fireworks2, length(names(trajBins)), type="quantitative")
names(tgroup_cmap) <- names(trajBins)

pList <- list()
for(gn in plotGenes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per tgroup
  set.seed(1)
  df <- df %>% group_by(tgroup) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$tgroup %in% names(trajBins),]

  covarLabel <- "tgroup"  

  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=tgroup, y=gene, fill=tgroup))
    + geom_violin(aes(fill=tgroup), adjust = 1.5, scale='width', position=dodge)
    + scale_color_manual(values=tgroup_cmap, limits=names(tgroup_cmap), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=tgroup_cmap)
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

pdf(paste0(plotDir, "/Expression_Violin_Tgroup_g20.pdf"), width=10, height=4)
pList
dev.off()


########################################################
# Dot plots of marker genes
########################################################

fineClustOrder <- c(
    # Keratinocytes
    "aKc1", # "Basal.Kc_1",
    "aKc2", # "Spinous.Kc_2",
    "aKc3", # "Spinous.Kc_1",
    "aKc4", # "HF.Kc_1", # SOX9, DKK3
    "aKc5", # "HF.Kc_2", # Lhx2, LGR5 high
    "aKc7", # "HF.Kc_4", # Lhx2, LGR5 high
    "aKc6", # "HF.Kc_3", 
    "aKc8", # "HF.Kc_5", # CD200 high
    "aKc9", # "HF.Kc_6", # Matrix?
    "aKc10" # "Glandular_1"
    #"aKc11" # "Unknown", # Mix of random stuff (doublet?)
)
LFineClustOrder <- unlist(atac.FineClust)[fineClustOrder]

markerGenes  <- c(
  # https://www.proteinatlas.org/humanproteome/tissue/skin
  "KRT14", "KRT15", "COL17A1", # Basal epithelia
  "KRT1", "KRT10", # Spinous epithelia
  #"DSC1", "KRT2", "IVL", # Granular epithelia 
  "ITGA6", "ITGB8", "LGR5","LHX2", "FZDB", "SOX9", # HFSCs/TAC
  "KRT83", "HOXC13", "LEF1", # Matrix hair keratins/genes
  #"KRT71", # IRS keratins / genes
  "KRT75", # ORS keratins / genes
  "ELOVL3", # Sebaceous 
  "CD200", "AR", # CD200 positive cluster?
  #"MKI67", "CDK1", "TOP2A", # Cycling
  "HLA-A", "HLA-B", # MHC1
  "WNT3", "WNT5A", "WNT10A", "WNT10B", # WNTs
  #"FOXC1", # Highly-regulated genes
  "RUNX1", "RUNX2", "RUNX3", # RUNX
  "KRT7", "KRT8", "AQP5" # Glandular
) %>% unique()

# Dot plot of GeneIntegrationMatrix cluster markers
GIM_se <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix")
GIM_mat <- assays(GIM_se)$GeneIntegrationMatrix
rownames(GIM_mat) <- rowData(GIM_se)$name

avgPctMat <- avgAndPctExpressed(GIM_mat[,getCellNames(atac_proj)], atac_proj$FineClust, feature_normalize=TRUE, min_pct=0)

# Subset to genes we care about:
avgPctMat <- avgPctMat[avgPctMat$feature %in% markerGenes,]
avgPctMat <- avgPctMat[avgPctMat$grp %in% fineClustOrder,]

# Assign labels
avgPctMat$grp <- unlist(atac.FineClust)[as.character(avgPctMat$grp)]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")
wide_df <- prettyOrderMat(wide_df[,LFineClustOrder], clusterCols=FALSE)

grp_order <- colnames(wide_df$mat)
gene_order <- rownames(wide_df$mat) #%>% rev() # Reverse this if planning on using plot vertically

pdf(paste0(plotDir, "/GIM_FineClust_markers_dot_plot_keratinocytes.pdf"), width=6, height=9)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", 
  xorder=grp_order, yorder=gene_order, cmap=cmaps_BOR$sunrise, aspectRatio=1.6)
dev.off()


####################################################################
