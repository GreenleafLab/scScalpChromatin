#!/usr/bin/env Rscript

##########################################################################
# Link genes to gkmSVM model prioritized SNPs
##########################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(parallel)
})

# Set Threads to be used
ncores <- 8
addArchRThreads(threads = ncores)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
#source(paste0(scriptPath, "/GO_wrappers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered"
gkm_res_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/snp_results"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

# Color Maps
broadClustCmap <- readRDS(paste0(scriptPath, "/scalpClusterColors.rds")) %>% unlist()
atacNamedClustCmap <- readRDS(paste0(scriptPath, "/scATAC_NamedClust_cmap.rds")) %>% unlist()
sample_cmap <- readRDS(paste0(scriptPath, "/sample_cmap.rds"))
atac_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(atac_proj$Sample2)] %>% unlist()

# Get label cmaps
source(paste0(scriptPath, "/cluster_labels.R"))
broadLabelCmap <- broadClustCmap
names(broadLabelCmap) <- unlist(BroadClust)[names(broadLabelCmap)]
atacLabelClustCmap <- atacNamedClustCmap
names(atacLabelClustCmap) <- unlist(atac.NamedClust)[names(atacNamedClustCmap)]

disease_cmap <- head(cmaps_BOR$stallion,3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

# Correlation cutoff for identifying p2g linkages
corrCutoff <- 0.5

# Retrieve GImat prior to re-assigning p2g links
GImat <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix")

##########################################################################################
# Links fine-mapped SNPs to candidate genes using Peak-to-Gene links
##########################################################################################

fmGR <- readRDS("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/snp_centered_250bp_fastas_v3/250bpSNPCentered.rds")

# Load full project p2g links, plot loops, etc.
full_p2gGR <- readRDS(file="/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered/multilevel_p2gGR.rds") # NOT merged or correlation filtered
plot_loop_list <- readRDS(file="/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered/multilevel_plot_loops.rds")

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

# Get full merged p2g links
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)

pol <- findOverlaps(p2gGR, fmGR, type="any", maxgap=-1L, ignore.strand=TRUE)
expandFMGR <- fmGR[to(pol)]
expandFMGR$linkedGene <- p2gGR[from(pol)]$symbol
expandFMGR$SNP_to_gene <- paste(expandFMGR$linked_SNP, expandFMGR$linkedGene, sep="_")

##########################################################################################
# Read in results of gkmSVM models
##########################################################################################

# Load intermediate results
full_sig_res <- readRDS(paste0(gkm_res_dir, "/significant_hits_table.rds"))
snp_gkmexplain <- readRDS(paste0(gkm_res_dir, "/snp_gkmexplain_results.rds"))

# Filter fine-mapping gene links by those that were in gkmSVM significant hits
keep_dis <- c("Male-pattern baldness", "Balding_Type4", "Eczema", "Hair color")
dis_FM_GR <- expandFMGR[expandFMGR$disease_trait %in% keep_dis]

# Add linked genes to gkmSVM results
snp_to_gene_df <- mcols(expandFMGR) %>% as.data.frame() %>% group_by(linked_SNP) %>% summarize(genes=paste(unique(linkedGene), collapse=";")) %>% as.data.frame()
gene_vec <- snp_to_gene_df$genes
names(gene_vec) <- snp_to_gene_df$linked_SNP

full_sig_res$linkedGenes <- gene_vec[full_sig_res$snp]

##########################################################################################
# Plot seqlogos using original gkmexplain matrix
##########################################################################################

library(ggseqlogo)

plotSNPimportance <- function(snp, snp_table, gkm_explain_output, celltype, gr, indices=101:151){
  # Plot the ref, alt, and delta importance scores for a given SNP in a given cell type
  # snp = the snp to plot
  # snp_table = df of snp hits with the following columns: (region, snp, ref_score, alt_score, score_delta,...)
  # gkm_explain_output = giant list of full gkmexplain output 
  # celltype = celltype to plot
  # gr = genomic range of snps

  # Get required names for accessing data
  snp_info <- snp_table[(snp_table$snp == snp & snp_table$cluster == celltype),]
  ref_group <- paste0(celltype, "-ref_snp_seqs")
  alt_group <- paste0(celltype, "-alt_snp_seqs")
  region <- snp_info$region[1]
  region <- paste0(region, "_", snp, "_125") # Need to adjust this for different relative SNP position

  # Get SNP region
  snp_gr <- gr[gr$linked_SNP == snp] %>% resize(width=50, fix="center")
  true_region <- (snp_gr %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})
  ref_alt <- snp_gr$linked_refalt[1]
  disease <- snp_gr$disease_trait[1]

  # Get info on suspected motif match, linked genes, etc
  linked_genes <- snp_info$linkedGenes[1]
  top_motifs <- snp_info$top_motifs[1]
  fm_prob <- snp_info$fm_probs[1]

  # Get matrices for plotting
  ref_matrix <- gkm_explain_output$gkmexplain_output[[ref_group]]$seq_matrices[[region]][,indices]
  alt_matrix <- gkm_explain_output$gkmexplain_output[[alt_group]]$seq_matrices[[region]][,indices]
  delta_matrix <- ref_matrix - alt_matrix

  # Generate plot for all matrices
  upper_lim <- max(rbind(ref_matrix, alt_matrix, delta_matrix)) * 1.2
  lower_lim <- min(rbind(ref_matrix, alt_matrix, delta_matrix)) * 1.2
  mat_list <- list("ref_allele"=ref_matrix, "alt_allele"=alt_matrix, "delta"=delta_matrix)

  p <- (
    ggseqlogo(mat_list, method='custom', seq_type='dna', ncol=1) 
    + ylab('gkmexplain importance')
    + ylim(lower_lim, upper_lim)
    + ggtitle(paste0(disease, " - ", celltype, " - ", snp, " | ", ref_alt, " - ", true_region, "\n", 
      linked_genes, "\n", top_motifs, "\n", "FM probabiltiy: ", fm_prob))
    )
  p
}

valid_hc_snps <- fmGR[fmGR$disease_trait %in% c("Hair color")]$linked_SNP
valid_aga_snps <- fmGR[fmGR$disease_trait %in% c("Balding_Type4", "Male-pattern baldness")]$linked_SNP
valid_ecz_snps <- fmGR[fmGR$disease_trait %in% c("Eczema")]$linked_SNP

hc_sig_res <- full_sig_res[full_sig_res$snp %in% valid_hc_snps,]
hc_sig_res <- hc_sig_res[order(abs(hc_sig_res$prominence), decreasing=TRUE),]

aga_sig_res <- full_sig_res[full_sig_res$snp %in% valid_aga_snps,]
aga_sig_res <- aga_sig_res[order(abs(aga_sig_res$prominence), decreasing=TRUE),]

ecz_sig_res <- full_sig_res[full_sig_res$snp %in% valid_ecz_snps,]
ecz_sig_res <- ecz_sig_res[order(abs(ecz_sig_res$prominence), decreasing=TRUE),]

plot_hc <- hc_sig_res[!duplicated(hc_sig_res$snp),]
plot_aga <- aga_sig_res[!duplicated(aga_sig_res$snp),]
plot_ecz <- ecz_sig_res[!duplicated(ecz_sig_res$snp),]

# Only plot hits that have a linked gene
plot_hc <- plot_hc[!is.na(plot_hc$linkedGenes),]
plot_aga <- plot_aga[!is.na(plot_aga$linkedGenes),]
plot_ecz <- plot_ecz[!is.na(plot_ecz$linkedGenes),]

# Most prominent, high-effect SNP hits
nplot <- min(nrow(plot_hc), 50)
pdf(paste0(gkm_res_dir, "/most_prominent_hc_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in 1:nplot){
  snp <- plot_hc$snp[i]
  ct <- plot_hc$cluster[i]
  plotList[[i]] <- plotSNPimportance(snp, plot_hc, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()

nplot <- min(nrow(plot_aga), 50)
pdf(paste0(gkm_res_dir, "/most_prominent_aga_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in 1:nplot){
  snp <- plot_aga$snp[i]
  ct <- plot_aga$cluster[i]
  plotList[[i]] <- plotSNPimportance(snp, plot_aga, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()

nplot <- min(nrow(plot_ecz), 50)
pdf(paste0(gkm_res_dir, "/most_prominent_ecz_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in 1:nplot){
  snp <- plot_ecz$snp[i]
  ct <- plot_ecz$cluster[i]
  plotList[[i]] <- plotSNPimportance(snp, plot_ecz, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()


## Candidate SNPs of interest
all_candidate_snps <- c(
  "rs2058622", # IL18RAP linked, eczema
  "rs72966077", # WNT10A linked, AGA
  "rs10769041", # ALX4 linked, AGA
  "rs12350739" # BNC2 linked, HC
)

# Plot all model results for candidate SNPs:
sub_sig_res <- full_sig_res[sapply(full_sig_res$snp, function(x) grepl(paste(all_candidate_snps, collapse="|"), x)),]

pdf(paste0(gkm_res_dir, "/all_candidate_snp_model_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in 1:nrow(sub_sig_res)){
  snp <- sub_sig_res$snp[i]
  ct <- sub_sig_res$cluster[i]
  plotList[[i]] <- plotSNPimportance(snp, sub_sig_res, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()


# Tracks of genes:
disPicsGR <- expandFMGR[expandFMGR$disease_trait %in% keep_dis]
promoterGR <- promoters(getGenes(atac_proj))
candidate_GR <- disPicsGR[disPicsGR$linked_SNP %in% all_candidate_snps]
index_df <- data.frame(index=candidate_GR$index_SNP, candidate=candidate_GR$linked_SNP)
index_df <- index_df[!duplicated(index_df$candidate),]
index_map <- index_df$index
names(index_map) <- index_df$candidate

# The original index SNP is often not in a peak, so in order to find it we need to load the unfiltered fm gr
fm_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gwas/PICS2"
full_fm_gr <- readRDS(paste0(fm_dir, "/unfiltered_finemapping_genomic_range.rds"))
dis_full_fm_gr <- full_fm_gr[full_fm_gr$disease_trait  %in% keep_dis]
index_GR <- dis_full_fm_gr[dis_full_fm_gr$linked_SNP %in% index_map] %>% unique()

# Marker Genes
markerGenes <- candidate_GR$linkedGene %>% unique()

# Combine plot loops from multi-level p2g linking
plotLoops <- unlist(as(plot_loop_list, "GRangesList"))
plotLoops <- plotLoops[order(plotLoops$value, decreasing=TRUE)] %>% unique() %>% sort()

mPromoterGR <- promoterGR[promoterGR$symbol %in% markerGenes]
mP2G_GR <- p2gGR[p2gGR$symbol %in% markerGenes]

flist <- list()
flist[["peaks"]] <- getPeakSet(atac_proj)
flist[["index_SNPs"]] <- unique(index_GR) %>% resize(250, fix="center")
flist[["linked_SNPs"]] <- unique(dis_full_fm_gr[dis_full_fm_gr$index_SNP %in% index_map]) %>% resize(250, fix="center")
flist[["candidate_SNPs"]] <- candidate_GR[!duplicated(candidate_GR)] %>% resize(250, fix="center")

sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops <- plotLoops[width(plotLoops) > 500]

# Bracket plot regions around SNPs
plotRegions <- lapply(markerGenes, function(x){
  candidate_snps <- candidate_GR[candidate_GR$linkedGene == x]$linked_SNP
  gr <- c(
    range(candidate_GR[candidate_GR$linkedGene == x]), 
    range(index_GR[index_GR$linked_SNP %in% index_map[candidate_snps]]), 
    resize(range(mPromoterGR[mPromoterGR$symbol == x]), width=2000))
  lims <- grLims(gr)
  message(sprintf("Trying %s...", x))
  gr <- GRanges(
      seqnames = seqnames(gr)[1],
      ranges = IRanges(start=lims[1], end=lims[2])
    )
  gr
  }) %>% as(., "GRangesList") %>% unlist()
plotRegions <- resize(plotRegions, 
  width=width(plotRegions) + 0.1*width(plotRegions), 
  fix="center")

plotRegions <- resize(plotRegions, width=ifelse(width(plotRegions) > 100000, width(plotRegions), 100000), fix="center")

# Map FineClust colors to NamedClust
cM <- as.matrix(confusionMatrix(atac_proj$FineClust, atac_proj$NamedClust))
allFineClust <- unique(atac_proj$FineClust)
map_colors <- apply(cM, 1, function(x)colnames(cM)[which.max(x)])
atacFineClustCmap <- atacNamedClustCmap[map_colors[allFineClust]]
names(atacFineClustCmap) <- allFineClust


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
    "aFb1", # "D.Fib" 
    "aFb2", # "D.Sheath"
    "aMu1", # "Muscle_1"
    "aMu2", # "Muscle_2" 
    "aVe1", # "Vas.Endo_1"
    "aVe2", # "Vas.Endo_2"
    "aLe1", # "Lymph.Endo"
    "aMe1" # "Melanocytes"
)

# NamedClust
p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "NamedClust",
    useGroups = atacOrder,
    features = flist,
    pal = atacNamedClustCmap,
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    sizes = c(8, 0.8, 1.25, 0.5),
    region = plotRegions, 
    loops = getPeak2GeneLinks(atac_proj, corCutOff=corrCutoff), # All peak-to-gene loops
    tileSize=250,
    minCells=100,
    title = markerGenes
)

plotPDF(plotList = p, 
    name = "aga_candidate_snp_tracks_100k_width_allClust_allLoops.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 8, height = 7)


##################################################
# Violin plots of RNA expression for select genes
##################################################

markerGenes <- c("IL18RAP", "WNT10A", "ALX4", "BNC2")

data_mat <- assays(GImat)[[1]]
rownames(data_mat) <- rowData(GImat)$name
sub_mat <- data_mat[markerGenes,]

grouping_data <- data.frame(cluster=factor(atac_proj$NamedClust, 
  ordered=TRUE, levels=atacOrder))
rownames(grouping_data) <- getCellNames(atac_proj)
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

pList <- list()
for(gn in markerGenes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  set.seed(1)
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$cluster %in% atacOrder,] %>% as.data.frame()

  covarLabel <- "cluster"  

  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=gene, fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
    + scale_color_manual(values=atacNamedClustCmap, limits=names(atacNamedClustCmap), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=atacNamedClustCmap)
    + guides(fill=guide_legend(title=covarLabel), 
      colour=guide_legend(override.aes = list(size=5)))
    + ggtitle(gn)
    + xlab("")
    + ylab("Integrated RNA Expression")
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Expression_Violin_gkmsvm_model_genes_NamedClust_unclipped.pdf"), width=10, height=4)
pList
dev.off()
