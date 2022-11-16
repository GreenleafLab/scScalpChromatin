#!/usr/bin/env Rscript

##########################################################################################
# TF-regulator analysis Keratinocytes
##########################################################################################

#Load ArchR (and associated libraries)
library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggrastr)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory
subgroup <- "Keratinocytes"
wd <- sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/subclustered_%s", subgroup)


#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 1.0
corrCutoff <- 0.45 # Used in labeling peak2gene links

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)

# Non-ArchR plots:
plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

#############################################################
# Identify Correlated TF Motifs and TF Gene Score/Expression
#############################################################

# To identify 'Positive TF regulators', i.e. TFs whose motif accessibility 
# is correlated with with their own gene activity (either by gene score or gene expression)

seGroupMotif <- getGroupSE(ArchRProj=atac_proj, useMatrix="MotifMatrix", groupBy="FineClust")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",] # Subset to just deviation z-scores

# identify the maximum delta in z-score between all clusters
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


corGSM_MM <- correlateMatrices(
    ArchRProj = atac_proj,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)

corGIM_MM <- correlateMatrices(
    ArchRProj = atac_proj,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)

# For each correlation analyses, we annotate each motif with the maximum delta observed between clusters
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

# "we consider positive regulators as those TFs whose correlation between motif and gene score 
# (or gene expression) is greater than 0.5 with an adjusted p-value less than 0.01 and a maximum 
# inter-cluster difference in deviation z-score that is in the top quartile (Max TF Motif Delta)."
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
#corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
corGSM_MM$TFRegulator[which(abs(corGSM_MM$cor) > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])


p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  ggrepel::geom_label_repel(
          data = data.frame(corGSM_MM[corGSM_MM$TFRegulator=="YES",]), aes(x = cor, y = maxDelta, label = GeneScoreMatrix_name), 
          size = 1.5,
          #nudge_x = 2,
          color = "black") +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
)

pdf(paste0(plotDir, "/corGSM_MM_posTFregulators.pdf"), width=5, height=5)
p
dev.off()

# Same thing for RNA:
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
#corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
corGIM_MM$TFRegulator[which(abs(corGIM_MM$cor) > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"

sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])


p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  ggrepel::geom_label_repel(
          data = data.frame(corGIM_MM[corGIM_MM$TFRegulator=="YES",]), aes(x = cor, y = maxDelta, label = GeneIntegrationMatrix_name), 
          size = 1.5,
          #nudge_x = 2,
          color = "black") +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
)

pdf(paste0(plotDir, "/corGIM_MM_posTFregulators.pdf"), width=5, height=5)
p
dev.off()


#############################################################
# Identify regulatory targets of TFs 
#############################################################

# ChromVAR deviations matrix: (rows motif names x cols cell names)
motifMatrix <- getMatrixFromProject(atac_proj, useMatrix="MotifMatrix")
deviationsMatrix <- assays(motifMatrix)$deviations

# GeneIntegration Matrix: (rows gene names x cols cell names)
GIMatrix <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix")
GImat <- assays(GIMatrix)$GeneIntegrationMatrix
rownames(GImat) <- rowData(GIMatrix)$name
GImat <- as(GImat[Matrix::rowSums(GImat) > 0,], "sparseMatrix") # Remove unexpressed genes

# Use only motifs that are 'TFRegulators' as determined by analysis above
GSMreg <- rownames(motifMatrix)[corGSM_MM[corGSM_MM$TFRegulator == "YES",]$MotifMatrix_idx]
GIMreg <- rownames(motifMatrix)[corGIM_MM[corGIM_MM$TFRegulator == "YES",]$MotifMatrix_idx]
regulators <- unique(c(GSMreg, GIMreg))
regulators <- c(
  "POU2F3_613", "RORA_658", "REL_721", "CEBPD_152", "CEBPA_155", "KLF10_826", 
  "EGR1_195", "EGR2_196", "EGR3_259", "EGR4_207",
  regulators
  )
deviationsMatrix <- deviationsMatrix[regulators,]

# Identify pseudobulks for performing matrix correlations
knn_groups <- getLowOverlapAggregates(atac_proj, target.agg=500, k=100, 
  overlapCutoff=0.8, dimReduc="Harmony")

kgrps <- unique(knn_groups$group)

# GeneIntegrationMatrix
GIMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(GImat[,use_cells])
  }) %>% do.call(cbind,.)
colnames(GIMatPsB) <- kgrps

# In rare instances, we can get pseudo-bulked genes that have zero averages
GIMatPsB <- GIMatPsB[Matrix::rowSums(GIMatPsB) > 0,]

# DeviationsMatrix
DevMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(deviationsMatrix[,use_cells])
  }) %>% do.call(cbind,.)
colnames(DevMatPsB) <- kgrps

# Perform chromVAR deviations to Integrated RNA correlation analysis:
start <- Sys.time()
geneCorMat <- cor2Matrices(DevMatPsB, GIMatPsB)
colnames(geneCorMat) <- c("motifName", "symbol", "Correlation", "FDR")
end <- Sys.time()
message(sprintf("Finished correlations in %s minutes.", round((end  - start)/60.0, 2)))

allGenes <- rownames(GIMatPsB) %>% sort() # Already filtered to only expressed genes

# Get locations of motifs of interest:
motifPositions <- getPositions(atac_proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

# Get peak to gene GR
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)


calculateLinkageScore <- function(motifLocs, p2gGR){
  # Calculate Linkage Score (LS) for each gene in p2gGR with regards to a motif location GR
  ###################################
  # For a given gene, the LS = sum(corr peak R2 * motifScore)
  ol <- findOverlaps(motifLocs, p2gGR, maxgap=0, type=c("any"), ignore.strand=TRUE)
  olGenes <- p2gGR[to(ol)]
  olGenes$motifScore <- motifLocs[from(ol)]$score
  olGenes$R2 <- olGenes$Correlation**2 # All p2g links here are already filtered to only be positively correlated
  LSdf <- mcols(olGenes) %>% as.data.frame() %>% group_by(symbol) %>% summarise(LS=sum(R2*motifScore)) %>% as.data.frame()
  LSdf <- LSdf[order(LSdf$LS, decreasing=TRUE),]
  LSdf$rank <- 1:nrow(LSdf)
  return(LSdf)
}

calculateMotifEnrichment <- function(motifLocs, p2gGR){
  # Calculate Motif enrichment per gene
  ###################################
  # For a given gene, calculate the hypergeometric enrichment of motifs in 
  # linked peaks (generally will be underpowered)
  motifP2G <- p2gGR[overlapsAny(p2gGR, motifLocs, maxgap=0, type=c("any"), ignore.strand=TRUE)]
  m <- length(motifP2G) # Number of possible successes in background
  n <- length(p2gGR) - m # Number of non-successes in background

  motifLinks <- motifP2G$symbol %>% getFreqs()
  allLinks <- p2gGR$symbol %>% getFreqs()
  df <- data.frame(allLinks, motifLinks=motifLinks[names(allLinks)])
  df$motifLinks[is.na(df$motifLinks)] <- 0
  df$mLog10pval <- apply(df, 1, function(x) -phyper(x[2]-1, m, n, x[1], lower.tail=FALSE, log.p=TRUE)/log(10))
  df <- df[order(df$mLog10pval, decreasing=TRUE),]
  df$symbol <- rownames(df)
  return(df)
}

# plot all TF regulators
regPlotDir <- paste0(plotDir, "/TFregulatorPlots")
dir.create(regPlotDir, showWarnings = FALSE, recursive = TRUE)

markerGenes  <- c(
  # https://www.proteinatlas.org/humanproteome/tissue/skin
  "KRT14", "KRT5", "KRT15", "COL17A1", # Basal epithelia
  "KRT1", "KRT10", # Spinous epithelia
  "DSC1", "KRT2", "IVL", # Granular epithelia 
  "ITGA6", "ITGB1", "CD200", "LGR5","LHX2", "FRZB", "FZD1", "FZD5", "FZD10", # HFSCs
  "CD34", "CDH3", "LGR5", "RUNX1",
  "KRT81", "KRT83", "HOXC13", "LEF1", # Matrix hair keratins/genes
  "KRT75", # ORS keratins / genes
  "SOX9", "LHX2", "NFATC1", "TCF3", # Key HFSC TFs
  "WNT3", "WNT5A", "WNT6", "WNT10A", "WNT10B", # WNTs
  "PERP", "DUSP6", "S100A2", "IL1B"
) %>% unique()

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

markerGenes <- c(cornification_genes, hemidesmosome_genes, hfdev_genes, markerGenes) %>% unique() %>% sort()

# Store results for each TF
res_list <- list()

###########################################
for(motif in regulators){
  motif_short <- strsplit(motif,"_")[[1]][1]
  # First get motif positions
  motifLocs <- motifGR[motifGR$motifName == motif]
  # Calculate Linkage Score for motif
  LS <- calculateLinkageScore(motifLocs, p2gGR)
  # Get just genes correlated to motif
  motifGeneCorDF <- geneCorMat[geneCorMat$motifName == motif,]
  plot_df <- merge(LS, motifGeneCorDF, by="symbol", all.x=TRUE)
  # Calculate motif enrichment per gene
  ME <- calculateMotifEnrichment(motifLocs, p2gGR)
  plot_df <- merge(plot_df, ME, by="symbol", all.x=TRUE)
  plot_df <- plot_df[,c("symbol", "LS", "Correlation", "FDR", "mLog10pval")]
  plot_df$toLabel <- "NO"
  topN <- 5
  plot_df <- plot_df[order(plot_df$LS, decreasing=TRUE),]
  plot_df$rank_LS <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$Correlation, decreasing=TRUE),]
  plot_df$rank_Corr <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=TRUE),]
  plot_df$rank_Pval <- 1:nrow(plot_df)
  plot_df$toLabel[1:10] <- "YES"
  plot_df$meanRank <- apply(plot_df[,c("rank_LS", "rank_Corr", "rank_Pval")], 1, mean)
  plot_df <- plot_df[order(plot_df$meanRank, decreasing=FALSE),]
  plot_df$toLabel[1:topN] <- "YES"
  # Label any marker genes in window of interest
  LS_window <- quantile(plot_df$LS, 0.8)
  corr_window <- 0.25
  pos_top_genes <- plot_df[plot_df$LS > LS_window & plot_df$Correlation > corr_window,]$symbol
  neg_top_genes <- plot_df[plot_df$LS > LS_window & -plot_df$Correlation > corr_window,]$symbol
  if(nrow(plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]) > 0){
    plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]$toLabel <- "YES"
  }
  res_list[[motif_short]] <- pos_top_genes # Save regulatory targets
  # Save dataframe of results
  save_df <- plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes),c(1:5)]
  save_df <- save_df[order(save_df$Correlation, decreasing=TRUE),]
  saveRDS(save_df, paste0(regPlotDir, sprintf("/regulatory_targets_%s.rds", motif_short)))
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=FALSE),]
  # Label motif as well
  plot_df$toLabel[which(plot_df$symbol == motif_short)] <- "YES"
  plot_df$symbol[which(plot_df$toLabel == "NO")] <- ""
  # Threshold pvalue for plotting
  maxPval <- 5
  plot_df$mLog10pval <- ifelse(plot_df$mLog10pval > maxPval, maxPval, plot_df$mLog10pval)
  #Plot results
  p <- (
    ggplot(plot_df, aes(x=Correlation, y=LS, color=mLog10pval)) 
      #+ geom_point(size = 2)
      + ggrastr::geom_point_rast(size=2)
      + ggrepel::geom_text_repel(
          data=plot_df[plot_df$toLabel=="YES",], aes(x=Correlation, y=LS, label=symbol), 
          #data = plot_df, aes(x=Correlation, y=LS, label=symbol), #(don't do this, or the file will still be huge...)
          size=2,
          point.padding=0, # additional pading around each point
          box.padding=0.5,
          min.segment.length=0, # draw all line segments
          max.overlaps=Inf, # draw all labels
          #nudge_x = 2,
          color="black"
      ) 
      + geom_vline(xintercept=0, lty="dashed") 
      + geom_vline(xintercept=corr_window, lty="dashed", color="red")
      + geom_vline(xintercept=-corr_window, lty="dashed", color="red")
      + geom_hline(yintercept=LS_window, lty="dashed", color="red")
      + theme_BOR(border=FALSE)
      + theme(panel.grid.major=element_blank(), 
              panel.grid.minor= element_blank(), 
              plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
              aspect.ratio=1.0,
              #legend.position = "none", # Remove legend
              axis.text.x = element_text(angle=90, hjust=1))
      + ylab("Linkage Score") 
      + xlab("Motif Correlation to Gene") 
      + scale_color_gradientn(colors=cmaps_BOR$zissou, limits=c(0, maxPval))
      + scale_y_continuous(expand = expansion(mult=c(0,0.05)))
      + ggtitle(sprintf("%s putative targets", motif_short))
      )
  # Positively regulated genes:
  upGO <- rbind(
    calcTopGo(allGenes, interestingGenes=pos_top_genes, nodeSize=5, ontology="BP"), 
    calcTopGo(allGenes, interestingGenes=pos_top_genes, nodeSize=5, ontology="MF")
    )
  upGO <- upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
  up_go_plot <- topGObarPlot(upGO, cmap=cmaps_BOR$comet, nterms=6, border_color="black", 
    barwidth=0.9, title=sprintf("%s putative targets (%s genes)", motif_short, length(pos_top_genes)), enrichLimits=c(0, 6))
  # Negatively regulated genes:
  downGO <- rbind(
    calcTopGo(allGenes, interestingGenes=neg_top_genes, nodeSize=5, ontology="BP"), 
    calcTopGo(allGenes, interestingGenes=neg_top_genes, nodeSize=5, ontology="MF")
    )
  downGO <- downGO[order(as.numeric(downGO$pvalue), decreasing=FALSE),]
  down_go_plot <- topGObarPlot(downGO, cmap=cmaps_BOR$comet, nterms=6, border_color="black", 
    barwidth=0.9, title=sprintf("%s putative targets (%s genes)", motif_short, length(neg_top_genes)), enrichLimits=c(0, 6))
  pdf(paste0(regPlotDir, sprintf("/%s_LS.pdf", motif_short)), width=8, height=6)
  print(p)
  print(up_go_plot)
  print(down_go_plot)
  dev.off()
}
###########################################

#################################################################################################################################
# Test for enrichment of differential genes from orthogonal datasets of mutant TFs or TF knockdowns
#################################################################################################################################

### Qu et al (p63 mutant keratinocytes RNA-seq) ###
qu_dt <- fread("/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gene_sets/Qu_2018_SuppTable1d.csv")
colnames(qu_dt) <- c("ENSEMBL_ID", "gene_name", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj")
qu_dt <- qu_dt[qu_dt$gene_name %in% allGenes,] # Filter genes that are not expressed in our dataset
qu_neg <- qu_dt[qu_dt$log2FC < 0]$gene_name
qu_pos <- qu_dt[qu_dt$log2FC > 0]$gene_name

### Fortunel et al (KLF4 knockdown keratinocytes RNA-seq) ###

# Need to re-run DESeq2:

# Load unnormalized counts table
fortunel_raw_dt <- fread("/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gene_sets/GSE111786_counts_raw.csv.gz")

count_mat <- as.matrix(fortunel_raw_dt[,2:ncol(fortunel_raw_dt)])
rownames(count_mat) <- fortunel_raw_dt %>% pull(1)
colnames(count_mat) <- c(
  "KLFWT_1", "KLFWT_2", "KLFWT_3",
  "KLFWTGFP_1", "KLFWTGFP_2", "KLFWTGFP_3",
  "KLFKDGFP_1", "KLFKDGFP_2", "KLFKDGFP_3"
)

# Perform light pre-filtering
count_mat <- count_mat[rowSums(count_mat) > 10,]

# Assign conditions
samples <- data.frame(condition=factor(gsub('_[[:digit:]]+', '', colnames(count_mat))))
rownames(samples) <- colnames(count_mat)

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=samples, design= ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "KLFKDGFP", "KLFWTGFP"))
res <- res[order(res$padj, decreasing=FALSE),]

# Translate ensemble IDs to symbols
library(org.Hs.eg.db)
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl_to_symbol <- select(
  org.Hs.eg.db, 
  keys=unique(rownames(res)), 
  columns=c("ENSEMBL", "SYMBOL"), 
  keytype="ENSEMBL"
  )
ens_to_sym <- ensembl_to_symbol$SYMBOL
names(ens_to_sym) <- ensembl_to_symbol$ENSEMBL
res$symbol <- ens_to_sym[rownames(res)]
res <- res[!is.na(res$symbol),]
res <- res[!is.na(res$padj),]
sig_res <- res[res$padj <= 0.05,]
sig_res <- sig_res[sig_res$symbol %in% allGenes,] # Filter genes that are not expressed in our dataset
fortunel_pos <- sig_res[sig_res$log2FoldChange > 0,]$symbol %>% unname()
fortunel_neg <- sig_res[sig_res$log2FoldChange < 0,]$symbol %>% unname()

### Test for enrichment in regulatory targets
reg_targets <- list()
reg_targets[["TP63"]] <- res_list$TP63
reg_targets[["KLF4"]] <- res_list$KLF4

diff_genes <- list()
diff_genes[["TP63_up"]] <- qu_pos
diff_genes[["TP63_down"]] <- qu_neg
diff_genes[["KLF4_up"]] <- fortunel_pos
diff_genes[["KLF4_down"]] <- fortunel_neg

enrichment_res <- lapply(names(reg_targets), function(tf){
  tf_genes <- reg_targets[[tf]]
  ntf_genes <- allGenes[allGenes %ni% tf_genes]
  diff_res <- lapply(names(diff_genes), function(dg){
    dgenes <- diff_genes[[dg]] # differential genes
    tfnol <- tf_genes[tf_genes %in% dgenes] %>% length() # n TF genes that are also diff genes
    tfnnol <- length(tf_genes) - tfnol # n TF genes that are not diff genes
    dgnol <- length(dgenes) - tfnol # n non-TF genes that are also diff genes
    dgnnol <- length(ntf_genes) - dgnol # n non-TF genes that are also not diff genes
    OR <- (tfnol/tfnnol)/(dgnol/dgnnol)
    pval <- fisher.test(matrix(c(tfnol, tfnnol, dgnol, dgnnol),2,2), alternative="greater")$p.value
    list(OR=OR, pval=pval)
    }) %>% do.call(rbind,.) %>% as.data.frame()
  rownames(diff_res) <- names(diff_genes)
  diff_res
  }) 
names(enrichment_res) <- names(reg_targets)


