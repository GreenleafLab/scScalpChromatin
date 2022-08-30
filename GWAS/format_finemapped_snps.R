#!/usr/bin/env Rscript

##########################################################################
# Format finemapped SNPs for later analyses
##########################################################################


#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(parallel)
  library(rtracklayer)
})

# Set Threads to be used
ncores <- 8
addArchRThreads(threads = ncores)

# Get additional functions, etc.:
scriptPath <- "/home/users/boberrey/git_clones/scScalpChromatin/"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered"
pics_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gwas/PICS2"
fin_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gwas/finucane_finemapping/release1.1"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

atac_proj <- loadArchRProject(wd, force=TRUE)
peaks_gr <- getPeakSet(atac_proj)
# Peaks do not have unique names by default...
peaks_gr$peakName <- (peaks_gr %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})

##########################################################################################
# Preparing Data
##########################################################################################

# Parser for PICS2 output files

parsePICSfile <- function(pics_file, extra_cols=NULL, check_header=TRUE){
  # Read in a PICS2 output file and return a genomic range
  ######################################################
  # pics_file: path/to/PICS/file.txt (should be a tab-delimited text file)
  # extra_cols: extra columns we want to keep in returned GR
  message(sprintf("Reading in file %s...", pics_file))
  dt <- fread(pics_file, fill=TRUE, na.strings="")

  # Some PICS output files have an annoying quirk where they repeat the header for each index SNP
  
  header <- colnames(dt)
  if(check_header){
    message("Removing intra-table header rows...")
    dt <- dt[!apply(dt, 1, function(x) identical(as.vector(x), header)),]
  }

  # Match similar column names for the ones we want to keep
  message("Relabeling column names...")
  mlist <- list(
      "index_SNP" = c("index_SNP", "indexsnp", "Index_SNP", "IndexSNP"),
      "linked_SNP" = c("linked_SNP", "linked_snp", "Linked_SNP", "LinkedSNP"),
      "linked_pos" = c("linked_pos", "Linked_position_hg38 (hg19)", "linked_position"),
      "rsq" = c("rsq", "Rsquare", "rsquare"),
      "dprime" = c("dprime", "Dprime"),
      "linked_refalt" = c("linked_refalt", "Linked_RefAlt"),
      "fm_prob" = c("PICS_prob", "PICS Prob", "pics_probability", "PICS_probability"),
      "consequence" = c("consequence", "consequences", "Consequence"),
      "phase" = c("phase", "Phase"),
      "disease_trait" = c("disease_trait", "disease", "trait", "Disease/Trait")
    )
  imlist <- invertList(mlist)
  setnames(dt, names(imlist), as.vector(unlist(imlist)), skip_absent=TRUE)

  # Now extract just the columns we want
  keep_cols <- unique(c(names(mlist),extra_cols))
  keep_cols <- keep_cols[keep_cols %in% colnames(dt)]
  dt <- dt[, ..keep_cols]

  # Keep only hg38 linked_pos
  message("Preparing genomicRange object...")
  dt$linked_pos <- sapply(strsplit(dt$linked_pos, split=" "), `[`, 1) # (This gets the first element of each split)
  dt$chr <- paste0("chr", sapply(strsplit(dt$linked_pos, split=":"), `[`, 1)) # Add chr to chromosome numbers
  dt$start <- as.integer(sapply(strsplit(dt$linked_pos, split=":"), `[`, 2))
  dt <- dt[,linked_pos:=NULL] # Drop unnecessary column now
  dt$fm_prob <- as.numeric(dt$fm_prob)

  # Rarely, there will not be a position for a SNP (e.g. if SNP was not in PICS database)
  dt <- dt[!is.na(dt$start),]
  
  # Construct Genomic Range
  gr <- makeGRangesFromDataFrame(dt, keep.extra.columns=TRUE, ignore.strand=TRUE, start.field="start", end.field="start") %>% sort()
  message("Done")
  gr
}

# All autoimmune PICS (this dataset doesn't have the linked_refalt column available)
# We can get around this for most by looking up the linked_SNP in the GWAS catalog set and copying the ref:alt there
pics_ai_gr <- parsePICSfile(paste0(pics_dir, "/PICS2-AIsnplist-8mar19-annotated.txt"))
pics_ai_gr$source <- "pics_autoimmune"
# GWAS catalog PICS
gwas_cat_gr <- parsePICSfile(paste0(pics_dir, "/PICS2-GWAScat-2021-06-11.txt.gz"), check_header=FALSE)
gwas_cat_gr$source <- "pics_GWC"

snpLookup <- as.character(gwas_cat_gr$linked_refalt)
names(snpLookup) <- as.character(gwas_cat_gr$linked_SNP)
pics_ai_gr$linked_refalt <- snpLookup[as.character(pics_ai_gr$linked_SNP)]
pics_ai_gr <- pics_ai_gr[!is.na(pics_ai_gr$linked_refalt)]

# Combine PICS results and then filter against peakset
full_pics_gr <- suppressWarnings(
  c(
    pics_ai_gr, # Pre-computed AA dataset
    gwas_cat_gr # full GWAS catalog PICS2 dataset
    )
  ) %>% sort()
rm(gwas_cat_gr); gc()

# Parse Finucane finemapped SNPs from 94 UKBB traits
fin_dt <- fread(paste0(fin_dir, "/UKBB_94traits_release1.bed.gz"))
colnames(fin_dt) <- c("chr_hg19", "start_hg19", "end_hg19", "varID", "linked_SNP", # (linked_SNP = rsid)
  "ref_allele", "alt_allele", "minor_allele", "cohort", "model_marginal", "fm_method", 
  "disease_trait", "region", "maf", "beta_marginal", "se_marginal", "chisq_marginal", "fm_prob", "cs_id", # (fm_prob = pip)
  "beta_post", "sd_post", "LD_HWE", "LD_SV")

# Keep only subset of columns and reformat to match syntax where possible
fin_dt$linked_refalt <- paste0(fin_dt$ref_allele, ":", fin_dt$alt_allele)
keep_cols <- c("chr_hg19", "start_hg19", "end_hg19", "linked_SNP", "linked_refalt", 
  "model_marginal", "fm_method", "disease_trait", "fm_prob")
fin_dt <- fin_dt[,..keep_cols]
fin_dt$source <- "finacune_UKBB"

# Convert to genomic range
fin_gr <- makeGRangesFromDataFrame(fin_dt, keep.extra.columns=TRUE, ignore.strand=TRUE, 
  seqnames.field="chr_hg19", start.field="start_hg19", end.field="end_hg19") %>% sort()

# LiftOver SNPs from hg19 to hg38
chain <- import.chain("/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/liftover/hg19ToHg38.over.chain")
fin_hg38_gr <- liftOver(fin_gr, chain) %>% unlist()

# Resize to match format of PICS GR
fin_hg38_gr <- resize(fin_hg38_gr, width=1, fix="end")

# Combine finemapping GRs
full_finemapping_gr <- c(full_pics_gr, fin_hg38_gr)

# Save a completely unfiltered GR
saveRDS(full_finemapping_gr, file = paste0(pics_dir, "/unfiltered_finemapping_genomic_range.rds"))
full_finemapping_gr <- subsetByOverlaps(full_finemapping_gr, peaks_gr, type="any", maxgap=-1L, ignore.strand=TRUE)

# In the filtered GR, keep only PICS SNPs that pass PICS probability cutoff or have rsq > 0.8
filt_finemapping_gr <- full_finemapping_gr[(full_finemapping_gr$fm_prob >= 0.01) | 
    ((full_finemapping_gr$rsq >= 0.8) & !is.na(full_finemapping_gr$rsq)) | 
    (full_finemapping_gr$source == "finacune_UKBB")]

# Save finemapped results
saveRDS(filt_finemapping_gr, file = paste0(pics_dir, "/filtered_finemapping_genomic_range.rds"))

###################################################################