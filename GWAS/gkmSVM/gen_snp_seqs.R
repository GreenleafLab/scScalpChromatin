#!/usr/bin/env Rscript

##########################################################################
# Create fasta files containing SNP-modified peak sequences
##########################################################################


#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
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
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/scATAC_preprocessing/fine_clustered"
fm_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gwas/PICS2"

#Set/Create Working Directory to Folder
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

atac_proj <- loadArchRProject(wd, force=TRUE)

peaks_gr <- getPeakSet(atac_proj)
# Peaks do not have unique names by default
peaks_gr$peakName <- (peaks_gr %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})

##########################################################################################
# Preparing Data
##########################################################################################

# Load previously generated fine-mapping GR
full_fm_gr <- readRDS(paste0(fm_dir, "/unfiltered_finemapping_genomic_range.rds"))

# Some of the fine-mapped SNPs are duplicated (i.e. the Finacune SNPs sometimes have both FINEMAP and SuSiE finemapping results)
# Deduplicate trait-SNP pairs prior to proceeding with enrichment analyses:
full_fm_gr <- full_fm_gr[order(full_fm_gr$fm_prob, decreasing=TRUE)]
full_fm_gr$trait_snp <- paste0(full_fm_gr$disease_trait, "_", full_fm_gr$linked_SNP)
full_fm_gr <- full_fm_gr[!duplicated(full_fm_gr$trait_snp)] %>% sort()

# Subset to traits of interest
pos_traits <- c(
  "Male-pattern baldness", # PICS Finemapping (29950 total SNPs)
  "Balding_Type4", # Finacune Finemapping (21948 total SNPs)
  "Eczema", # PICS Finemapping (7898 total SNPs)
  "Hair color" # PICS Finemapping (15630 total SNPs)
)

pos_fm_gr <- full_fm_gr[full_fm_gr$disease_trait %in% pos_traits]

# Sort pos traits by fine-mapping posterior prob and deduplicate
pos_fm_gr <- pos_fm_gr[order(pos_fm_gr$fm_prob, decreasing=TRUE)]
pos_fm_gr <- pos_fm_gr[!duplicated(pos_fm_gr$linked_SNP)]
# (don't use any overlapping SNPs)
neg_fm_gr <- full_fm_gr[full_fm_gr$linked_SNP %ni% pos_fm_gr$linked_SNP]
neg_fm_gr <- neg_fm_gr[neg_fm_gr$source == "pics_GWC"] # Some SNPs from the Finacune finemapping data have the 'wrong' ref:alt allele
neg_fm_gr <- neg_fm_gr[order(neg_fm_gr$fm_prob, decreasing=TRUE)]
neg_fm_gr <- neg_fm_gr[!duplicated(neg_fm_gr$linked_SNP)]
neg_fm_gr$disease_trait <- "random" # Negative traits are random PICS fine-mapped SNPs from the GWAS catalog

# Keep only positive SNPs that pass fine-mapping probability cutoff
pos_fm_gr <- pos_fm_gr[pos_fm_gr$fm_prob >= 0.01]

# Recombine groups
full_fm_gr <- c(pos_fm_gr, neg_fm_gr) %>% sort()

# Male-pattern baldness:  6732
# Balding_Type4:          4473
# Eczema:                 1890
# Hair color:             3798
# Random:                 1983172

# Get ref length and set actual size
reflens <- sapply(strsplit(full_fm_gr$linked_refalt, split=":"), `[`, 1) %>% nchar() %>% as.numeric()
altlens <- sapply(strsplit(full_fm_gr$linked_refalt, split=":"), `[`, 2) %>% nchar() %>% as.numeric()

# Keep only true SNPs (i.e. single base swaps)
full_fm_gr <- full_fm_gr[(reflens == 1 & altlens == 1)]

# Male-pattern baldness:  6133
# Balding_Type4:          3908
# Eczema:                 1732
# Hair color:             3493
# Random:                 1771521

# Subset to only snps fully contained within peaks
full_fm_gr <- subsetByOverlaps(full_fm_gr, peaks_gr, type="within", maxgap=-1L, ignore.strand=TRUE)

# Male-pattern baldness:  946
# Balding_Type4:          685
# Eczema:                 365
# Hair color:             612
# Random:                 189076

# Keep only 2500 of negative SNP set
ngr <- full_fm_gr[full_fm_gr$disease_trait == "random"]
full_fm_gr <- full_fm_gr[full_fm_gr$disease_trait != "random"]
set.seed(1)
ngr <- ngr[sample(1:length(ngr), size=2500, replace=FALSE)]
full_fm_gr <- c(full_fm_gr, ngr) %>% sort()

# Save SNP start position
full_fm_gr$snp_start <- start(full_fm_gr)

# identify containing peak and save SNP position relative to peak start
ol <- findOverlaps(full_fm_gr, peaks_gr, type="within", maxgap=0, ignore.strand=TRUE)
full_fm_gr$peak <- peaks_gr[to(ol)]$peakName
full_fm_gr$snp_rel_pos <- full_fm_gr$snp_start - start(peaks_gr[to(ol)])

# SNP rs28971643 (Balding_Type4) does not have the correct ref/alt allele.
full_fm_gr <- full_fm_gr[full_fm_gr$linked_SNP != "rs28971643"]

# 250bp window centered on SNP position
snps_250bp_center_gr <- resize(full_fm_gr, width=251, fix="center")
names(snps_250bp_center_gr) <- (snps_250bp_center_gr %>% {paste0(seqnames(.), "_", start(.), "_", end(.), "_", .$linked_SNP)})
snps_250bp_center_gr$refseq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, names=snps_250bp_center_gr)
snps_250bp_center_gr$snp_rel_pos <- snps_250bp_center_gr$snp_start - start(snps_250bp_center_gr)

# Prepare list of genomic ranges to run
gr_list <- list(
  "250bpSNPCentered" = snps_250bp_center_gr
)

########################################
# Functions
########################################

shuffleSNPseqs <- function(gr, nseqs=10){
  # Shuffle the dinucleotide content of a sequence while maintaining the original
  # snp position and allele
  # gr = genomic range with the following mcols:
  #   refseq = sequence to be shuffled
  #   linked_refalt = the ref and alt alleles separated by a ':'
  #   snp_start = genomic location of snp start

  # First get all ref alleles
  refs <- sapply(strsplit(gr$linked_refalt, split=":"), `[`, 1)
  ref_start <- gr$snp_start - start(gr) + 1
  ref_end <- ref_start + unname(nchar(refs)) - 1

  # Get sequences with ref allele removed
  message("Generating sequences without reference allele...")
  no_ref_seqs <- sapply(1:length(gr$refseq), function(i){
    replaceAt(gr$refseq[i], IRanges(ref_start[i], ref_end[i]), "")
    }) %>% do.call(c,.)

  # Get shuffled sequences
  shuffled_seqs <- getDinucShuffledSeqs(no_ref_seqs, nseqs=nseqs)

  # Replace SNP allele
  expanded_starts <- rep(ref_start, each=nseqs)
  expanded_refs <- rep(refs, each=nseqs)
  expanded_gr <- rep(gr, each=nseqs)

  message("Reinserting reference allele into shuffled sequences...")
  expanded_gr$refseq <- sapply(1:length(shuffled_seqs), function(i){
    replaceAt(shuffled_seqs[i], expanded_starts[i], expanded_refs[i]) # Insert allele back in
    }) %>% do.call(c,.)
  return(expanded_gr)
}


getDinucShuffledSeqs <- function(seqs, nseqs=1, meme_path="/home/users/boberrey/software/meme-5.4.1/scripts"){
  # Create shuffled sequences using MEME suite's dinucleotide shuffle
  # seqs = vector of sequences to shuffle
  # nseqs = how many shuffled sequences to generate for each sequence in seqs
  # meme_path = path to meme fasta-dinucleotide-shuffle script
  #
  # Use MEME suite's dinucleotide shuffle:
  # USAGE:
  #     software/meme/lib/meme-5.4.1/python/fasta-dinucleotide-shuffle.py [options]
  #     -f <filename>   file name (required)
  #     -t <tag>        added to shuffled sequence names
  #     -s <seed>       random seed; default: 1
  #     -c <n>          make <n> shuffled copies of each sequence; default: 1
  #     -a <filename>   alphabet file to use non-DNA alphabets
  #     -h              print this usage message
  #
  # First save sequences to temporary fasta
  in_file <- "./_temp_meme_input.fasta"
  out_file <- "./_temp_meme_output.fasta"
  writeXStringSet(seqs, file=in_file)

  # Build MEME command
  cmd <- sprintf("python %s/fasta-dinucleotide-shuffle.py -f %s -t _shuff -s 1 -c %s > %s", meme_path, in_file, nseqs, out_file)
  message("Running MEME fasta-dinucleotide-shuffle with the following command:")
  message(cmd)

  # Run meme
  system(cmd)
  message("Done.")

  # Collect results
  shuffled_seqs <- readDNAStringSet(out_file)

  # Cleanup temporary files
  system(sprintf("rm %s %s", in_file, out_file))

  return(shuffled_seqs)
}


getAllAltSeqs <- function(gr){
  # Get all SNP alt seqs 
  ##########################
  # gr: the SNP genomic range

  alt_seqs <- sapply(seq_along(gr), function(x){
  snp <- gr[x]
  out <- tryCatch(
    createRefAltseq(
    refalt=snp$linked_refalt,
    sequence=snp$refseq[[1]],
    snp_pos=snp$snp_start,
    seq_start=start(snp),
    rsid=snp$linked_SNP
    ),
    error=function(cond){
      message(sprintf("Choked on %s", snp$linked_SNP))
      message("Here's the original error message:")
      message(cond)
      return(snp$refseq[[1]])
    }
  )
  out
  }) %>% DNAStringSet()
  gr$altseq <- alt_seqs
  gr
}


# Create 'ref/alt' versions of each peak sequence
createRefAltseq <- function(refalt, sequence, snp_pos, seq_start, rsid=NULL){
  # Return the 'alt' sequence of a peak based on the ref and alt alleles of a 
  # provided SNP
  ############################
  # refalt: The ref and alt alleles as a character in the form of 'ret:alt'
  # sequence: the DNAString sequence of the sequence containing the snp
  # snp_pos: the genomic position of the SNP
  # seq_start: the genomic position of the start of the sequence

  # First split ref and alt
  splt <- strsplit(refalt, split=":")[[1]]
  ref <- splt[1]
  alt <- splt[2]

  # Make sure ref actually corresponds to the sequence seen at the expected position
  idx1 <- snp_pos - seq_start + 1
  idx2 <- idx1 + nchar(ref) - 1

  # If the ref is not fully contained within the peak, we will just discard that SNP
  if(any(c(idx1, idx2) > length(sequence))){
    message(sprintf("Warning, rsid %s ref allele was not contained within provided seq!", rsid))
    return(NA)
  }

  # If ref doesn't match native, we have a problem
  native <- sequence[idx1:idx2]
  if(ref != as.character(native)){
    message(sprintf("Warning, rsid %s did not have a matching ref and native!", rsid))
    message(sprintf("ref:alt = %s; native = %s", refalt, native))
    return(NA)
  }

  # Now create the new sequence
  newseq <- replaceAt(sequence, IRanges(idx1:idx2), alt)
  newseq
}

##################################################################
# Create all shuffled sequences
##################################################################

# We also want to create 'null' sequences that will have shuffled dinucleotide sequences
# per sequence but with the same ref allele at the same position as the 
# original sequences

# Prepare list of genomic ranges to run
gr_names <- names(gr_list)

for(grn in gr_names){
  gr <- gr_list[[grn]]
  sgr <- shuffleSNPseqs(gr, nseqs=3)
  sgrn <- paste0(grn, "_shuff")
  gr_list[[sgrn]] <- sgr
}

##################################################################
# Create all alternate alleles for original sequence contexts
##################################################################

fasta_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/snp_fastas"
dir.create(fasta_dir, showWarnings = FALSE, recursive = TRUE)

# Get new gr_names with shuffled sequences
gr_names <- names(gr_list)

for(grn in gr_names){
  message(sprintf("Getting sequences for %s...", grn))

  original_gr <- gr_list[[grn]]

  # Get all alt sequences in parallel
  ##################################################################
  by_chr <- split(original_gr, seqnames(original_gr))[1:23]
  chr_names <- names(by_chr)
  final_gr <- mclapply(
    names(by_chr), 
    function(x){
        gr <- by_chr[[x]]
        gr <- getAllAltSeqs(gr)
        message(sprintf("Finished %s...", x))
        gr
      }, 
    mc.cores=ncores) %>% as(., "GRangesList")
  names(final_gr) <- chr_names
  ###################################################################

  # Save sequences for use in model fitting
  # (Save split by chromosome for easier parallelization with gkmexplain)
  for(chr in names(by_chr)){
    gr <- final_gr[[chr]]
    ref_seqs <- gr$refseq
    alt_seqs <- gr$altseq

    names(ref_seqs) <- paste(gr$peak, gr$linked_SNP, gr$snp_rel_pos, sep="_")
    names(alt_seqs) <- paste(gr$peak, gr$linked_SNP, gr$snp_rel_pos, sep="_")

    writeXStringSet(ref_seqs, file=paste0(fasta_dir, sprintf("/ref_snp_seqs.%s.%s.fasta", grn, chr)))
    writeXStringSet(alt_seqs, file=paste0(fasta_dir, sprintf("/alt_snp_seqs.%s.%s.fasta", grn, chr)))
  }

  final_gr <- unlist(unname(final_gr))
  saveRDS(final_gr, file = paste0(fasta_dir, sprintf("/%s.rds", grn)))
}

##################################################################

