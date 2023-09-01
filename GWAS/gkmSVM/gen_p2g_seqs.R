#!/usr/bin/env Rscript

##########################################################################
# Create fasta files containing SNP-modified peak sequences
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
  #library(BSgenome.Hsapiens.UCSC.hg38.masked)
  library(SuscrofaTxdb.11.108.july)
  library(parallel)
  library(ComplexHeatmap)
})

# Set Threads to be used
ncores <- 5
addArchRThreads(threads = ncores)

# Get additional functions, etc.:
scriptPath <- "/home/adufour/work/scScalpChromatin"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/fine_clustered"
fm_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/gwas/PICS2"

#Set/Create Working Directory to Folder
setwd(wd)

load(file = "/home/adufour/work/rds_storage/omics/archr_all_v6.RData")
corrCutoff <- 0.5
p2gGR <- getP2G_GR(archrproj_sub, corrCutoff=corrCutoff)

##########################################################################################
# Preparing Data
##########################################################################################

# 250bp window centered on SNP position
snps_250bp_center_gr <- resize(p2gGR, width=251, fix="center")
names(snps_250bp_center_gr) <- (snps_250bp_center_gr %>% {paste0(seqnames(.), "_", start(.), "_", end(.), "_", .$linked_SNP)})
snps_250bp_center_gr$refseq <- BSgenome::getSeq(SuscrofaTxdb.11.108.july, names=snps_250bp_center_gr)
snps_250bp_center_gr$snp_rel_pos <- start(snps_250bp_center_gr)

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
  refs <- gr$refseq

  # Get shuffled sequences
  shuffled_seqs <- getDinucShuffledSeqs(refs, nseqs=nseqs)
  return(shuffled_seqs)
}


getDinucShuffledSeqs <- function(seqs, nseqs=1, meme_path="/home/adufour/work/programms/meme-5.5.3/scripts"){
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
  cmd <- sprintf("python %s/fasta-dinucleotide-shuffle-py3.in -f %s -t _shuff -s 1 -c %s > %s", meme_path, in_file, nseqs, out_file)
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

fasta_dir <- "/home/adufour/work/gskm/snp_fasta"
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
        #gr <- getAllAltSeqs(gr)
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
    #alt_seqs <- gr$altseq

    names(ref_seqs) <- paste(seqnames(gr), start(gr), end(gr), sep="_")
    #names(alt_seqs) <- paste(gr$peak, gr$linked_SNP, gr$snp_rel_pos, sep="_")

    writeXStringSet(ref_seqs, file=paste0(fasta_dir, sprintf("/ref_snp_seqs.%s.%s.fasta", grn, chr)))
    #writeXStringSet(alt_seqs, file=paste0(fasta_dir, sprintf("/alt_snp_seqs.%s.%s.fasta", grn, chr)))
  }

  final_gr <- unlist(unname(final_gr))
  saveRDS(final_gr, file = paste0(fasta_dir, sprintf("/%s.rds", grn)))
}

##################################################################

