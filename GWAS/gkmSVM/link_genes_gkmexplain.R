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
  library(SuscrofaTxdb.11.108.july)
  library(parallel)
})

# Set Threads to be used
ncores <- 5
addArchRThreads(threads = ncores)

# Get additional functions, etc.:
scriptPath <- "/home/adufour/work/scScalpChromatin"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
#source(paste0(scriptPath, "/GO_wrappers.R"))

# set working directory (The directory of the full preprocessed archr project)
gkm_res_dir <- "/home/adufour/work/gskm/snp_results"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

##########################################################################################
# Preparing Data
##########################################################################################

load(file = "/home/adufour/work/rds_storage/omics/archr_all_v6.RData")
plotDir <- "/home/adufour/work/notebook/plots/gskm"

# Get all peaks
allPeaksGR <- getPeakSet(archrproj_sub)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

# Correlation cutoff for identifying p2g linkages
corrCutoff <- 0.5

# Retrieve GImat prior to re-assigning p2g links
GImat <- getMatrixFromProject(archrproj_sub, useMatrix="GeneExpressionMatrix")

##########################################################################################
# Links fine-mapped SNPs to candidate genes using Peak-to-Gene links
##########################################################################################

fmGR <- readRDS("/home/adufour/work/gskm/snp_fasta/250bpSNPCentered.rds")

# Load full project p2g links, plot loops, etc.
full_p2gGR <- readRDS(file="/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/scATAC_preprocessing/fine_clustered/multilevel_p2gGR.rds") # NOT merged or correlation filtered
plot_loop_list <- readRDS(file="/oak/stanford/groups/wjg/boberrey/hairATAC/results/scATAC_preprocessing/fine_clustered/multilevel_plot_loops.rds")

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
p2gGR <- getP2G_GR(archrproj_sub, corrCutoff=corrCutoff)

pol <- findOverlaps(p2gGR, fmGR, type="any", maxgap=-1L, ignore.strand=TRUE)
expandFMGR <- fmGR[to(pol)]
expandFMGR$linkedGene <- p2gGR[from(pol)]$symbol
expandFMGR$SNP_to_gene <- paste(expandFMGR$peakName, expandFMGR$linkedGene, sep="_")

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
snp_to_gene_df <- mcols(expandFMGR) %>% as.data.frame() %>% group_by(peakName) %>% summarize(genes=paste(unique(linkedGene), collapse=";")) %>% as.data.frame()
gene_vec <- snp_to_gene_df$genes
names(gene_vec) <- snp_to_gene_df$peakName

full_sig_res$linkedGenes <- gene_vec[full_sig_res$snp]

mcols(expandFMGR) %>% as.data.frame() %>% group_by(peakName) %>% summarize(genes=paste(unique(linkedGene))) %>% as.data.frame()

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
  region <- snp_info$region[1]

  # Get SNP region
  snp_gr <- gr[gr$peakName == snp] %>% resize(width=50, fix="center")
  true_region <- (snp_gr %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})

  # Get info on suspected motif match, linked genes, etc
  #linked_genes <- snp_info$linkedGenes[1]
  #top_motifs <- snp_info$top_motifs[1]
  #fm_prob <- snp_info$fm_probs[1]

  # Get matrices for plotting
  ref_matrix <- gkm_explain_output$gkmexplain_output[[ref_group]]$seq_matrices[[region]][,indices]

  # Generate plot for all matrices
  upper_lim <- max(ref_matrix) * 1.2
  lower_lim <- min(ref_matrix) * 1.2
  mat_list <- list("ref_allele"=ref_matrix)

  p <- (
    ggseqlogo(mat_list, method='custom', seq_type='dna', ncol=1) 
    + ylab('gkmexplain importance')
    + ylim(lower_lim, upper_lim)
    + ggtitle(paste0(celltype, " - ", snp, " | ", true_region
      ))
    )
  p
}

snp_gkmexplain$gkmexplain_output[["C5-ref_snp_seqs"]]$seq_matrices[["3_126328970_126329220"]][,101:151]

gene <- "GRHL1"

gene_peak_position <- p2gGR[p2gGR$symbol == gene]

snp_list <- paste0(seqnames(gene_peak_position), "_", start(gene_peak_position)+125, "_", end(gene_peak_position)-125)

pdf(paste0(gkm_res_dir, "/most_prominent_hc_logos.pdf"), width=12, height=7)
plotList <- list()
snp <- snp_list[1]
ct <- "C5"
plotList <- plotSNPimportance(snp, significant_hits, snp_gkmexplain, celltype=ct, gr=fmGR)
plotList
dev.off()


# Most prominent, high-effect SNP hits
pdf(paste0(gkm_res_dir, "/most_prominent_hc_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in snp_list){
  snp <- i
  ct <- "C2"
  plotList[[i]] <- plotSNPimportance(snp, significant_hits, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()

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
disPicsGR <- expandFMGR
promoterGR <- promoters(getGenes(archrproj_sub))
candidate_GR <- disPicsGR
index_df <- data.frame(index=candidate_GR$index_SNP, candidate=candidate_GR$linked_SNP)
index_df <- index_df[!duplicated(index_df$candidate),]
index_map <- index_df$index
names(index_map) <- index_df$candidate

# The original index SNP is often not in a peak, so in order to find it we need to load the unfiltered fm gr
#index_GR <- dis_full_fm_gr[dis_full_fm_gr$linked_SNP %in% index_map] %>% unique()

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

# NamedClust
p <- plotBrowserTrack(
    ArchRProj = archrproj_sub, 
    groupBy = "Clusters",
    geneSymbol = gene,
    pal = c("#B6E0EE", "#00b894", "#78c4ce", "#fdcb6e", "#996ea5", "#ce8787"),
    #plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    sizes = c(8, 0.8, 1.25, 0.5),
    #region = plotRegions, 
    loops = getPeak2GeneLinks(archrproj_sub, corCutOff=corrCutoff), # All peak-to-gene loops
    tileSize=250,
    minCells=100,
    #title = markerGenes
)

archrproj_sub@projectMetadata$outputDirectory <- "/home/adufour/work/notebook/plots/gskm"

plotPDF(plotList = p, 
    name = paste0(gene, "_tracks.pdf"), 
    ArchRProj = archrproj_sub, 
    addDOC = FALSE, 
    width = 8, height = 7)


##################################################
# Violin plots of RNA expression for select genes
##################################################

markerGenes <- c("GRHL1", "GRHL2")

data_mat <- assays(GImat)[[1]]
rownames(data_mat) <- rowData(GImat)$name
sub_mat <- data_mat[markerGenes,]

grouping_data <- data.frame(cluster=factor(archrproj_sub$Clusters, 
  ordered=TRUE))
rownames(grouping_data) <- getCellNames(archrproj_sub)
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

pList <- list()
for(gn in markerGenes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  set.seed(1)
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df %>% as.data.frame()

  covarLabel <- "cluster"  

  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=gene, fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
    + scale_color_manual(values=c("#B6E0EE", "#00b894", "#78c4ce", "#fdcb6e", "#996ea5", "#ce8787"), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=c("#B6E0EE", "#00b894", "#78c4ce", "#fdcb6e", "#996ea5", "#ce8787"))
    + guides(fill=guide_legend(title=covarLabel), 
      colour=guide_legend(override.aes = list(size=5)))
    + ggtitle(gn)
    + xlab("")
    + ylab("RNA Expression")
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
