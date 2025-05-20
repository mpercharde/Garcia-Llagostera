# plotting motif distribution using FIMO

#################################################################################
# load packages
#################################################################################

### load packages
packages = c("tidyverse", "magrittr", "glue", "here", "data.table", "furrr", 
             "BiocParallel", "BiocManager", "knitr", "readr", "kableExtra", 
             "DT", "pheatmap", "RColorBrewer", "GenomicFeatures", "DESeq2", 
             "PCAtools", "IHW", "ashr", "variancePartition", "AnnotationDbi", 
             "org.Mm.eg.db", "viridis", "RUVSeq", "ComplexHeatmap", "Cairo", 
             "clusterProfiler", "eulerr", "openxlsx", "biomaRt", "BSgenome", 
             "BSgenome.Mmusculus.UCSC.mm39", "memes", "universalmotif", 
             "MotifDb", "ggtree", "motifStack", "cowplot", "ggsignif", 
             "rstatix", "ggpubr", "ROCR")

for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only=TRUE))
}

rm(p, packages)

#################################################################################
# set functions and more
#################################################################################

## set paths
motif_dists_outputs <- list(
  fimo_res = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss3kb/results/fimo", # for fimo results
  downst_res = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss3kb/results/downstream", # for downstream results
  fig_motif_dists_tss3kb = "../analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss3kb/figures", # for ggsaver() function
  fig_motif_dists_tss3kb.do = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss3kb/figures", # for dev.off() operation
  fig_motif_dists_tss3kb.sig_test.do = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss3kb/figures/significance_testing", # for significance testing
  res_motif_dists_tss3kb = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss3kb/results/downstream"
)

"%ni%" = Negate("%in%")

pts_to_cm = function(pts) {
  return(pts * 0.0352778)
}

cm_to_inch = function(cm) {
  return(cm / 2.54)
}

# load custom functions
ggsaver = function(
    basename, plot=ggplot2::last_plot(), device=c("png", "svg"), ...) {
  # wrapper around ggplot2::ggsave automating saving of multiple image formats
  args = c(as.list(environment()), list(...))
  device = c(device)
  args$basename = NULL
  written = c()
  for (d in device) {
    args$filename = glue("{basename}.{d}")
    args$device = d
    do.call(ggplot2::ggsave, args)
    written = c(written, args$filename)
  }
  return(written)
}


#################################################################################
# load in and process data
#################################################################################

## initialize empty list
fimo_pre <- list()
## read in each tsv and add to list
c_pre <- read.table(glue("{motif_dists_outputs$fimo_res}/c_genes/fimo-c_genes.pre_proc.anno.bed"), 
                    header = FALSE, sep = "\t", stringsAsFactors = FALSE)
imi_pre <- read.table(glue("{motif_dists_outputs$fimo_res}/imi_genes/fimo-imi_genes.pre_proc.anno.bed"), 
                      header = FALSE, sep = "\t", stringsAsFactors = FALSE)
ir_pre <- read.table(glue("{motif_dists_outputs$fimo_res}/ir_genes/fimo-ir_genes.pre_proc.anno.bed"), 
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)
mi_pre <- read.table(glue("{motif_dists_outputs$fimo_res}/mi_genes/fimo-mi_genes.pre_proc.anno.bed"), 
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)
### add to list
fimo_pre$c <- c_pre
fimo_pre$imi <- imi_pre
fimo_pre$ir <- ir_pre
fimo_pre$mi <- mi_pre
## process each tsv
### make a vector of column names
col_names <- c("chr", "start_motif", "end_motif", "strand_motif", "motif_score", "motif_pvalue", "motif_qvalue", 
               "motif_id", "matched_sequence", "chr_input", "start_input", "end_input", "gene_id", "gene_name", 
               "strand_gene", "bp_overlap")
### initialize empty list to store processed data
fimo_proc <- list()
### process fimo results
for (i in seq_along(fimo_pre)) {
  a <- fimo_pre[[i]]
  
  n <- names(fimo_pre)[i]
  print(n)
  
  colnames(a) <- col_names
  
  a <- a %>% select(-chr_input)
  
  a$tf_name <- sub("\\..*", "", a$motif_id)
  
  a$gene_set <- n
  
  a$relative_orientation <- ifelse(a$strand_gene == a$strand_motif, "same", "reverse")
  
  a <- a[, c("chr", "start_motif", "end_motif", "strand_motif", "start_input", "end_input", "strand_gene", "motif_score", 
             "motif_pvalue", "motif_qvalue", "matched_sequence", "bp_overlap", "motif_id", "tf_name", "gene_id", "gene_name", 
             "gene_set", "relative_orientation")]
  
  fimo_proc[[n]] <- a
  
  ## report completion
  print("finished")
}
### check output
names(fimo_proc)
colnames(fimo_proc[[3]])
head(fimo_proc[[1]])
tail(fimo_proc[[2]])

## export processed dfs as csvs
  ### check names
  names(fimo_proc)
### export
write.csv(fimo_proc[["c"]], glue("{motif_dists_outputs$fimo_res}/c_genes/fimo-c_genes.post_proc.anno.csv"), 
          row.names = FALSE)
write.csv(fimo_proc[["imi"]], glue("{motif_dists_outputs$fimo_res}/imi_genes/fimo-imi_genes.post_proc.anno.csv"), 
          row.names = FALSE)
write.csv(fimo_proc[["ir"]], glue("{motif_dists_outputs$fimo_res}/ir_genes/fimo-ir_genes.post_proc.anno.csv"), 
          row.names = FALSE)
write.csv(fimo_proc[["mi"]], glue("{motif_dists_outputs$fimo_res}/mi_genes/fimo-mi_genes.post_proc.anno.csv"), 
          row.names = FALSE)

f_dists <- list()
for (i in seq_along(fimo_proc)) {
  a <- fimo_proc[[i]]
  n <- names(fimo_proc)[i]
  print(n)
  
  a <- a %>% 
    dplyr::filter(motif_pvalue <= 0.00001)
  
  f <- a %>% 
    dplyr::filter(tf_name %in%
                    c("KLF4", "NFKB1", "IRF3", "IRF9", "STAT1", "STAT2"))
  
  f_dists[[n]] <- f
  ## report completion
  print("finished")
}


#################################################################################
# plot motif distributions over sequences
#################################################################################


## calculate motif locations to map distribution across input sequence and normalize for strandedness
dist_calc_all <- list() ## irrespective of orientation
dist_calc_same <- list() ## same sense orientation only for motif with gene
### calculate the number of occurrences of each motif over each gene in the gene sets
for (i in seq_along(f_dists)) {
  a <- f_dists[[i]]
  ## extract and print df name
  n <- names(f_dists)[i]
  print(n)
  
  ## calculate center motif position in new column, 'center_motif'
  a <- a %>% 
    mutate(center_motif = (start_motif + end_motif) / 2)
  ### specify as integer to round center_motif down to the nearest whole number
  a$center_motif <- as.integer(a$center_motif)
  
  ## calculate position of each motif relative to TSS oriented upstream and downstream for all genes
  a <- a %>%
    mutate(orient_position = case_when(
      strand_gene == "+" ~ end_input - center_motif,        ## for + strand genes
      strand_gene == "-" ~ (6000 - (end_input - center_motif))        ## for - strand genes
    ))
  
  ## add all results to relevant list
  dist_calc_all[[n]] <- a
  
  ## filter for same sense orientation and add to list
  dist_calc_same[[n]] <- a %>%
    filter(relative_orientation == "same")
  
  ## report completion
  print("finished")
}

## setting colors
dist_col <- brewer.pal(6, "Dark2") 

## plot motif distributions across +/-3kb region around TSS - same
for (i in seq_along(dist_calc_same)) {
  a <- dist_calc_same[[i]]
  ## extract and print df name
  n <- names(dist_calc_same)[i]
  print(n)
  
  ### plot stacked histograms of raw motif counts for each gene set
  #### fixed scale
  g1 <- ggplot(a, aes(x = orient_position)) +
    geom_histogram(bins = 100, alpha = 0.6, fill = "steelblue", color = NA) + 
    facet_wrap(~tf_name, ncol = 1, scales = "fixed") +  
    theme_minimal() +
    labs(title = glue("Motif distribution - {n} genes (fixed scale)"),
         caption = "3kbTSS, p-val <= 1e-5",
         x = "Position Along Sequence (bp)",
         y = "Count of Occurrences") +
    scale_x_continuous(limits = c(0, 6000)) +  
    #scale_fill_brewer(palette = "Dark2") + 
    theme(strip.text = element_text(face = "bold", size = 12),
          text = element_text(size = 14))
  #### independent scales
  g2 <- ggplot(a, aes(x = orient_position)) +
    geom_histogram(bins = 100, alpha = 0.6, fill = "steelblue", color = NA) + 
    facet_wrap(~tf_name, ncol = 1, scales = "free_y") +  
    theme_minimal() +
    labs(title = glue("Motif distribution - {n} genes (unfixed scale)"),
         caption = "3kbTSS, p-val <= 1e-5",
         x = "Position Along Sequence (bp)",
         y = "Count of Occurrences") +
    scale_x_continuous(limits = c(0, 6000)) +  
    #scale_fill_brewer(palette = "Dark2") + 
    theme(strip.text = element_text(face = "bold", size = 12),
          text = element_text(size = 14))
  ### plot density curves of each TF
  #### using curve adjust 1
  g3.1 <- ggplot(a, aes(x = orient_position, color = tf_name)) +
    geom_density(adjust = 1, size = 1.2) +  
    theme_minimal() +
    labs(title = glue("Motif distribution - {n} genes"),
         caption = "3kbTSS, p-val <= 1e-5",
         x = "Position Along Sequence (bp)",
         y = "Density of Occurrences",
         color = "Motif") +
    scale_x_continuous(limits = c(0, 6000)) + 
    #scale_fill_brewer(palette = "Dark2") + 
    scale_fill_manual(values = dist_col) +
    theme(text = element_text(size = 14))
  #### using curve adjust 2
  g3.2 <- ggplot(a, aes(x = orient_position, color = tf_name)) +
    geom_density(adjust = 2, size = 1.2) +  
    theme_minimal() +
    labs(title = glue("Motif distribution - {n} genes"),
         caption = "3kbTSS, p-val <= 1e-5",
         x = "Position Along Sequence (bp)",
         y = "Density of Occurrences",
         color = "Motif") +
    scale_x_continuous(limits = c(0, 6000)) + 
    #scale_fill_brewer(palette = "Dark2") + 
    scale_fill_manual(values = dist_col) +
    theme(text = element_text(size = 14))
  ### plot histogram with density curves
  g4 <- ggplot(a, aes(x = orient_position, color = tf_name)) +
    geom_histogram(aes(y = ..density.., fill = tf_name), bins = 50, alpha = 0.4, position = "identity", color = NA) +  
    geom_density(size = 1.2, adjust = 1, fill = NA) +  
    theme_minimal() +
    labs(title = glue("Motif distribution - {n} genes"),
         caption = "3kbTSS, p-val <= 1e-5",
         x = "Position Along Sequence (bp)",
         y = "Density of Occurrences",
         fill = "Motif",
         color = "Motif") +
    scale_x_continuous(limits = c(0, 6000)) +  
    #scale_fill_brewer(palette = "Dark2") + 
    scale_fill_manual(values = dist_col) +
    scale_color_manual(values = dist_col) + 
    theme(text = element_text(size = 14))
  
  ### export plots
  ggsaver(
    glue("fin_v1_filtered-motif_distribution_over_TSSregion-histogram_stacked-raw_motif_counts-fixed_scale-same_sense.p0.00001.{n}_genes.tss3kb.fimo"), 
    g1, path=motif_dists_outputs$fig_motif_dists_tss3kb.do, 
    height=cm_to_inch(26), width=cm_to_inch(25), units="in",
    device=c("png", "pdf"))
  ggsaver(
    glue("fin_v1_filtered-motif_distribution_over_TSSregion-histogram_stacked-raw_motif_counts-unfixed_scale-same_sense.p0.00001.{n}_genes.tss3kb.fimo"), 
    g2, path=motif_dists_outputs$fig_motif_dists_tss3kb.do, 
    height=cm_to_inch(26), width=cm_to_inch(25), units="in",
    device=c("png", "pdf"))
  ggsaver(
    glue("fin_v1_filtered-motif_distribution_over_TSSregion-density_curves_adj_1-same_sense.p0.00001.{n}_genes.tss3kb.fimo"), 
    g3.1, path=motif_dists_outputs$fig_motif_dists_tss3kb.do, 
    height=cm_to_inch(10), width=cm_to_inch(15), units="in",
    device=c("png", "pdf"))
  ggsaver(
    glue("fin_v1_filtered-motif_distribution_over_TSSregion-density_curves_adj_2-same_sense.p0.00001.{n}_genes.tss3kb.fimo"), 
    g3.2, path=motif_dists_outputs$fig_motif_dists_tss3kb.do, 
    height=cm_to_inch(10), width=cm_to_inch(15), units="in",
    device=c("png", "pdf"))
  ggsaver(
    glue("fin_v1_filtered-motif_distribution_over_TSSregion-histogram_density_curves-same_sense.p0.00001.{n}_genes.tss3kb.fimo"), 
    g4, path=motif_dists_outputs$fig_motif_dists_tss3kb.do, 
    height=cm_to_inch(10), width=cm_to_inch(15), units="in",
    device=c("png", "pdf"))
  
  ## report completion
  print("finished")
}

