# motif analysis using FIMO

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
  fimo_res = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss1kb/results/fimo", # for fimo results
  downst_res = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss1kb/results/downstream", # for downstream results
  fig_motif_dists_tss1kb = "../analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss1kb/figures", # for ggsaver() function
  fig_motif_dists_tss1kb.do = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss1kb/figures", # for dev.off() operation
  fig_motif_dists_tss1kb.sig_test.do = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss1kb/figures/significance_testing", # for significance testing
  res_motif_dists_tss1kb = "./analysis/3.rna_downstream/mm39/motif_analysis/venn_motif_distances/tss1kb/results/downstream"
)

"%ni%" = Negate("%in%")

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
## check output and structure
head(c_pre)
dim(c_pre)
### add to list
fimo_pre$c <- c_pre
fimo_pre$imi <- imi_pre
fimo_pre$ir <- ir_pre
fimo_pre$mi <- mi_pre
## process each tsv (extract motif names, reorder columns, add column names, extract gene information (name and strandedness), add further gene annotation, load in bed files and add information on TSS1kb start and end, add a column specifying same sense orientation (as the gene) or not, add a column for "analysis" or "gene_set" (to record c_genes, imi_genes, calculate motif coordinates, etc.), in case I ever combine all of these)
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

## export processed dfs as csvs
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


#################################################################################
# measure and plot motif occurrences in gene sets
#################################################################################

# calculate the number of each motif per gene and the percentage of genes with each motif, for all and separately for same-stranded
### initialize empty lists to store results of the calculations
motif_counts <- list()
motif_perc_ss <- list()
motif_perc_a <- list()
### calculate the number of occurrences of each motif over each gene in the gene sets
for (i in seq_along(fimo_proc)) {
  a <- fimo_proc[[i]]
  n <- names(fimo_proc)[i]
  print(n)
  ## determine the number of unique genes in the set
  x <- length(unique(a$gene_id))  
    print(glue("There are {x} unique genes in gene set {n}"))
  a <- a %>% 
    dplyr::filter(motif_qvalue <= 0.25)
  ## calculate the number of occurrences of each motif per gene SEPARATELY for each strand orientation
    ## (output will have the specified columns in group_by, and one for count)
  g_s <- a %>%
    group_by(gene_id, gene_name, tf_name, gene_set, relative_orientation) %>%
    summarise(count = n(), .groups = "drop")
  ## calculate the TOTAL number of occurrences of each motif per gene, independently of strand orientation (aka total)
  g_a <- a %>%
    group_by(gene_id, gene_name, tf_name, gene_set) %>%
    summarise(count = n(), .groups = "drop")
  ### add a column to specify number irrespective of strand orientation
  g_a$relative_orientation <- "total"
  ## combine the results of both calculations
  g_c <- bind_rows(g_s, g_a)
  ## sort alphabetically by the value of gene_id and then tf_name
  g_sort <- g_c %>% arrange(gene_id, tf_name)
  ## add results df to motif_counts list
  motif_counts[[n]] <- g_sort
  
  ## calculate the percentage of genes in each set containing motifs for each TF - SAME STRANDED ONLY
  ### filter for same strand orientation
  ss <- g_sort %>%
    filter(relative_orientation == "same")
  ### calculate the percentage of genes containing each motif
  mp_s <- ss %>%
    dplyr::select(gene_id, tf_name) %>%  # Select relevant columns
    distinct() %>%  # Keep only unique gene-motif pairs
    group_by(tf_name) %>%  # Group by motif type (tf_name)
    summarise(genes_with_motif = n_distinct(gene_id)) %>%
    mutate(percentage = (genes_with_motif / x) * 100)
  ### add gene_set information to the tibble
  mp_s <- mp_s %>%
    mutate(gene_set = n)
  ### save the results tibble to list
  motif_perc_ss[[n]] <- mp_s
  ## calculate the percentage of genes in each set containing motifs for each TF - TOTAL ONLY
  ### filter for same strand orientation
  t <- g_sort %>%
    filter(relative_orientation == "total")
  ### calculate the percentage of genes containing each motif
  mp_a <- t %>%
    dplyr::select(gene_id, tf_name) %>%  # Select relevant columns
    distinct() %>%  # Keep only unique gene-motif pairs
    group_by(tf_name) %>%  # Group by motif type (tf_name)
    summarise(genes_with_motif = n_distinct(gene_id)) %>% 
    mutate(percentage = (genes_with_motif / x) * 100)
  ### add gene_set information to the tibble
  mp_a <- mp_a %>%
    mutate(gene_set = n)
  ### save the results tibble to list
  motif_perc_a[[n]] <- mp_a
  
  ## report completion
  print("finished")
}
### check list output
names(motif_counts)
head(motif_counts[[1]])
tail(motif_counts[[3]])
names(motif_perc_a)
names(motif_perc_ss)
head(motif_perc_ss[[1]])
tail(motif_perc_ss[[2]])
## collect the results of motif counts per gene into one df for export
counts_df <- bind_rows(motif_counts[[1]], motif_counts[[2]], motif_counts[[3]], motif_counts[[4]])
### export as csv
write.csv(counts_df, glue("{motif_dists_outputs$downst_res}/fimo-motif_counts_per_gene.main_venn.q0.25.tss1kb.csv"), 
          row.names = FALSE)
## collect the results of motif percentage calculations into one df for export (and plotting)
### collect results
perc_df_ss <- bind_rows(motif_perc_ss[[1]], motif_perc_ss[[2]], motif_perc_ss[[3]], motif_perc_ss[[4]])
perc_df_a <- bind_rows(motif_perc_a[[1]], motif_perc_a[[2]], motif_perc_a[[3]], motif_perc_a[[4]])
### export as csv
write.csv(perc_df_ss, glue("{motif_dists_outputs$downst_res}/fimo-percent_genes_w_motifs.main_venn.same_strand.q0.25.tss1kb.csv"), 
          row.names = FALSE)
write.csv(perc_df_a, glue("{motif_dists_outputs$downst_res}/fimo-percent_genes_w_motifs.main_venn.total.q0.25.tss1kb.csv"), 
          row.names = FALSE)

# plot the results of the calculations of motif counts per gene (bubble plot and grouped bar chart) - qval 0.25
## plot a bubble plot of the percentage of genes containing motifs - SAME SENSE ONLY
### check details
unique(perc_df_ss$gene_set)
unique(perc_df_ss$tf_name)
### plot
gp_ss <- ggplot(perc_df_ss, aes(x = gene_set, y = tf_name, colour=gene_set)) +
  theme_classic() + 
  #xlab('') +
  #ylab('') +
  labs(title = "Percent of genes with motif",
       x = "",
       y = "",
       caption = "TSS1kb, qval 0.25, same strand") +
  theme(plot.title = element_text(size = 11)) +
  geom_point(aes(size=percentage), stroke=0) +
  scale_size_continuous(range = c(1.4,8)) +  ## adjust the range of point sizes, ensuring each is visually noticeable
  ylim(sort(unique(perc_df_ss$tf_name), decreasing=TRUE)) + ## sort motif order
  guides(colour='none')
### view output
gp_ss
### export plot
ggsaver(
  "bubble_plot-percent_genes_with_motif.main_venn_sets.same_orientation.tss1kb.q0.25.fimo", 
  gp_ss, path=motif_dists_outputs$fig_motif_dists_tss1kb.do, 
  height=cm_to_inch(12), width=cm_to_inch(8), units="in",
  device=c("png", "pdf"))
#### cleanup
rm(gp_ss)
## plot a bubble plot of the percentage of genes containing motifs - TOTAL
### check details
unique(perc_df_a$gene_set)
unique(perc_df_a$tf_name)
### plot
gp_a <- ggplot(perc_df_a, aes(x = gene_set, y = tf_name, colour=gene_set)) +
  theme_classic() + 
  #xlab('') +
  #ylab('') +
  labs(title = "Percent of genes with motif",
       x = "",
       y = "",
       caption = "TSS1kb, qval 0.25, total hits") +
  theme(plot.title = element_text(size = 11)) +
  geom_point(aes(size=percentage), stroke=0) +
  scale_size_continuous(range = c(1.4,8)) +  ## adjust the range of point sizes, ensuring each is visually noticeable
  ylim(sort(unique(perc_df_ss$tf_name), decreasing=TRUE)) + ## sort motif order
  guides(colour='none')
### view output
gp_a
### export plot
ggsaver(
  "bubble_plot-percent_genes_with_motif.main_venn_sets.total.tss1kb.q0.25.fimo", 
  gp_a, path=motif_dists_outputs$fig_motif_dists_tss1kb.do, 
  height=cm_to_inch(12), width=cm_to_inch(8), units="in",
  device=c("png", "pdf"))
#### cleanup
rm(gp_a)
## plotting bubble plot and grouped bar chart (in order to compare to results of previous analysis) of the average number of motifs per gene - SAME SENSE ONLY
### prepare gene_set and relative_orientation as factors
counts_pl <- counts_df %>%
  mutate(relative_orientation = factor(relative_orientation, levels = c("reverse", "same", "total")),
         gene_set = factor(gene_set))  ## ensure gene_set is treated as a categorical variable
### compute the average motif occurrences per gene in each gene set -> same sense orientation ONLY
#### filter for 'same' category
counts_s <- counts_pl %>%
  filter(relative_orientation == "same")
#### count the number of unique genes in each gene set
gene_counts <- counts_s %>%
  group_by(gene_set) %>%
  summarise(unique_genes = n_distinct(gene_id), .groups = "drop")
#### calculate the sum of each motif's occurrences per gene set
motif_sums <- counts_s %>%
  group_by(gene_set, tf_name) %>%
  summarise(total_count = sum(count), .groups = "drop")
#### merge the gene count info and calculate average motif occurrences per gene
average_counts_s <- motif_sums %>%
  left_join(gene_counts, by = "gene_set") %>%
  mutate(avg_count_per_gene = total_count / unique_genes)
### plot the results in a bubble plot
### plotting bubble plot of the average number of motifs per gene - SAME SENSE ONLY
gp_ss <- ggplot(average_counts_s, aes(x = gene_set, y = tf_name, colour=gene_set)) +
  theme_classic() + 
  #xlab('') +
  #ylab('') +
  labs(title = "Average motif occurrences per gene",
       x = "",
       y = "",
       caption = "TSS1kb, qval 0.25, same-sense") +
  theme(plot.title = element_text(size = 10)) +
  geom_point(aes(size=avg_count_per_gene), stroke=0) +
  scale_size_continuous(range = c(1.4,8)) +  ## adjust the range of point sizes, ensuring each is visually noticeable
  ylim(sort(unique(average_counts_s$tf_name), decreasing=TRUE)) + ## sort motif order
  guides(colour='none')
#### view output
gp_ss
#### export the plot
ggsaver(
  "bubble_plot-avg_motif_counts_per_gene.main_venn_sets.same_orientation.tss1kb.q0.25.fimo", 
  gp_ss, path=motif_dists_outputs$fig_motif_dists_tss1kb.do, 
  height=cm_to_inch(14), width=cm_to_inch(10), units="in",
  device=c("png", "pdf"))
### plot the results in a grouped bar chart - SAME SENSE ONLY
g_bar <- ggplot(average_counts_s, aes(x = tf_name, y = avg_count_per_gene, fill = gene_set)) +
  geom_bar(stat = "identity", position = "dodge") +  # Grouped bars
  theme_bw() +
  labs(title = "Average motif occurrences per gene (TSS1kb, same-sense, qval 0.25)",
       x = "Motif",
       y = "Average Count per Gene",
       fill = "Gene Set") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5))
#### check the output
g_bar
#### export the plot
ggsaver(
  "barplot-avg_motif_counts_per_gene.main_venn_sets.same_orientation.tss1kb.q0.25.fimo", 
  g_bar, path=motif_dists_outputs$fig_motif_dists_tss1kb.do, 
  height=cm_to_inch(14), width=cm_to_inch(16), units="in",
  device=c("png", "pdf"))
#### cleanup
rm(g_bar, gp_ss)
## plotting bubble plot and grouped bar chart (in order to compare to results of previous analysis) of the average number of motifs per gene - TOTAL
### prepare gene_set and relative_orientation as factors
counts_pl <- counts_df %>%
  mutate(relative_orientation = factor(relative_orientation, levels = c("reverse", "same", "total")),
         gene_set = factor(gene_set))  ## ensure gene_set is treated as a categorical variable
### compute the average motif occurrences per gene in each gene set -> same sense orientation ONLY
#### filter for 'same' category
counts_s <- counts_pl %>%
  filter(relative_orientation == "total")
#### count the number of unique genes in each gene set
gene_counts <- counts_s %>%
  group_by(gene_set) %>%
  summarise(unique_genes = n_distinct(gene_id), .groups = "drop")
#### calculate the sum of each motif's occurrences per gene set
motif_sums <- counts_s %>%
  group_by(gene_set, tf_name) %>%
  summarise(total_count = sum(count), .groups = "drop")
#### merge the gene count info and calculate average motif occurrences per gene
average_counts_s <- motif_sums %>%
  left_join(gene_counts, by = "gene_set") %>%
  mutate(avg_count_per_gene = total_count / unique_genes)
### plot the results in a bubble plot
### plotting bubble plot of the average number of motifs per gene - SAME SENSE ONLY
gp_a <- ggplot(average_counts_s, aes(x = gene_set, y = tf_name, colour=gene_set)) +
  theme_classic() + 
  #xlab('') +
  #ylab('') +
  labs(title = "Average motif occurrences per gene",
       x = "",
       y = "",
       caption = "TSS1kb, qval 0.25, total hits") +
  theme(plot.title = element_text(size = 10)) +
  geom_point(aes(size=avg_count_per_gene), stroke=0) +
  scale_size_continuous(range = c(1.4,8)) +  ## adjust the range of point sizes, ensuring each is visually noticeable
  ylim(sort(unique(average_counts_s$tf_name), decreasing=TRUE)) + ## sort motif order
  guides(colour='none')
#### view output
gp_a
#### export the plot
ggsaver(
  "bubble_plot-avg_motif_counts_per_gene.main_venn_sets.total.tss1kb.q0.25.fimo", 
  gp_a, path=motif_dists_outputs$fig_motif_dists_tss1kb.do, 
  height=cm_to_inch(14), width=cm_to_inch(10), units="in",
  device=c("png", "pdf"))
### plot the results in a grouped bar chart - SAME SENSE ONLY
g_bar <- ggplot(average_counts_s, aes(x = tf_name, y = avg_count_per_gene, fill = gene_set)) +
  geom_bar(stat = "identity", position = "dodge") +  # Grouped bars
  theme_bw() +
  labs(title = "Average motif occurrences per gene (TSS1kb, total hits, qval 0.25)",
       x = "Motif",
       y = "Average Count per Gene",
       fill = "Gene Set") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5))
#### check the output
g_bar
#### export the plot
ggsaver(
  "barplot-avg_motif_counts_per_gene.main_venn_sets.total.tss1kb.q0.25.fimo", 
  g_bar, path=motif_dists_outputs$fig_motif_dists_tss1kb.do, 
  height=cm_to_inch(14), width=cm_to_inch(16), units="in",
  device=c("png", "pdf"))
#### cleanup
rm(g_bar, gp_a)

# version
## plot a bubble plot of the percentage of genes containing motifs - SAME SENSE ONLY
### check details
unique(perc_df_ss$gene_set)
unique(perc_df_ss$tf_name)
### plot
gp_ss <- ggplot(perc_df_ss, aes(x = gene_set, y = tf_name, colour=gene_set)) +
  theme_classic() + 
  #xlab('') +
  #ylab('') +
  labs(title = "Percent of genes with motif",
       x = "",
       y = "",
       caption = "TSS1kb, qval 0.25, same strand") +
  theme(plot.title = element_text(size = 11)) +
  geom_point(aes(size=percentage), stroke=0) +
  scale_size_continuous(range = c(1.4,8)) +  ## adjust the range of point sizes, ensuring each is visually noticeable
  ylim(sort(unique(perc_df_ss$tf_name), decreasing=TRUE)) + ## sort motif order
  guides(colour='none')
### view output
gp_ss
### export plot
ggsaver(
  "bubble_plot-percent_genes_with_motif.main_venn_sets.same_orientation.tss1kb.q0.25.wide_version.fimo", 
  gp_ss, path=motif_dists_outputs$fig_motif_dists_tss1kb.do, 
  height=cm_to_inch(11), width=cm_to_inch(12), units="in",
  device=c("png", "pdf"))
#### cleanup
rm(gp_ss)
## plot a bubble plot of the percentage of genes containing motifs - TOTAL
### check details
unique(perc_df_a$gene_set)
unique(perc_df_a$tf_name)
### plot
gp_a <- ggplot(perc_df_a, aes(x = gene_set, y = tf_name, colour=gene_set)) +
  theme_classic() + 
  #xlab('') +
  #ylab('') +
  labs(title = "Percent of genes with motif",
       x = "",
       y = "",
       caption = "TSS1kb, qval 0.25, total hits") +
  theme(plot.title = element_text(size = 11)) +
  geom_point(aes(size=percentage), stroke=0) +
  scale_size_continuous(range = c(1.4,8)) +  ## adjust the range of point sizes, ensuring each is visually noticeable
  ylim(sort(unique(perc_df_ss$tf_name), decreasing=TRUE)) + ## sort motif order
  guides(colour='none')
### view output
gp_a
### export plot
ggsaver(
  "bubble_plot-percent_genes_with_motif.main_venn_sets.total.tss1kb.q0.25.wide_version.fimo", 
  gp_a, path=motif_dists_outputs$fig_motif_dists_tss1kb.do, 
  height=cm_to_inch(11), width=cm_to_inch(12), units="in",
  device=c("png", "pdf"))
#### cleanup
rm(gp_a)

