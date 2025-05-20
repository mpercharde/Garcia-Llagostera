# chromatin state metrics

#################################################################################
# load packages
#################################################################################
packages = c("tidyverse", "magrittr", "glue", "here", "data.table", "furrr", 
             "BiocParallel", "BiocManager", "cowplot", "ggsignif", 
             "rstatix", "ggpubr", "ROCR")

for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only=TRUE))
}

rm(p, packages)

#################################################################################
# set functions and more
#################################################################################

# specify output locations
outputs <- list(
  data = "./analysis/5.chip_downstream/mm39/new_dat/chromHMM/initial_eval/new_ptm_tar/state_properties_analysis/data",
  figures = "./analysis/5.chip_downstream/mm39/new_dat/chromHMM/initial_eval/new_ptm_tar/state_properties_analysis/figures",
  results = "./analysis/5.chip_downstream/mm39/new_dat/chromHMM/initial_eval/new_ptm_tar/state_properties_analysis/results",
  oe_dat = "./analysis/5.chip_downstream/mm39/new_dat/chromHMM/initial_eval/new_ptm_tar/overlap_enrichment/results/model_14",
  oe_res = "./analysis/5.chip_downstream/mm39/new_dat/chromHMM/initial_eval/new_ptm_tar/overlap_enrichment/results/model_14/r_analysis_results/results",
  oe_fig = "./analysis/5.chip_downstream/mm39/new_dat/chromHMM/initial_eval/new_ptm_tar/overlap_enrichment/results/model_14/r_analysis_results/figures"
)

# specify org.db
org_db = org.Mm.eg.db

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

pts_to_cm = function(pts) {
  return(pts * 0.0352778)
}

cm_to_inch = function(cm) {
  return(cm / 2.54)
}

#################################################################################
# data import and processing - chromHMM segments
#################################################################################

## read in the bed file
seg <- read.table("./analysis/5.chip_downstream/mm39/new_dat/chromHMM/initial_eval/new_ptm_tar/models/model_14/esc_14_chrHMM.new_ptm_tar_segments.bed", 
                    header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  ### check output and structure
  head(seg)
  tail(seg)
  dim(seg)
## set column names
colnames(seg) <- c("chr", "start", "end", "state")
## filtering contigs
seg <- seg %>%
  filter(chr %in% c("1", "2", "3", "4", "5", "6", 
                    "7", "8", "9", "10", "11", "12", 
                    "13", "14", "15", "16", "17", "18", 
                    "19", "X", "Y"))
## set order of factors for plotting
unique(seg$state)
seg$state <- factor(seg$state, levels = c("E1", "E2", "E3", "E4",
                                          "E5", "E6", "E7", "E8",
                                          "E9", "E10", "E11", "E12",
                                          "E13", "E14"))
  
#################################################################################
# calculate and plot chromHMM domain metrics - number and lengths
#################################################################################
  
# calculate, record, and plot total number of contiguous segments (domains) of each state
d_counts <- seg %>%
  dplyr::count(state)
## process and export the results
### specify column names
colnames(d_counts) <- c("state", "count")
### set order of factors for export and plotting
d_counts$state <- factor(d_counts$state, levels = c("E1", "E2", "E3", "E4",
                                                    "E5", "E6", "E7", "E8",
                                                    "E9", "E10", "E11", "E12",
                                                    "E13", "E14"))
### export the results as csv
write.csv(d_counts, glue("{outputs$results}/state_domain_counts.new_ptm_tar.model_14.csv"), 
          row.names = FALSE)
## plot the results
g <- ggplot(d_counts, aes(x = state, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_classic() +
  labs(title = "Number of segments per state (genome-wide)",
       x = "State", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 14)) +
  scale_y_continuous(expand = c(0, 0)) ## removes padding before 0
g
### export the plot to figures directory
ggsaver(
  "barplot-state_domain_counts.new_ptm_tar.model_14", g, path=outputs$figures,
  height=15, width=16, units="cm")
## cleanup
rm(g)

# calculate, record, and plot the length of each state's segments
## calculate the length of each state and place into a new column
seg <- seg %>%
  mutate(length = end - start)
  ### check the output
  head(seg)
  tail(seg)
  dim(seg)
## calculate the median length of each state's segments and export
### calculate median length
length_df <- seg %>%
  group_by(state) %>%
  summarise(median_length = median(length, na.rm = TRUE))
### check the results
print(length_df)
### export the results as csv
write.csv(length_df, glue("{outputs$results}/state_domain_median_lengths.new_ptm_tar.model_14.csv"), 
          row.names = FALSE)
## plot a boxplot of the lengths of each segment (w/o quiescent (E2, E4))
### remove the quiescent states (E2, E4), keep states of interest
seg_noq <- seg %>%
  filter(!state %in% c("E2", "E4"))
seg_noq$state <- droplevels(seg_noq$state)
### check output
unique(seg_noq$state)
### plot
g <- ggplot(seg_noq, aes(x = state, y = length)) +
  geom_boxplot(fill = "skyblue") +
  labs(x = "State", y = "Segment Length (bp)", title = "Segment lengths per state") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                 size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 14))
### view the plot
g
### export the plot to figures directory
ggsaver(
  "barplot-state_domain_lengths-no_quies.new_ptm_tar.model_14", g, path=outputs$figures,
  height=15, width=16, units="cm")
## cleanup
rm(g)


#################################################################################
# convert counts to tpm, and prepare to extract top expressed genes
#################################################################################

# load in sample metadata
## metadata table and location: `./sample_sheets/deseq-sample_table_all.csv`
metadata <- read_csv("./sample_sheets/deseq-sample_table_all.add_for_ppr.csv")
head(metadata)
## convert column values to factors
metadata <- metadata %>% 
  mutate_all(as.factor)
head(metadata)
## subset metadata
metadata_old <- metadata %>% 
  dplyr::filter(dataset %in% c("old", "old_add")) 
metadata_add <- metadata_old %>% 
  dplyr::filter(treatment != "irf35d") 
metadata_irf <- metadata[metadata$dataset == "old", ]
metadata_new <- metadata[metadata$dataset == "new", ]
### drop excess factor levels
metadata_old <- droplevels(metadata_old)
metadata_add <- droplevels(metadata_add)
metadata_irf <- droplevels(metadata_irf)
metadata_new <- droplevels(metadata_new)
# import and process counts tables
all_counts = 
  metadata %>% 
  dplyr::pull(sample) %>% 
  ## read in each file, marking row names as feature identifiers, defining columns, and combine them into one large df
  furrr::future_map_dfc(function(f) {
    readr::read_tsv(
      here("./analysis/2.telocal.collected", glue("{f}.telocal_out.cntTable")),
      col_names = c("feature", f), col_types = "ci", skip = 1) %>%
      tibble::column_to_rownames("feature") %>%
      purrr::set_names(f)  # Ensure the column name is set to the sample name
  })
rownames(all_counts) = str_remove(rownames(all_counts), ":.*") # cleanup feature names (keeps all before first ':')
## extract, process, and filter raw gene count information
gene_counts <- all_counts[stringr::str_starts(rownames(all_counts), "ENS"), ]
### remove genes with a zero count across all samples
gene_counts <- gene_counts[rowSums(gene_counts) > 0, ]
### filter to only include E14_noIFN
gene_counts <- gene_counts %>% 
  select(E14_noIFN_1, E14_noIFN_2, E14_noIFN_3)

# load in gene lengths information
## read in the bed file
datf <- read.table(glue("{outputs$data}/gene_coords_with_ids_strand_and_lengths.all.tab.bed"), 
                  header = FALSE, sep = "\t", stringsAsFactors = FALSE)
## set column names
colnames(datf) <- c("chr", "start", "end", "strand", "gene_id", "gene_length")
  ### check output
  head(datf)
  tail(datf)
  dim(datf)
### save to seg for processing
seg <- datf
  ### check output and structure
  head(seg)
  tail(seg)
  dim(seg)
## filter to only include chromosomes of interest
unique(seg$chr)
seg <- seg %>%
  filter(chr %in% c("1", "2", "3", "4", "5", "6", 
                    "7", "8", "9", "10", "11", "12", 
                    "13", "14", "15", "16", "17", "18", 
                    "19", "X", "Y"))
  ### check output
  head(seg)
  dim(seg)
  unique(seg$chr)
## filter genes in the counts_df to only include those in the chromosomes of interest
genes_chr <- seg$gene_id
gene_counts <- gene_counts[rownames(gene_counts) %in% genes_chr, ]

# calculate TPM
## convert counts to matrix
counts_mat <- as.matrix(gene_counts)
## filter the gene lengths to only include those in the counts matrix and keep the same order too
### extract gene vector
genes_in_counts <- rownames(counts_mat)
### filter gene lengths to only include those genes
seg <- seg[seg$gene_id %in% genes_in_counts, ]
### reorder to match the counts matrix
seg <- seg[match(genes_in_counts, seg$gene_id), ]
### sanity check
all.equal(seg$gene_id, genes_in_counts) 
## extract gene lengths as a vector
lens <- seg$gene_length
## calculate tpm (Michael Love's method, see https://www.biostars.org/p/335187/ and https://support.bioconductor.org/p/91218/#91256 )
x <- counts_mat / lens
tpm.mat <- t( t(x) * 1e6 / colSums(x) )

# define the mean TPM and record in a new column, and add gene information and filter
## convert to df
tpm.df <- as.data.frame(tpm.mat)
## calculate mean and place into a new column
tpm.df$mean_tpm <- rowMeans(tpm.df[, c("E14_noIFN_1", "E14_noIFN_2", "E14_noIFN_3")])
## add gene_id rownames into a new column too
tpm.df$gene_id <- rownames(tpm.df)
## add gene annotations
tpm.df = annotate_genes(tpm.df, org_db)

# extract top 20% expressed genes
## sort from largest to smallest mean tpm value
tpm.df <- tpm.df[order(tpm.df$mean_tpm, decreasing = TRUE), ]
## extract top 20% expressed genes
thresh <- quantile(tpm.df$mean_tpm, probs = 0.80)
tpm_top20 <- tpm.df[tpm.df$mean_tpm >= thresh, ]
### check output
head(tpm_top20)
tail(tpm_top20)
dim(tpm_top20)

# process coordinates to calculate +/- 2kb around TSS, and export for OE analysis
## filter to only include genes in top20
### save new data for processing
dat <- datf
### filter columns to remove length
dat <- dat %>% 
  dplyr::select(!gene_length)
### make a vector of the top 20% of genes
top20_g <- tpm_top20$gene_id
### filter to include only those top20 genes
dat <- dat[dat$gene_id %in% top20_g, ]
  ### check output
  dim(dat)
  length(top20_g)
  head(dat)
## make a TSS column (extracting 'start' if + and 'end' value if -)
dat$tss <- ifelse(dat$strand == "+", dat$start, dat$end)
  ### check the output
  head(dat)
  tail(dat)
  ### sanity check
  table(dat$strand, dat$tss == dat$start)
## calculate + and - 2kb around TSS
### calculate +2kb around TSS
dat$tss_plus2kb <- dat$tss + 2000
### calculate -2kb around TSS
dat$tss_min2kb <- dat$tss - 2000
  ### check output
  head(dat)
  tail(dat)
## set up bed file of +/-2kb around TSS of highly expressed genes for export
dat <- dat %>% 
  select(chr, tss_min2kb, tss_plus2kb, gene_id)
  ### check output
  head(dat)
  tail(dat)
  dim(dat)
## export the results as bed file
write.table(dat, file = glue("{outputs$results}/top20p_exp_genes_coords.e14_noIFN.bed"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## cleanup
rm(metadata, metadata_add, metadata_irf, metadata_new, metadata_old, seg, x, 
   all_counts, counts_mat, dat, lens, thresh, genes_in_counts, genes_chr,
   tpm.mat, tpm.df)

#################################################################################
# transforming (log10) and plotting OE analysis of gene sets
#################################################################################

## load in fold enrichment values from overlap enrichment
res <- read.table(glue("{oe_outputs$oe_dat}/model_14-main_venn_sets-2kbTSS.no_center.txt"), header = T, 
                  row.names = 1, ## setting states column as rownames
                  sep = "\t")
  ### check output
  head(res)
  tail(res)
  colnames(res)

## process data frame
### rename columns
#### create a vector of column names
colnames(res)
col_names <- c("BV", "C", "IMI", "IR", "MI", "RAN.n1106", "RAN.n123", "RAN.n2408")
colnames(res) <- col_names
### remove base row (dont need it for log2)
rownames(res)
res_f <- res[!(rownames(res) %in% "Base"), ]
### preparing with only main venn sets
res_v <- res_f %>% 
  select(-BV, -RAN.n1106, 
         -RAN.n123, -RAN.n2408)
#### set order of levels for display
res_v <- res_v[,c(1,2,3,4)]
#### convert the dataframe to a numeric matrix
mat_v <- as.matrix(res_v, mode = "numeric")
  #### check if numeric
  is.numeric(mat_v)
  class(mat_v)
### preparing for main venn sets with all random sets
res_r <- res_f %>% 
  select(-BV)
#### set order of levels for display
res_r <- res_r[,c(1,2,3,4,6,5,7)]
#### convert the dataframe to a numeric matrix
mat_r <- as.matrix(res_r, mode = "numeric")
  #### check if numeric
  is.numeric(mat_r)
  class(mat_r)

## perform log10 transformation of data and export - gene sets
  ##mat_v_log <- log10(mat_v + 1e-8)  ## pseudocount
mat_v_log <- log10(mat_v)
### check output
summary(mat_v_log)
any(is.na(mat_v_log))
## convert to dataframe
df.v.log <- as.data.frame(mat_v_log)
### add chromatin state info
df.v.log$state <- rownames(mat_v_log)
## export results as csv
write.csv(df.v.log, glue("{outputs$oe_res}/oe-main_venn_sets-log10_trans_vals.no_center.new_ptm_tar.model_14.csv"), 
          row.names = FALSE)

## plot heatmap of log10 transformed fold enrichment values - gene sets
oeCol_r = circlize::colorRamp2(c(-1.6, 0, 1.6), 
                               c("blue", "white", "red"))
### plot without number labels
pl1_v_log <- ComplexHeatmap::Heatmap(mat_v_log, 
                                     name = "Log10 FE", 
                                     col = oeCol_r, 
                                     cluster_rows = F, 
                                     show_row_dend = FALSE,
                                     cluster_columns = FALSE,
                                     row_names_side = "left",
                                     rect_gp = gpar(col = "darkgrey", lwd = 0.1),
                                     column_title = "Log2FE of states in gene sets (w/ random genes, no_center)")
pl1_v_log
#### export the plot
pdf(here(outputs$oe_fig, glue("oe-main_venn_sets-log10_trans_vals.no_center.new_ptm_tar.model_14.no_num.pdf")),
    width = cm_to_inch(13), height = cm_to_inch(22)
)
draw(pl1_v_log)
dev.off()
png(here(outputs$oe_fig, glue("oe-main_venn_sets-log10_trans_vals.no_center.new_ptm_tar.model_14.no_num.png")),
    width = 13, height = 22, units = "cm", res = 150
)
draw(pl1_v_log)
dev.off()
### plot with number value labels
pl2_v_log <- ComplexHeatmap::Heatmap(mat_v_log, 
                                     name = "Log10 FE", 
                                     col = oeCol_r, 
                                     cluster_rows = F, 
                                     show_row_dend = FALSE,
                                     cluster_columns = FALSE,
                                     row_names_side = "left",
                                     rect_gp = gpar(col = "darkgrey", lwd = 0.1),
                                     column_title = "Log10FE of states in gene sets (no_center)",
                                     cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.text(sprintf("%.1f", mat_v_log[i, j]), x, y, gp = gpar(fontsize = 10))}
)
pl2_v_log
#### export the plot
pdf(here(outputs$oe_fig, glue("oe-main_venn_sets-log10_trans_vals.no_center.new_ptm_tar.model_14.num.pdf")),
    width = cm_to_inch(13), height = cm_to_inch(22)
)
draw(pl2_v_log)
dev.off()
png(here(outputs$oe_fig, glue("oe-main_venn_sets-log10_trans_vals.no_center.new_ptm_tar.model_14.num.png")),
    width = 13, height = 22, units = "cm", res = 150
)
draw(pl2_v_log)
dev.off()

