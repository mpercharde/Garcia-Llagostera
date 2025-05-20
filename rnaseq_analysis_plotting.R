# RNAseq analysis and plotting

#################################################################################
# load packages
#################################################################################

packages = c("tidyverse", "magrittr", "glue", "here", "data.table", "furrr", 
             "BiocParallel", "BiocManager", "knitr", "readr", "kableExtra", 
             "DT", "pheatmap", "RColorBrewer", "GenomicFeatures", "DESeq2", 
             "PCAtools", "IHW", "ashr", "variancePartition", "AnnotationDbi", 
             "org.Mm.eg.db", "viridis", "RUVSeq", "ComplexHeatmap", "Cairo", 
             "clusterProfiler", "eulerr", "openxlsx", "biomaRt", "BSgenome", 
             "BSgenome.Mmusculus.UCSC.mm39", "memes", "universalmotif", 
             "MotifDb", "ggtree", "motifStack", "cowplot", "ggsignif", 
             "rstatix", "ggpubr", "ROCR", "enrichplot", "fgsea", "ggbreak")

for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only=TRUE))
}

rm(p, packages)

#################################################################################
# set functions and more
#################################################################################

# specify output locations
outputs <- list(
  figures = "./analysis/diff_exp/figures", 
  results_g = "./analysis/diff_exp/results/genes",
  res_contrasts_g = "./analysis/diff_exp/results/genes/contrasts",
  res_go_g = "./analysis/diff_exp/results/genes/go"
)

# specify org.db
org_db = org.Mm.eg.db

# specify statistical thresholds for differential expression testing
stats <- list(
  de_alpha = 0.01,  # padj threshold for DESeq2
  de_l2FC_old = 0,  # log2 fold change threshold: old dataset
  de_l2FC_E14 = 0,  # log2 fold change threshold: E14 (mESCs)
  de_l2FC_MEF = 0,  # log2 fold change threshold: MEFs
  enrich_alpha = 0.05,  # q-value for GO analysis w/ clusterProfiler
  go_minimum = 15  # minimum number of genes required to pass gene set onto clusterProfiler for GO analysis
)

# stylistic parameters for plotting
hCol = colorRampPalette(rev(brewer.pal(11, "RdBu"))) 

# functions and other parameters/settings

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

annotate_genes = function(
    x, org_db,
    columns=c("gene_symbol"="SYMBOL", "gene_name"="GENENAME"),
    keycol="gene_id", keytype="ENSEMBL") {
  # wraps AnnotationDb::mapIds to add gene annotations to a table
  # keycol can be either a column name within `x` or "rownames"
  if (keycol %ni% colnames(x) && keycol != "rownames") stop(glue::glue("Column '{keycol}' not in table"))
  if (keytype %ni% columns(org_db)) stop(glue::glue("Column type '{keytype}' not present within provided database"))
  if (!is.data.frame(x)) x = data.frame(x)
  if (keycol == "rownames") {
    keys = rownames(x)
  } else {
    keys = x %>% dplyr::pull(var=c(glue::glue("{keycol}")))
  }
  purrr::map2(
    names(columns), columns,
    function(name, column) {
      x[[name]] <<- unname(AnnotationDbi::mapIds(org_db, keys=keys, column=column, keytype=keytype, multiVals="first"))
    }
  )
  x
}

contrast_and_save = function(
    contrast_name, feature_type, ## feature_type = genes, teloci, tefam, etc.
    output_location_tab, output_location_fig,  ## output_location_tab = contrasts tsv output location ; output_location_fig = MAplot output location
    filter_result = T, plot_result=T, write_table=T, add_names=T, ...) {
  ## run and return the output of DESeq2::results for a given contrast and 
  ## parameter set, additionally filtering and saving the results table
  
  args = c(as.list(environment()), list(...))
  args$contrast_name = NULL
  args$feature_type = NULL
  args$output_location_tab = NULL
  args$output_location_fig = NULL
  args$filter_result = NULL
  args$plot_result = NULL
  args$write_table = NULL
  args$add_names = NULL
  args$vst_obj = NULL ## use to set vst separately for sample subsets
  
  res = rlang::exec(DESeq2::results, !!!args, filterFun=IHW::ihw)
  res = DESeq2::lfcShrink(dds=args$object, res=res, type="ashr")
  
  if (plot_result) {
    ma_data = DESeq2::plotMA(res, alpha=args$alpha, returnData=T)
    y_limit = ceiling(max(abs(quantile(ma_data$lfc, c(0.001, 0.999)))))
    g = 
      ggplot2::ggplot(ma_data, aes(mean, lfc, color=isDE)) +
      ggplot2::geom_point(alpha=0.33, size=1, stroke=0) +
      ggplot2::geom_hline(yintercept=c(-args$lfcThreshold, args$lfcThreshold)) +
      ggplot2::scale_color_manual(values=c(`TRUE`="blue", `FALSE`="grey75")) +
      ggplot2::coord_cartesian(ylim=c(-y_limit, y_limit)) +  ## computes the 0.001 and 0.999 quantiles of the lfc values, returning the values that correspond to the 0.1% and 99.9% percentiles of the data
      ggplot2::scale_x_continuous(trans="log10") +   ## apply log10 transformation to mean of normalized counts
      ggplot2::labs(
        x="Mean of normalized counts",
        y="Shrunken log2 Fold-Change",
        title=glue("{contrast_name} ({feature_type})"),
        caption=glue('DESeq results: padj {args$alpha}, logFC {args$lfcThreshold}')) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="none")
    # plot the result
    plot(g)
    # export
    ggsaver(
      glue("de-MAplot_{feature_type}-{contrast_name}.add_for_ppr"), g, path=output_location_fig, 
      height=cm_to_inch(10), width=cm_to_inch(10), units="in",
      device=c("png", "pdf"))
  }
  
  res = data.frame(res)
  res_all <- res ## extract a version with all, not just significant, results
  if (filter_result) {
    res = dplyr::filter(res, padj < args$alpha, abs(log2FoldChange) >= args$lfcThreshold)
  }
  res = 
    res %>%
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    tibble::rownames_to_column(var="gene_id")  %>%
    dplyr::select(gene_id, log2FoldChange, lfcSE, padj)
  
  res_all = res_all %>%
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    tibble::rownames_to_column(var="gene_id") %>%
    dplyr::select(gene_id, log2FoldChange, lfcSE, padj)
  
  if (add_names) {
    res = annotate_genes(res, org_db)
    res_all = annotate_genes(res_all, org_db)
  }
  
  if (write_table) readr::write_tsv(
    res, here(output_location_tab, glue("de_{feature_type}-{contrast_name}.sig_res.add_for_ppr.tsv")))
  
  if (write_table) readr::write_tsv(
    res_all, here(output_location_tab, glue("de_{feature_type}-{contrast_name}.all_res.add_for_ppr.tsv")))
  
  return(res)
}

pts_to_cm = function(pts) {
  return(pts * 0.0352778)
}

cm_to_inch = function(cm) {
  return(cm / 2.54)
}

# calculate covariance
co.var <- function(x) {
  return(100 * sd(x) / mean(x))  ## covariance = (sd/mean) * 100
}


#################################################################################
# data import and processing
#################################################################################

# load in sample metadata
metadata <- read_csv("./sample_sheets/deseq-sample_table_all.csv")
metadata <- metadata %>% 
  mutate_all(as.factor)
## subset metadata
metadata_old <- metadata %>% 
  dplyr::filter(dataset %in% c("old", "old_add"))
metadata_add <- metadata_old %>% 
  dplyr::filter(treatment != "irf35d") ## subset metadata to only contain recently added old datasets (aka sting gof, sting lof, and the bfp control)
metadata_irf <- metadata[metadata$dataset == "old", ]  ## subset metadata for old data already analyzed (aka irf35d + bfp only)
metadata_new <- metadata[metadata$dataset == "new", ]  ## subset metadata for new datasets
### drop excess factor levels
metadata_old <- droplevels(metadata_old)
metadata_add <- droplevels(metadata_add)
metadata_irf <- droplevels(metadata_irf)
metadata_new <- droplevels(metadata_new)


# import and process TElocal counts tables
  ## lower expression cutoff thresholds:
    ## Genes = removing genes with a zero count across all samples
    ## TE locus-level = keeping TEs with a non-zero count in \>=3 samples
    ## TE subfamily = removing subfamilies with a zero count across all samples

## read in all TElocal counts tables, and stitch them together into one df (rows as features)
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
## check the output
head(all_counts)
tail(all_counts)
dim(all_counts)

## extract and process gene information
dds_genes = 
    DESeq2::DESeqDataSetFromMatrix( # create a deseq dataset only containing genes
      all_counts[stringr::str_starts(rownames(all_counts), "ENS"), ],
      metadata,
      ~0) # placeholder design, no specific design set yet
  # check output
  #dds_genes
### remove genes with a zero count across all samples
dds_genes = dds_genes[rowSums(assay(dds_genes)) > 0, ]
  # check output
  #dds_genes
    ## from 56941 to 38748 rows (68.1% maintained)
dds_genes = DESeq2::estimateSizeFactors(dds_genes) # estimate size factors
vsd_genes = assay(DESeq2::vst(dds_genes)) # VST
### extract gene information specifically for each dataset (running each DESeq modeling and analysis separately due to differences in within-group-variability)
  ### subset counts data
  old_counts <- all_counts[, colnames(all_counts) %in% metadata_old$sample]
  add_counts <- all_counts[, colnames(all_counts) %in% metadata_add$sample]
  irf_counts <- all_counts[, colnames(all_counts) %in% metadata_irf$sample]
  new_counts <- all_counts[, colnames(all_counts) %in% metadata_new$sample]
  ### extract and process gene information for old datasets
    dds_genes_old =  ### prepare dds_genes to only include relevant samples
      DESeq2::DESeqDataSetFromMatrix(  # create a deseq dataset only containing genes from relevant samples
        old_counts[stringr::str_starts(rownames(old_counts), "ENS"), ],
        metadata_old,
        ~0) # placeholder design, no specific design set yet
    dds_genes_old = DESeq2::estimateSizeFactors(dds_genes_old) # estimate size factors
      ### double check output
      dds_genes_old@colData@rownames  ## all good
    vsd_genes_old = assay(DESeq2::vst(dds_genes_old)) # VST of subset
  ### extract and process gene information for added (old_add) datasets
    dds_genes_add =  ### prepare dds_genes to only include relevant samples
      DESeq2::DESeqDataSetFromMatrix(  # create a deseq dataset only containing genes from relevant samples
        add_counts[stringr::str_starts(rownames(add_counts), "ENS"), ],
        metadata_add,
        ~0) # placeholder design, no specific design set yet
    dds_genes_add = DESeq2::estimateSizeFactors(dds_genes_add) # estimate size factors
    ### double check output
    dds_genes_add@colData@rownames  ## all good
    vsd_genes_add = assay(DESeq2::vst(dds_genes_add)) # VST of subset
  ### extract and process gene information for IRF3 datasets
    ### subset dds_genes to only include relevant samples
    dds_genes_irf = 
      DESeq2::DESeqDataSetFromMatrix(  # create a deseq dataset only containing genes from relevant samples
        irf_counts[stringr::str_starts(rownames(irf_counts), "ENS"), ],
        metadata_irf,
        ~0) # placeholder design, no specific design set yet
    dds_genes_irf = DESeq2::estimateSizeFactors(dds_genes_irf) # estimate size factors
    ### double check output
    dds_genes_irf@colData@rownames  ## all good
    vsd_genes_irf = assay(DESeq2::vst(dds_genes_irf)) # VST of subset
  ### extract and process gene information for new datasets
    ### subset dds_genes to only include relevant samples
    dds_genes_new = 
      DESeq2::DESeqDataSetFromMatrix(  # create a deseq dataset only containing genes from relevant samples
        new_counts[stringr::str_starts(rownames(new_counts), "ENS"), ],
        metadata_new,
        ~0) # placeholder design, no specific design set yet
    dds_genes_new = DESeq2::estimateSizeFactors(dds_genes_new) # estimate size factors
      ### double check output
      dds_genes_new@colData@rownames  ## all good
    vsd_genes_new = assay(DESeq2::vst(dds_genes_new)) # VST of subset

## extract and process TE information (locus-specific)
dds_teloci = 
    DESeq2::DESeqDataSetFromMatrix( # create a deseq dataset only containing TEs
      all_counts[!stringr::str_starts(rownames(all_counts), "ENS"), ], # remove genes
      metadata,
      ~0) # placeholder design, no specific design set yet
  # check output
  #dds_teloci
### remove low-count TEs, keeping those with a non-zero count in >=3 samples
dds_teloci = dds_teloci[rowSums(counts(dds_teloci)>0) >= 3, ]
  # check output
  #dds_teloci
    ## from 4467522 to 334623 rows (7.5% maintained)
dds_teloci = DESeq2::estimateSizeFactors(dds_teloci) # estimate size factors
### extract TE information (locus-specific) separately for each dataset
  ### extract and process TE information (locus-specific) for old datasets
    dds_teloci_old =  ### prepare dds_teloci to only include relevant samples
      DESeq2::DESeqDataSetFromMatrix(  # create a deseq dataset only containing teloci from relevant samples
        old_counts[!stringr::str_starts(rownames(old_counts), "ENS"), ],
        metadata_old,
        ~0) # placeholder design, no specific design set yet
    ### remove low-count TEs, keeping those with a non-zero count in >=3 samples
    dds_teloci_old = dds_teloci_old[rowSums(counts(dds_teloci_old)>0) >= 3, ]
    dds_teloci_old = DESeq2::estimateSizeFactors(dds_teloci_old) # estimate size factors
  ### extract and process TE information (locus-specific) for new datasets
    ### subset dds_teloci to only include relevant samples
    dds_teloci_new =  ### prepare dds_teloci to only include relevant samples
      DESeq2::DESeqDataSetFromMatrix(  # create a deseq dataset only containing teloci from relevant samples
        new_counts[!stringr::str_starts(rownames(new_counts), "ENS"), ],
        metadata_new,
        ~0) # placeholder design, no specific design set yet
    ### remove low-count TEs, keeping those with a non-zero count in >=3 samples
    dds_teloci_new = dds_teloci_new[rowSums(counts(dds_teloci_new)>0) >= 3, ]
    dds_teloci_new = DESeq2::estimateSizeFactors(dds_teloci_new) # estimate size factors

## extract and process TE information (subfamily level)
dds_tefam = 
    DESeq2::DESeqDataSetFromMatrix( # create a deseq dataset for TE subfamilies
      all_counts %>% 
        tibble::rownames_to_column(var="feature") %>%  # move rownames into a column
        dplyr::filter(!stringr::str_starts(feature, "ENS")) %>%  # remove genes
        tidyr::separate_wider_delim(feature, names=c("family", "element", "chrom", "start", "stop"), delim="|") %>% # split the feature names into their specific elements
        dplyr::mutate(te_type=glue("{family}_{element}")) %>%  # combine family and element and move into a new column
        dplyr::select(-c("family", "element", "chrom", "start", "stop")) %>%  # remove unnecessary columns, keeping only subfamily TE info
        dplyr::group_by(te_type) %>%  # group by TE type, for subsequent summarization
        dplyr::summarise(across(everything(), sum)) %>%  # for each subfamily, takes the sum of the counts, thereby summarizing the counts for each TE type across samples
        tibble::column_to_rownames(var="te_type"), # convert te_type column to rownames
      metadata,
      ~0)
  # check output
  #dds_tefam
### remove subfamilies with a zero count across all samples
dds_tefam = dds_tefam[rowSums(assay(dds_tefam)) > 0, ]
  # check output
  #dds_tefam
    ## from 1238 to 1218 rows (98.3% maintained)
dds_tefam = DESeq2::estimateSizeFactors(dds_tefam)
### extract TE information (subfamily level) separately for each dataset
  ### extract and process TE information (subfamily level) for old datasets
    dds_tefam_old =  ### prepare dds_tefam to only include relevant samples
         DESeq2::DESeqDataSetFromMatrix( # create a deseq dataset for TE subfamilies
            old_counts %>% 
              tibble::rownames_to_column(var="feature") %>%  # move rownames into a column
              dplyr::filter(!stringr::str_starts(feature, "ENS")) %>%  # remove genes
              tidyr::separate_wider_delim(feature, names=c("family", "element", "chrom", "start", "stop"), delim="|") %>% # split the feature names into their specific elements
              dplyr::mutate(te_type=glue("{family}_{element}")) %>%  # combine family and element and move into a new column
              dplyr::select(-c("family", "element", "chrom", "start", "stop")) %>%  # remove unnecessary columns, keeping only subfamily TE info
              dplyr::group_by(te_type) %>%  # group by TE type, for subsequent summarization
              dplyr::summarise(across(everything(), sum)) %>%  # for each subfamily, takes the sum of the counts, thereby summarizing the counts for each TE type across samples
              tibble::column_to_rownames(var="te_type"), # convert te_type column to rownames
            metadata_old,
            ~0)  # placeholder design, no specific design set yet
    ## remove subfamilies with a zero count across all samples
    dds_tefam_old = dds_tefam_old[rowSums(assay(dds_tefam_old)) > 0, ]
    dds_tefam_old = DESeq2::estimateSizeFactors(dds_tefam_old) # estimate size factors
  ### extract and process TE information (subfamily level) for new datasets
    ### subset dds_tefam to only include relevant samples
    dds_tefam_new =  ### prepare dds_tefam to only include relevant samples
         DESeq2::DESeqDataSetFromMatrix( # create a deseq dataset for TE subfamilies
            new_counts %>% 
              tibble::rownames_to_column(var="feature") %>%  # move rownames into a column
              dplyr::filter(!stringr::str_starts(feature, "ENS")) %>%  # remove genes
              tidyr::separate_wider_delim(feature, names=c("family", "element", "chrom", "start", "stop"), delim="|") %>% # split the feature names into their specific elements
              dplyr::mutate(te_type=glue("{family}_{element}")) %>%  # combine family and element and move into a new column
              dplyr::select(-c("family", "element", "chrom", "start", "stop")) %>%  # remove unnecessary columns, keeping only subfamily TE info
              dplyr::group_by(te_type) %>%  # group by TE type, for subsequent summarization
              dplyr::summarise(across(everything(), sum)) %>%  # for each subfamily, takes the sum of the counts, thereby summarizing the counts for each TE type across samples
              tibble::column_to_rownames(var="te_type"), # convert te_type column to rownames
            metadata_new,
            ~0)  # placehnewer design, no specific design set yet
    ## remove subfamilies with a zero count across all samples
    dds_tefam_new = dds_tefam_new[rowSums(assay(dds_tefam_new)) > 0, ]
    dds_tefam_new = DESeq2::estimateSizeFactors(dds_tefam_new) # estimate size factors

## cleanup
rm(old_counts, add_counts, irf_counts, new_counts)


#################################################################################
# initial quality assessment
#################################################################################

# plotting barplot of the number of genes recovered in each sample
## plotting no. of genes recovered in each sample
g = ggplot2::ggplot(
  metadata %>% 
    dplyr::mutate(
      sample=forcats::as_factor(sample), # ensure sample is a factor variable
      num_genes=colSums(assay(dds_genes)>10)), # determine the no. of genes with a count > 10 in each sample
  ggplot2::aes(x=sample, y=num_genes)) +
  ggplot2::geom_col() +
  ggplot2::ylab("Number of genes recovered") +
  ggplot2::theme_bw() + 
  ggplot2::theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
## view + assess the plot
g
## save the output to figures directory
ggsaver("num_genes_recovered.add_for_ppr", g, path=outputs$figures,
        height=10, width=12, units="cm")
## cleanup
rm(g)

# plotting log2-transformed gene expression
  ## labeling by treatment
## plotting
g = ggplot2::ggplot(
      log2(counts(dds_genes, normalized=FALSE)+1) %>%  # extract raw, unnormalized counts and perform log2 transformation
        data.table() %>%
        tidyr::gather(key="sample", value="log2_counts", factor_key=T) %>%
        dplyr::left_join(., metadata %>% dplyr::select(sample, treatment)),
      ggplot2::aes(x=log2_counts, group=sample, color=treatment)) +
    ggplot2::geom_density(show.legend = TRUE) +
    ggplot2::scale_colour_viridis_d(option="turbo", na.value="grey") +
    ggplot2::xlab("log2 gene expression") +
    ggplot2::theme_bw()
## view + assess the plot
g
## save the output to figures directory
ggsaver("log2_gene_expression-no_norm-treatment.add_for_ppr", g, path=outputs$figures,
    height=15, width=20, units="cm")
## labeling by dataset
g = ggplot2::ggplot(
      log2(counts(dds_genes, normalized=FALSE)+1) %>%  # extract raw, unnormalized counts and perform log2 transformation
        data.table() %>%
        tidyr::gather(key="sample", value="log2_counts", factor_key=T) %>%
        dplyr::left_join(., metadata %>% dplyr::select(sample, dataset)),
      ggplot2::aes(x=log2_counts, group=sample, color=dataset)) +
    ggplot2::geom_density(show.legend = TRUE) +
    #ggplot2::scale_colour_viridis_d(option="turbo", na.value="grey") +
    ggplot2::xlab("log2 gene expression") +
    ggplot2::theme_bw()
## view + assess the plot
g
## save the output to figures directory
ggsaver("log2_gene_expression-no_norm-dataset.add_for_ppr", g, path=outputs$figures,
    height=15, width=20, units="cm")
## cleanup
rm(g)

# note: log2 gene expression plots display a shifted curve for old_add vs irf3 and bfp data, it's intermediate between irf3 and new (some irf3 appear to still be in the curve though) (18.3.25)

#################################################################################
# assessing correlation patterns - plotting PCAs and heatmaps of euclidean distances (vst normalized, pre-batch correction)
#################################################################################

# plotting PCAs of vst normalized data (pre-batch correction)

## plotting PCAs (PCs 1-4) for all samples (VST normalized), coloring by treatment
### perform PCA (all samples)
pca = PCAtools::pca(
  vsd_genes,
  tibble::column_to_rownames(metadata, "sample"),
  removeVar=0.1)  # remove the lower 10% of variables based on variance
### prepare a df with PCs 1 and 2
pcadf = data.frame(PC1=pca$rotated$PC1, PC2=pca$rotated$PC2, pca$metadata)
### prepare one with PCs 3 and 4
pcadf2 = data.frame(PC3=pca$rotated$PC3, PC4=pca$rotated$PC4, pca$metadata)
### generate PCA plot of all samples, labeled by treatment (NOT adding sample labels to the plots)
#### PC1 + PC2
g = ggplot2::ggplot(pcadf, aes(PC1, PC2, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata$sample), show.legend=FALSE, max.overlaps=nrow(metadata)) +
  ggplot2::xlab(glue('PC1: {round(pca$variance["PC1"])}% variance')) +
  ggplot2::ylab(glue('PC2: {round(pca$variance["PC2"])}% variance')) +
  ggplot2::labs(title = "PCA of all samples by treatment", caption = "VST normalized counts") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
##### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-all_data_treatment.pc1_2.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
#### PC3 + PC4
g = ggplot2::ggplot(pcadf2, aes(PC3, PC4, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata$sample), show.legend=FALSE, max.overlaps=nrow(metadata)) +
  ggplot2::xlab(glue('PC3: {round(pca$variance["PC3"])}% variance')) +
  ggplot2::ylab(glue('PC4: {round(pca$variance["PC4"])}% variance')) +
  ggplot2::labs(title = "PCA of all samples by treatment", caption = "VST normalized counts") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
##### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-all_data_treatment.pc3_4.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")

## generate PCA plot of all samples, labeled by replicate
### PC1 + PC2
g = ggplot2::ggplot(pcadf, aes(PC1, PC2, color=replicate)) +
  ggplot2::geom_point(size=2) +
  ggrepel::geom_text_repel(mapping=aes(label=metadata$sample), show.legend=FALSE, max.overlaps=nrow(metadata)) +
  ggplot2::xlab(glue('PC1: {round(pca$variance["PC1"])}% variance')) +
  ggplot2::ylab(glue('PC2: {round(pca$variance["PC2"])}% variance')) +
  ggplot2::labs(title = "PCA of all samples by replicate", caption = "VST normalized counts") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
### save the output to figures directory
ggsaver(
  "pca_genes-all_data_replicate.pc1_2.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
### PC3 + PC4
g = ggplot2::ggplot(pcadf2, aes(PC3, PC4, color=replicate)) +
  ggplot2::geom_point(size=2) +
  ggrepel::geom_text_repel(mapping=aes(label=metadata$sample), show.legend=FALSE, max.overlaps=nrow(metadata)) +
  ggplot2::xlab(glue('PC3: {round(pca$variance["PC3"])}% variance')) +
  ggplot2::ylab(glue('PC4: {round(pca$variance["PC4"])}% variance')) +
  ggplot2::labs(title = "PCA of all samples by replicate", caption = "VST normalized counts") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
### save the output to figures directory
ggsaver(
  "pca_genes-all_data_replicate.pc3_4.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
### cleanup
rm(pca, pcadf, pcadf2, g)

## perform PCA separately for each dataset, plotting for old data only
### using vst performed separately on each dataset
### perform PCA ('old' samples)
pca = PCAtools::pca(
  vsd_genes_old,
  tibble::column_to_rownames(metadata_old, "sample"), # use filtered metadata
  removeVar=0.1)
# prepare a df with PCs 1 and 2 ('old' samples)
pcadf = data.frame(PC1=pca$rotated$PC1, PC2=pca$rotated$PC2, pca$metadata)
# prepare one with PCs 3 and 4 ('old' samples)
pcadf2 = data.frame(PC3=pca$rotated$PC3, PC4=pca$rotated$PC4, pca$metadata)
### generate PCA plot of 'old' samples, labeled by treatment
#### PC1 + PC2
g = ggplot2::ggplot(pcadf, aes(PC1, PC2, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata_old$sample), show.legend=FALSE, max.overlaps=nrow(metadata_old)) +
  ggplot2::xlab(glue('PC1: {round(pca$variance["PC1"])}% variance')) +
  ggplot2::ylab(glue('PC2: {round(pca$variance["PC2"])}% variance')) +
  ggplot2::labs(title = "PCA of original dataset (old) samples by treatment", caption = "VST normalized counts") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  #ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-old_data_treatment.pc1_2.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
#### PC3 + PC4
g = ggplot2::ggplot(pcadf2, aes(PC3, PC4, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata_old$sample), show.legend=FALSE, max.overlaps=nrow(metadata_old)) +
  ggplot2::xlab(glue('PC3: {round(pca$variance["PC3"])}% variance')) +
  ggplot2::ylab(glue('PC4: {round(pca$variance["PC4"])}% variance')) +
  ggplot2::labs(title = "PCA of original dataset (old) samples by treatment", caption = "VST normalized counts") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  #ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-old_data_treatment.pc3_4.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
### cleanup
rm(pca, pcadf, pcadf2, g)

## perform PCA separately for each dataset, plotting for old_add data only
### using vst performed separately on each dataset
### perform PCA ('old_add' samples)
pca = PCAtools::pca(
  vsd_genes_add,
  tibble::column_to_rownames(metadata_add, "sample"), # use filtered metadata
  removeVar=0.1)
# prepare a df with PCs 1 and 2 ('old_add' samples)
pcadf = data.frame(PC1=pca$rotated$PC1, PC2=pca$rotated$PC2, pca$metadata)
# prepare one with PCs 3 and 4 ('old_add' samples)
pcadf2 = data.frame(PC3=pca$rotated$PC3, PC4=pca$rotated$PC4, pca$metadata)
### generate PCA plot of 'old_add' samples, labeled by treatment
#### PC1 + PC2
g = ggplot2::ggplot(pcadf, aes(PC1, PC2, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata_add$sample), show.legend=FALSE, max.overlaps=nrow(metadata_add)) +
  ggplot2::xlab(glue('PC1: {round(pca$variance["PC1"])}% variance')) +
  ggplot2::ylab(glue('PC2: {round(pca$variance["PC2"])}% variance')) +
  ggplot2::labs(title = "PCA of STING dataset (old_add) samples by treatment", caption = "VST normalized counts") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  #ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-sting_data_treatment.pc1_2.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
#### PC3 + PC4
g = ggplot2::ggplot(pcadf2, aes(PC3, PC4, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata_add$sample), show.legend=FALSE, max.overlaps=nrow(metadata_add)) +
  ggplot2::xlab(glue('PC3: {round(pca$variance["PC3"])}% variance')) +
  ggplot2::ylab(glue('PC4: {round(pca$variance["PC4"])}% variance')) +
  ggplot2::labs(title = "PCA of STING dataset (old_add) samples by treatment", caption = "VST normalized counts") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  #ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-sting_data_treatment.pc3_4.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
### cleanup
rm(pca, pcadf, pcadf2, g)


# plotting an annotated heatmap of pairwise (euclidean) distances between samples (vst normalized)
## compute sample-to-sample distances (pairwise Euclidean distances)
dists = dist(t(vsd_genes))
dist_mat = as.matrix(dists) ## convert to matrix
## process matrix, assign sample names
colnames(dist_mat) = dds_genes$sample
## check it
head(dist_mat)
## prepare heatmap annotation from metadata
annot = 
  metadata %>% 
  tibble::column_to_rownames("sample") %>%
  dplyr::select(treatment, cell_type, replicate) ## choose relevant columns
## generate the annotated heatmap containing all samples
g = pheatmap(
  dist_mat,
  annotation=annot, ## specify the annotation
  ## set annotation colors
  annotation_colors =
    list(
      "Treatment"=
        viridis::turbo(length(levels(metadata$treatment))) %>% 
        magrittr::set_names(levels(metadata$treatment)),
      "Cell Type"=
        viridis::turbo(length(levels(metadata$cell_type))) %>% 
        magrittr::set_names(levels(metadata$cell_type)),
      "Replicate"=
        viridis::turbo(length(levels(metadata$replicate))) %>% 
        magrittr::set_names(levels(metadata$replicate))),
  ## set hierarchical clustering by sample distances
  clustering_distance_rows=dists,
  clustering_distance_cols=dists,
  col=viridis::viridis(100),
  show_rownames=FALSE,
  show_colnames=TRUE,
  ## title
  main = "Pairwise Distances of All Samples",
  cellwidth=10,
  cellheight=10,
  fontsize=6)
## view + assess the plot
g
## save the output to figures directory
### as png
png(
  here(outputs$figures.do, "hmap_genes-all_data.add_for_ppr.png"),
  width=18, height=18, units="cm", res=150)
g
dev.off()
### as pdf
pdf(
  here(outputs$figures.do, "hmap_genes-all_data.add_for_ppr.pdf"),
  width=cm_to_inch(18), height=cm_to_inch(18))
g
dev.off()
## cleanup
rm(g, dists, dist_mat, annot)



#################################################################################
# variance partitioning (pre-batch correction)
#################################################################################

# calculate and plot proportion of variance explained by treatment vs replicate

## run for 'old_add' dataset
### fit variance partitioning model for old_add data
varpart_add = variancePartition::fitExtractVarPartModel(
    vsd_genes_add, 
    ~ treatment + replicate, ## specify model variables
    tibble::column_to_rownames(metadata_add, "sample"))
varpart_add = variancePartition::sortCols(varpart_add)
### plot the results
g = variancePartition::plotVarPart(varpart_add) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Variance explained (sting data)", 
                  caption = "VST normalized counts") +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
### cleanup the plot for visualization
g$layers[[2]] = NULL
### calculate the median proportion of variance explained and add it to the plot at the corresponding position along y
purrr::walk2(
    names(varpart_add), 1:length(varpart_add), function (name, x) {
      m = round(median(varpart_add[[name]])*100, 1)
      g <<- g + ggplot2::annotate("text", x=x, y=m, label=m)
    })
### view + assess the plot
g
### save the output to figures directory
ggsaver(
    "variance_partitioning-genes-sting_data.add_for_ppr", g, path=outputs$figures,
    height=10, width=15, units="cm")

## note: received following warning message when fitting variance partition model on old_add data (i.e. creating varpart_add) (18.3.25)
  # Warning message:
  #   In .fitExtractVarPartModel(exprObj, formula, data, REML = REML,  :
  #                                Model failed for 23730 responses.
  #                              See errors with attr(., 'errors')
                             
## run for 'old' dataset (shows potentially problematic replicate variability)
### fit variance partitioning model
varpart = variancePartition::fitExtractVarPartModel(
    vsd_genes_old, 
    ~ treatment + replicate, ## specify model variables
    tibble::column_to_rownames(metadata_old, "sample"))
varpart = variancePartition::sortCols(varpart)
### plot the results
g = variancePartition::plotVarPart(varpart) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Variance explained (irf3 data)", 
                caption = "VST normalized counts") +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
### cleanup the plot for visualization
g$layers[[2]] = NULL
### calculate the median proportion of variance explained and add it to the plot at the corresponding position along y
purrr::walk2(
    names(varpart), 1:length(varpart), function (name, x) {
      m = round(median(varpart[[name]])*100, 1)
      g <<- g + ggplot2::annotate("text", x=x, y=m, label=m)
    })
### view + assess the plot
g
### save the output to figures directory
ggsaver(
    "variance_partitioning-genes-old_data.add_for_ppr", g, path=outputs$figures,
    height=10, width=15, units="cm")

## cleanup
rm(g, varpart, varpart_add)


#################################################################################
# batch effect correction
#################################################################################

# batch effect correction using RUV-s method
  ## performing batch correction on irf3 data, on old_add data, and on all old data

## run RUV correction on old dataset samples
ruv_corrected =
  RUVSeq::RUVs(
    vsd_genes_old,
    rownames(dds_genes_old),
    1,
    RUVSeq::makeGroups(metadata_old$treatment),
    isLog=TRUE)
## extract unwanted variation
W_1 = ruv_corrected$W
### is the first factor in the unwanted variation, will be used to adjust the data in the following steps
### remove batch effects using linear model approach in limma
#### genes
dds_genes_old$W_1 = W_1
ruv_vsd_genes_old = 
  limma::removeBatchEffect(
    assay(DESeq2::vst(dds_genes_old)),
    design=model.matrix(~treatment, colData(dds_genes_old)), ## specify variation due to treatment to be preserved
    covariates=W_1) # set the covariate to adjust for in the batch removal process
#### te loci
# dds_teloci_old$W_1 = W_1
# ruv_vsd_teloci_old = 
#   limma::removeBatchEffect(
#     assay(DESeq2::vst(dds_teloci_old)),
#     design=model.matrix(~treatment, colData(dds_teloci_old)), ## specify variation due to treatment to be preserved
#     covariates=W_1) # set the covariate to adjust for in the batch removal process
# #### te families
# dds_tefam_old$W_1 = W_1
# ruv_vsd_tefam_old = 
#   limma::removeBatchEffect(
#     assay(DESeq2::vst(dds_tefam_old, nsub=(.70 * nrow(dds_tefam_old)))),  ## had to set nsub due to the large amount of rows with mean normalized count > 5 (pulling an error). set it to 75% of the total number of features in the te family dataset (total no features = 1199, 70% = 839.3). note: the default for nsub is 1000
#     design=model.matrix(~treatment, colData(dds_tefam_old)), ## specify variation due to treatment to be preserved
#     covariates=W_1) # set the covariate to adjust for in the batch removal process
### Note: The vst function uses the estimated dispersion-mean trend to stabilize the variance across the range of mean expressions. If nsub is too small, the estimation might not capture the true relationship, especially if the dataset is complex with different trends in different subsets of features.

## run RUV correction on irf3 dataset samples
ruv_corrected =
  RUVSeq::RUVs(
    vsd_genes_irf,
    rownames(dds_genes_irf),
    1,
    RUVSeq::makeGroups(metadata_irf$treatment),
    isLog=TRUE)
## extract unwanted variation
W_1 = ruv_corrected$W
### is the first factor in the unwanted variation, will be used to adjust the data in the following steps
### remove batch effects using linear model approach in limma
#### genes
dds_genes_irf$W_1 = W_1
ruv_vsd_genes_irf = 
  limma::removeBatchEffect(
    assay(DESeq2::vst(dds_genes_irf)),
    design=model.matrix(~treatment, colData(dds_genes_irf)), ## specify variation due to treatment to be preserved
    covariates=W_1) # set the covariate to adjust for in the batch removal process
#### te loci
# dds_teloci_irf$W_1 = W_1
# ruv_vsd_teloci_irf = 
#   limma::removeBatchEffect(
#     assay(DESeq2::vst(dds_teloci_irf)),
#     design=model.matrix(~treatment, colData(dds_teloci_irf)), ## specify variation due to treatment to be preserved
#     covariates=W_1) # set the covariate to adjust for in the batch removal process
# #### te families
# dds_tefam_irf$W_1 = W_1
# ruv_vsd_tefam_irf = 
#   limma::removeBatchEffect(
#     assay(DESeq2::vst(dds_tefam_irf, nsub=(.70 * nrow(dds_tefam_irf)))),  ## had to set nsub due to the large amount of rows with mean normalized count > 5 (pulling an error). set it to 75% of the total number of features in the te family dataset (total no features = 1199, 70% = 839.3). note: the default for nsub is 1000
#     design=model.matrix(~treatment, colData(dds_tefam_irf)), ## specify variation due to treatment to be preserved
#     covariates=W_1) # set the covariate to adjust for in the batch removal process

## run RUV correction on old_add dataset samples
ruv_corrected =
  RUVSeq::RUVs(
    vsd_genes_add,
    rownames(dds_genes_add),
    1,
    RUVSeq::makeGroups(metadata_add$treatment),
    isLog=TRUE)
## extract unwanted variation
W_1 = ruv_corrected$W
### is the first factor in the unwanted variation, will be used to adjust the data in the following steps
### remove batch effects using linear model approach in limma
#### genes
dds_genes_add$W_1 = W_1
ruv_vsd_genes_add = 
  limma::removeBatchEffect(
    assay(DESeq2::vst(dds_genes_add)),
    design=model.matrix(~treatment, colData(dds_genes_add)), ## specify variation due to treatment to be preserved
    covariates=W_1) # set the covariate to adjust for in the batch removal process
#### te loci
# dds_teloci_add$W_1 = W_1
# ruv_vsd_teloci_add = 
#   limma::removeBatchEffect(
#     assay(DESeq2::vst(dds_teloci_add)),
#     design=model.matrix(~treatment, colData(dds_teloci_add)), ## specify variation due to treatment to be preserved
#     covariates=W_1) # set the covariate to adjust for in the batch removal process
#### te families
# dds_tefam_add$W_1 = W_1
# ruv_vsd_tefam_add = 
#   limma::removeBatchEffect(
#     assay(DESeq2::vst(dds_tefam_add, nsub=(.70 * nrow(dds_tefam_add)))),  ## had to set nsub due to the large amount of rows with mean normalized count > 5 (pulling an error). set it to 75% of the total number of features in the te family dataset (total no features = 1199, 70% = 839.3). note: the default for nsub is 1000
#     design=model.matrix(~treatment, colData(dds_tefam_add)), ## specify variation due to treatment to be preserved
#     covariates=W_1) # set the covariate to adjust for in the batch removal process

## cleanup
rm(ruv_corrected)

#################################################################################
# batch corrected data exploration - plotting PCAs and variance partitioning
#################################################################################

# assessing the results of batch correction on old_add

# plotting PCAs of ruv-corrected old_add samples to check the impact of RUV correction on sample clustering
## perform PCA ('old_add', ruv samples)
pca = PCAtools::pca(
  ruv_vsd_genes_add,
  tibble::column_to_rownames(metadata_add, "sample"), # use filtered metadata
  removeVar=0.1)
# prepare a df with PCs 1 and 2 ('old_add' samples)
pcadf = data.frame(PC1=pca$rotated$PC1, PC2=pca$rotated$PC2, pca$metadata)
# prepare one with PCs 3 and 4 ('old_add' samples)
pcadf2 = data.frame(PC3=pca$rotated$PC3, PC4=pca$rotated$PC4, pca$metadata)
### generate PCA plot of 'old_add' samples, labeled by treatment
#### PC1 + PC2
g = ggplot2::ggplot(pcadf, aes(PC1, PC2, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata_add$sample), show.legend=FALSE, max.overlaps=nrow(metadata_add)) +
  ggplot2::xlab(glue('PC1: {round(pca$variance["PC1"])}% variance')) +
  ggplot2::ylab(glue('PC2: {round(pca$variance["PC2"])}% variance')) +
  ggplot2::labs(title = "PCA of STING dataset (old_add) samples by treatment", caption = "VST norm counts, RUV corrected") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  #ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-sting_data_treatment-ruv_corrected.pc1_2.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
#### PC3 + PC4
g = ggplot2::ggplot(pcadf2, aes(PC3, PC4, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata_add$sample), show.legend=FALSE, max.overlaps=nrow(metadata_add)) +
  ggplot2::xlab(glue('PC3: {round(pca$variance["PC3"])}% variance')) +
  ggplot2::ylab(glue('PC4: {round(pca$variance["PC4"])}% variance')) +
  ggplot2::labs(title = "PCA of STING dataset (old_add) samples by treatment", caption = "VST norm counts, RUV corrected") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  #ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-sting_data_treatment-ruv_corrected.pc3_4.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
### cleanup
rm(pca, pcadf, pcadf2, g)

# plotting PCAs of ruv-corrected 'old' (ALL OLD) samples for (potential) use in the paper
## perform PCA ('old_add', ruv samples)
pca = PCAtools::pca(
  ruv_vsd_genes_old,
  tibble::column_to_rownames(metadata_old, "sample"), # use filtered metadata
  removeVar=0.1)
# prepare a df with PCs 1 and 2 ('old_add' samples)
pcadf = data.frame(PC1=pca$rotated$PC1, PC2=pca$rotated$PC2, pca$metadata)
# prepare one with PCs 3 and 4 ('old_add' samples)
pcadf2 = data.frame(PC3=pca$rotated$PC3, PC4=pca$rotated$PC4, pca$metadata)
### generate PCA plot of 'old_add' samples, labeled by treatment
#### PC1 + PC2
g = ggplot2::ggplot(pcadf, aes(PC1, PC2, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata_add$sample), show.legend=FALSE, max.overlaps=nrow(metadata_add)) +
  ggplot2::xlab(glue('PC1: {round(pca$variance["PC1"])}% variance')) +
  ggplot2::ylab(glue('PC2: {round(pca$variance["PC2"])}% variance')) +
  ggplot2::labs(title = "PCA of original dataset (old) samples by treatment", caption = "VST norm counts, RUV corrected") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  #ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-old_data_treatment-ruv_corrected.pc1_2.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
#### PC3 + PC4
g = ggplot2::ggplot(pcadf2, aes(PC3, PC4, color=treatment)) +
  ggplot2::geom_point(size=2) +
  #ggrepel::geom_text_repel(mapping=aes(label=metadata_add$sample), show.legend=FALSE, max.overlaps=nrow(metadata_add)) +
  ggplot2::xlab(glue('PC3: {round(pca$variance["PC3"])}% variance')) +
  ggplot2::ylab(glue('PC4: {round(pca$variance["PC4"])}% variance')) +
  ggplot2::labs(title = "PCA of original dataset (old) samples by treatment", caption = "VST norm counts, RUV corrected") +
  ggplot2::scale_x_continuous(expand=expansion(mult=0.3)) +
  #ggplot2::scale_color_viridis_d(option="turbo") +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio=1)
#### view + assess the plot
g
#### save the output to figures directory
ggsaver(
  "pca_genes-old_data_treatment-ruv_corrected.pc3_4.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
### cleanup
rm(pca, pcadf, pcadf2, g)

# calculating and plotting the proportion of variance explained in ruv corrected 'old_add'
## run for 'old_add' ruv-corrected dataset
### fit variance partitioning model for old_add data
varpart_add = variancePartition::fitExtractVarPartModel(
  ruv_vsd_genes_add, 
  ~ treatment + replicate, ## specify model variables
  tibble::column_to_rownames(metadata_add, "sample"))
varpart_add = variancePartition::sortCols(varpart_add)
### plot the results
g = variancePartition::plotVarPart(varpart_add) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Variance explained (sting data)", 
                caption = "VST norm counts, RUV corrected") +
  ggplot2::theme(
    axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
### cleanup the plot for visualization
g$layers[[2]] = NULL
### calculate the median proportion of variance explained and add it to the plot at the corresponding position along y
purrr::walk2(
  names(varpart_add), 1:length(varpart_add), function (name, x) {
    m = round(median(varpart_add[[name]])*100, 1)
    g <<- g + ggplot2::annotate("text", x=x, y=m, label=m)
  })
### view + assess the plot
g
### save the output to figures directory
ggsaver(
  "variance_partitioning-genes-sting_data-ruv_corrected.add_for_ppr", g, path=outputs$figures,
  height=10, width=15, units="cm")
### cleanup
rm(g, varpart_add)
 

#################################################################################
# prepare summarized expression data, and import relevant gene lists
#################################################################################

# prepare and export summarized expression data 
## prepare a df of expression data and summarize for each treatment (median exp value)
exprs = as_tibble(round(vsd_genes, 3), rownames="gene_id")
for (grp in levels(metadata$treatment)) {
  exprs[[grp]] = apply(exprs[, metadata[metadata$treatment==grp,]$sample], 1, median)
}
## add gene annotations
exprs = annotate_genes(exprs, org_db)
### export table containing the expression data (including summarized for each treatment)
readr::write_tsv(
  exprs, here(outputs$results_g, "vst-genes-summarized_expression-all_data.add_for_ppr.tsv"))
### cleanup
rm(grp)


#################################################################################
# differential expression analysis
#################################################################################

# de analysis
## creating separate models for new data, irf3 data, and sting (aka old_add) data
## performing the following comparisons, each essentially comparing treatment to its control
   # IRF3-5D vs BFP (irf3 dataset model)
   #   Design: `~ W_1 + treatment` incorporates RUV correction by controlling for the latent factor, and tests the effect of treatment
   # STING GOF vs BFP (old_add dataset model)
   #   Design: `~ W_1 + treatment` incorporates RUV correction by controlling for the latent factor, and tests the effect of treatment
   # STING LOF vs BFP (old_add dataset model)
   #   Design: `~ W_1 + treatment` incorporates RUV correction by controlling for the latent factor, and tests the effect of treatment
   # E14_IFN vs E14_noIFN (new dataset model)
   #   Design: `~ treatment` measuring differential expression due to treatment, comparing treatment, treatment as the variable of interest
   # MEF_IFN vs MEF_noIFN (new dataset model)
   #   Design: `~ treatment` identifying differential expression due to treatment, treatment as the variable of interest

# significance thresholds
  # pvalue 0.01, no logFC threshold

# run DESeq analysis, set contrasts, analyze contrasts, plot MA plots and heatmaps

## run DESeq2
### for irf3 datasets
DESeq2::design(dds_genes_irf) = ~ W_1 + treatment ## specify design to account for RUV
dds_genes_irf = DESeq2::DESeq(dds_genes_irf)
### for old_add datasets
DESeq2::design(dds_genes_add) = ~ W_1 + treatment
dds_genes_add = DESeq2::DESeq(dds_genes_add)
### for new datasets
DESeq2::design(dds_genes_new) = ~ treatment
dds_genes_new = DESeq2::DESeq(dds_genes_new)
## define the DE contrasts to be tested
contrasts_genes = list(  ## initialize a list of contrasts
  ## contrast no 1
  "irf35d-vs-bfp"=list(  ## name the contrast
    contrast_args = list(  ## specify args for this specific contrast
      contrast = c("treatment", "irf35d", "bfp"),  ## the design of the contrast
      object = dds_genes_irf, ## the DESeqDataSet to use for this contrast (relevant as we're running DESeq calculation separately for new and old, see rmd for rationale)
      alpha = stats$de_alpha,  ## padj significance threshold
      lfcThreshold = stats$de_l2FC_old,  ## log2 FC threshold
      altHypothesis = "greaterAbs",  ## alternative hypothesis, set to two-sided test ('greater in absolute value'), so IDs genes as DE regardless of direction (pos or neg), as long as they pass the thresholds
      vst_obj = ruv_vsd_genes_irf
    ),
    contrast_samples=  ## retrieve the samples involved in the contrast by filtering metadata
      metadata %>%
      dplyr::filter(dataset=="old") %>% 
      dplyr::pull(sample)
  ),
  ## contrast no 2
  "sting_gof-vs-bfp"=list(  ## name the contrast
    contrast_args = list(  ## specify args for this specific contrast
      contrast = c("treatment", "sting_gof", "bfp"),  ## the design of the contrast
      object = dds_genes_add, ## the DESeqDataSet to use for this contrast (relevant as we're running DESeq calculation separately for new and old, see rmd for rationale)
      alpha = stats$de_alpha,  ## padj significance threshold
      lfcThreshold = stats$de_l2FC_old,  ## log2 FC threshold
      altHypothesis = "greaterAbs",  ## alternative hypothesis, set to two-sided test ('greater in absolute value'), so IDs genes as DE regardless of direction (pos or neg), as long as they pass the thresholds
      vst_obj = ruv_vsd_genes_add
    ),
    contrast_samples=  ## retrieve the samples involved in the contrast by filtering metadata
      metadata %>%
      dplyr::filter(treatment %in% c("sting_gof", "bfp")) %>% 
      dplyr::pull(sample)
  ),
  ## contrast no 3
  "sting_lof-vs-bfp"=list(  ## name the contrast
    contrast_args = list(  ## specify args for this specific contrast
      contrast = c("treatment", "sting_lof", "bfp"),  ## the design of the contrast
      object = dds_genes_add, ## the DESeqDataSet to use for this contrast (relevant as we're running DESeq calculation separately for new and old, see rmd for rationale)
      alpha = stats$de_alpha,  ## padj significance threshold
      lfcThreshold = stats$de_l2FC_old,  ## log2 FC threshold
      altHypothesis = "greaterAbs",  ## alternative hypothesis, set to two-sided test ('greater in absolute value'), so IDs genes as DE regardless of direction (pos or neg), as long as they pass the thresholds
      vst_obj = ruv_vsd_genes_add
    ),
    contrast_samples=  ## retrieve the samples involved in the contrast by filtering metadata
      metadata %>%
      dplyr::filter(treatment %in% c("sting_lof", "bfp")) %>% 
      dplyr::pull(sample)
  ),
  ## contrast no 4
  "E14_IFN-vs-E14_noIFN"=list(
    contrast_args = list(
      contrast = c("treatment", "mESC_IFN", "mESC_noIFN"),
      object = dds_genes_new, ## the DESeqDataSet to use for this contrast (relevant as we're running DESeq calculation separately for new and old, see rmd for rationale)
      alpha = stats$de_alpha,
      lfcThreshold = stats$de_l2FC_E14,
      altHypothesis = "greaterAbs",
      vst_obj = vsd_genes_new
    ),
    contrast_samples=  ## retrieve the samples involved in the contrast by filtering metadata
      metadata %>%
      dplyr::filter(treatment %in% c("mESC_IFN", "mESC_noIFN")) %>% 
      dplyr::pull(sample)
  ),
  ## contrast no 5
  "MEF_IFN-vs-MEF_noIFN"=list(  ## name the contrast
    contrast_args = list(
      contrast = c("treatment", "MEF_IFN", "MEF_noIFN"),
      object = dds_genes_new, ## the DESeqDataSet to use for this contrast (relevant as we're running DESeq calculation separately for new and old, see rmd for rationale)
      alpha = stats$de_alpha,
      lfcThreshold = stats$de_l2FC_MEF,
      altHypothesis = "greaterAbs",
      vst_obj = vsd_genes_new
    ),
    contrast_samples=  ## retrieve the samples involved in the contrast by filtering metadata
      metadata %>%
      dplyr::filter(treatment %in% c("MEF_IFN", "MEF_noIFN")) %>% 
      dplyr::pull(sample)
  )
)
## report the thresholds for DE analysis
cat("For differential expression analysis, significance threshold set to padj", stats$de_alpha, "and log2 FC threshold set to", stats$de_l2FC_old, "for irf3 (old) and old_add datasets,", stats$de_l2FC_MEF, "for MEFs and", stats$de_l2FC_E14, "for E14 (new dataset)", "\n")

## run the contrasts using custom function (contrast_and_save)
### as output, for each contrast it produces: a list containing results for each contrast, an MA plot output to figures directory, a table containing results of the contrast to results directory
purrr::map(
  names(contrasts_genes),
  function(contrast) {
    print(contrasts_genes[[contrast]]$contrast_args)
    contrasts_genes[[contrast]]$results <<- invisible(rlang::exec(
      contrast_and_save,
      contrast_name=contrast, feature_type="genes",
      #object=dds_genes,  ## NOW specifying within contrast_args, enabling different DDS objects for each contrast
      filter_result=T, add_names=T, write_table=T,
      output_location_tab=here(outputs$res_contrasts_g), ## specify output location for the DE results tables
      output_location_fig=here(outputs$figures.do), ## specify output location for the MA plots
      !!!contrasts_genes[[contrast]]$contrast_args)
    )})

## plot heatmaps of DEGs for each contrast
purrr::map(
  names(contrasts_genes),
  function(contrast) {
    if (nrow(contrasts_genes[[contrast]]$results) == 0) return(NULL)
    exprs = contrasts_genes[[contrast]]$contrast_args$vst_obj[
      rownames(contrasts_genes[[contrast]]$contrast_args$object) %in% contrasts_genes[[contrast]]$results$gene_id,
      colnames(contrasts_genes[[contrast]]$contrast_args$object) %in% contrasts_genes[[contrast]]$contrast_samples]
    if (is.null(nrow(exprs))) return(NULL)
    exprs = exprs[apply(exprs, 1, function(x) length(unique(x)) != 1), ]
    if (is.null(nrow(exprs))) return(NULL)
    annot = metadata %>%
      tibble::column_to_rownames("sample") %>%
      dplyr::select(treatment)
    ## define the heatmap plot
    g <- pheatmap::pheatmap(
      exprs,
      annotation = annot,
      annotation_colors = list(
        "treatment" = viridis::turbo(length(levels(metadata$treatment))) %>% 
          magrittr::set_names(levels(metadata$treatment))
      ),
      scale = "row",
      cluster_rows = TRUE, ## hierarchical clustering of genes
      cluster_cols = TRUE, ## hierarchical clustering of samples
      show_rownames = FALSE,
      show_colnames = TRUE,
      cellwidth = 10,
      color = viridis::viridis(100),
      border_color = NA,
      fontsize = 8
    )
    ## view the hmaps
    g
    ## export hmaps to figures directory
    pdf(
      here(outputs$figures.do, glue("de_hmap_genes-{contrast}_heatmap.add_for_ppr.pdf")),
      width = cm_to_inch(15), height = cm_to_inch(20)
    )
    print(g)
    dev.off()
    png(
      here(outputs$figures.do, glue("de_hmap_genes-{contrast}_heatmap.add_for_ppr.png")),
      width = 15, height = 20, units = "cm", res = 150
    )
    print(g)
    dev.off()
  }
)
## report the number of significant DEGs in each contrast from this analysis
# dim_irf3 <- nrow(contrasts_genes$`irf35d-vs-bfp`$results)
# dim_mef <- nrow(contrasts_genes$`MEF_IFN-vs-MEF_noIFN`$results)
# dim_e14 <- nrow(contrasts_genes$`E14_IFN-vs-E14_noIFN`$results)
# ### print the output
# cat(dim_irf3, "significant DEGs were identified in irf35d-vs-bfp;", dim_mef, "significant DEGs were identified in MEF_IFN-vs-MEF_noIFN;", dim_e14, "significant DEGs were identified in E14_IFN-vs-E14_noIFN", "\n")

## cleanup
rm(dim_irf3, dim_mef, dim_e14)


# plot deg counts
  # plot the numbers and directionalities of significant DEGs identified in each contrast
## initialize empty df
df = data.frame(contrast = character(), direction = character(), degs = integer(), stringsAsFactors = FALSE)
## fill df with info from each contrast
for (contrast in names(contrasts_genes)) {
  ## specify contrast results variable
  i = contrasts_genes[[contrast]]$results
  ## specify log2fc for each contrast
  lfc = contrasts_genes[[contrast]]$contrast_args$lfcThreshold
  ## extract contrast, direction, and no. of degs
  rowUp = data.frame(contrast = contrast, direction = "up", degs = nrow(i[i$log2FoldChange > lfc ,]))
  rowDown = data.frame(contrast = contrast, direction = "down", degs = nrow(i[i$log2FoldChange < lfc ,]))
  ## combine the results into df
  df = rbind(df, rowUp, rowDown)
}
## process df output: changing order of levels - to have upregulated on top
df$direction <- factor(df$direction, levels=c("up", "down"))
### check the output
print(df)
## plot a barplot of the number of DEGs in each comparison
g <- ggplot(df) +
  ggplot2::geom_bar(aes(
    x = reorder(contrast, -degs), ## ordering bars to plot in decreasing order
    y = degs, fill = direction),
    stat = "identity") +
  ## add number labels over the bars to display DEG counts
  ggplot2::geom_text(aes(
    x = reorder(contrast, -degs), 
    y = degs + 
      0.08 * max(degs), ## add offset so numbers are visible for small bars
    label = degs,
    group = direction
  ),
  position = position_stack(vjust = 0.5)
  ) +
  ggplot2::xlab(NULL) +
  ggplot2::ylab("DEG Count") +
  ggplot2::labs(title = "Differentially expressed genes per contrast", caption = glue('padj {stats$de_alpha}, logFC for: irf3+sting={stats$de_l2FC_old}, MEFs={stats$de_l2FC_MEF}, E14={stats$de_l2FC_E14}')) +
  ggplot2::theme_classic()+
  ggplot2::scale_fill_manual(values = c("#C44550", "#4D66A0")) +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust=1),
                 legend.position = c(0.95, 0.95), legend.justification = c("right", "top"))
### view the plot
g
### save the output to figures directory
ggsaver(
  "de_numbers_bar_genes-allcontrasts.add_for_ppr", g, path=outputs$figures,
  height=20, width=14, units="cm")
## cleanup
rm(df, g, rowUp, rowDown)




#################################################################################
# gene ontology analysis
#################################################################################

# gene ontology analysis using clusterprofiler, performing separately for each directionality (GO:BP)

## run go analysis of contrasts results
  ## running GO analysis using clusterProfiler, accounting for directionality, extracting BP terms
for (contrast in names(contrasts_genes)) {
  ## for UP genes
  cat("Analysis: goUP of contrast", contrast, "with qvalue threshold of", stats$enrich_alpha, "\n")
  ## ensure there are significant genes present in the relevant direction, then run GO
  geneUP = dplyr::filter(contrasts_genes[[contrast]]$results, log2FoldChange>0) %>% dplyr::pull(gene_id)
  if (length(geneUP) < stats$go_minimum) {  ## setting a requirement for ... or more genes to be passed onto GO analysis
    cat("Less than", stats$go_minimum, "significant UP genes found in contrast", contrast, "\n")
  } else {
    goUP = clusterProfiler::enrichGO(
      gene=geneUP,
      OrgDb=org.Mm.eg.db, keyType='ENSEMBL', ont="BP",
      pAdjustMethod="BH", qvalueCutoff=stats$enrich_alpha, readable=T)
    ## ensure there are significant results before plotting/exporting
    if (is.null(goUP) || nrow(goUP@result) == 0) {
      cat("No significant GO terms found for UP genes in contrast", contrast, 
          "at qval", stats$enrich_alpha, "\n")
    } else {
      ### export table of GO results
      readr::write_tsv(
        dplyr::filter(goUP@result, qvalue<stats$enrich_alpha),
        here(outputs$res_go_g, glue("de-go_genes-{contrast}_goUP.add_for_ppr.txt")))
      ### plot and export: previous way
      n_enriched = nrow(dplyr::filter(goUP@result, qvalue<stats$enrich_alpha))
      if (n_enriched) {
        g = clusterProfiler::dotplot(goUP, showCategory=25) +
          ggplot2::labs(title = glue("GO:BP of UP in {contrast}"), 
                        caption = glue('qvalue {stats$enrich_alpha}')) +
          ggplot2::theme(axis.text.y=element_text(size=8))
        #### visualize the plot
        plot(g)
        #### export
        ggsaver(
          glue("de-go_dotplot-genes-{contrast}_goUP.add_for_ppr"), g, path=outputs$figures,
          height=7.5+(0.8*(min(n_enriched, 25))), width=15, units="cm")
      ### plot and export: new way
      n_enriched = nrow(dplyr::filter(goUP@result, qvalue<stats$enrich_alpha))
      if (n_enriched) {
        g = clusterProfiler::dotplot(goUP, showCategory=15) +
          ggplot2::labs(title = glue("GO:BP of UP in {contrast}"), 
                        caption = glue('qvalue {stats$enrich_alpha}')) +
          ggplot2::theme(axis.text.y=element_text(size=11),
                         legend.text = element_text(size = 11),
                         legend.title = element_text(size = 13),
                         legend.key.size = unit(1, "cm"))
        #### visualize the plot
        plot(g)
        #### export
        ggsaver(
          glue("de-go_dotplot-genes-{contrast}_goUP.v2.add_for_ppr"), g, path=outputs$figures,
          height=7.5+(0.8*(min(n_enriched, 15))), width=15, units="cm")
      }
    }
    }
  }
  ## for DOWN genes
  cat("Analysis: goDOWN of contrast", contrast, "with qvalue threshold of", stats$enrich_alpha, "\n")
  ## ensure there are significant genes present in the relevant direction, then run GO
  geneDOWN = dplyr::filter(contrasts_genes[[contrast]]$results, log2FoldChange<0) %>% dplyr::pull(gene_id)
  if (length(geneDOWN) < stats$go_minimum) {  ## setting a requirement for ... or more genes to be passed onto GO analysis
    cat("Less than", stats$go_minimum, "significant DOWN genes found in contrast", contrast, "\n")
  } else {
    goDOWN = clusterProfiler::enrichGO(
      gene=geneDOWN,
      OrgDb=org.Mm.eg.db, keyType='ENSEMBL', ont="BP",
      pAdjustMethod="BH", qvalueCutoff=stats$enrich_alpha, readable=T)
    ## ensure there are significant results before plotting/exporting
    if (is.null(goDOWN) || nrow(goDOWN@result) == 0) {
      cat("No significant GO terms found for DOWN genes in contrast", contrast, "at qval",
          stats$enrich_alpha, "\n")
    } else {
      ### export table of GO results
      readr::write_tsv(
        dplyr::filter(goDOWN@result, qvalue<stats$enrich_alpha),
        here(outputs$res_go_g, glue("de-go_genes-{contrast}_goDOWN.add_for_ppr.txt")))
      ### plot and export: previous way
      n_enriched = nrow(dplyr::filter(goDOWN@result, qvalue<stats$enrich_alpha))
      if (n_enriched) {
        g = clusterProfiler::dotplot(goDOWN, showCategory=25) +
          ggplot2::labs(title = glue("GO:BP of DOWN in {contrast}"), 
                        caption = glue('qvalue {stats$enrich_alpha}')) +
          ggplot2::theme(axis.text.y=element_text(size=8))
        ### visualize the plot
        plot(g)
        ### export
        ggsaver(
          glue("de-go_dotplot-genes-{contrast}_goDOWN.add_for_ppr"), g, path=outputs$figures,
          height=7.5+(0.8*(min(n_enriched, 25))), width=15, units="cm")
      ### plot and export: new way
      n_enriched = nrow(dplyr::filter(goDOWN@result, qvalue<stats$enrich_alpha))
      if (n_enriched) {
        g = clusterProfiler::dotplot(goDOWN, showCategory=15) +
          ggplot2::labs(title = glue("GO:BP of DOWN in {contrast}"), 
                        caption = glue('qvalue {stats$enrich_alpha}')) +
          ggplot2::theme(axis.text.y=element_text(size=11),
                         legend.text = element_text(size = 11),
                         legend.title = element_text(size = 13),
                         legend.key.size = unit(1, "cm"))
        ### visualize the plot
        plot(g)
        ### export
        ggsaver(
          glue("de-go_dotplot-genes-{contrast}_goDOWN.v2.add_for_ppr"), g, path=outputs$figures,
          height=7.5+(0.8*(min(n_enriched, 15))), width=15, units="cm")
        }
      }
    }
  }
}
## cleanup
rm(contrast, goUP, goDOWN, geneUP, geneDOWN, g, n_enriched)


#################################################################################
# overlap analysis of upregulated genes
#################################################################################

## load in vectors of upregulated genes in desired comparisons
### initialize an empty list
venn_g <- list()
### extract gene IDs into the list of vectors
for (contrast in names(contrasts_genes)) {
  ## specify contrast results variable
  i = contrasts_genes[[contrast]]$results
  ## extract upregulated gene_ids (for all contrasts, I'm interested in the upregulated genes)
  gen = dplyr::filter(i, log2FoldChange>0) %>% dplyr::pull(gene_id)
    ## ensure it's a character vector
    gen = as.character(gen)
  ## store the gene IDs in the list of vectors, labeling with direction
  venn_g[[glue('{contrast}--UP')]] <- gen
}

# plotting Venn of all 3 contrasts: E14_IFN vs MEF_IFN vs IRF35D
## identify unique and conserved elements across all 3 contrasts
  
## plot and export the venn diagram
g = plot(euler(venn_g), 
     #fill = three_col,
     main = "Overlapping DEGs: upreg in E14_IFN, MEF_IFN, and IRF35D",
     quantities = TRUE)
### view the plot
g
### export the plot to figures directory
pdf(
  here(outputs$figures.do, "de_venn_genes-E14_IFN-vs-MEF_IFN-vs-irf35d--up.pdf"),
  width = cm_to_inch(20), height = cm_to_inch(15)
)
print(g)
dev.off()
png(
  here(outputs$figures.do, "de_venn_genes-E14_IFN-vs-MEF_IFN-vs-irf35d--up.png"),
  width = 20, height = 15, units = "cm", res = 150
)
print(g)
dev.off()

## extract and save conserved and unique genes from the venn comparison  (all 3 treatments)
### retrieve the gene_ids that are unique and conserved in E14 vs MEFs vs IRF35D
conserved_genes <- Reduce(intersect, venn_g) ## conserved across all
u_E14 <- setdiff(venn_g[["E14_IFN-vs-E14_noIFN--UP"]], 
                             c(venn_g[["MEF_IFN-vs-MEF_noIFN--UP"]], 
                               venn_g[["irf35d-vs-bfp--UP"]]))  ## unique in E14
u_MEF <- setdiff(venn_g[["MEF_IFN-vs-MEF_noIFN--UP"]], 
                             c(venn_g[["E14_IFN-vs-E14_noIFN--UP"]], 
                               venn_g[["irf35d-vs-bfp--UP"]]))  ## unique in MEF
u_irf3 <- setdiff(venn_g[["irf35d-vs-bfp--UP"]], 
                             c(venn_g[["E14_IFN-vs-E14_noIFN--UP"]], 
                               venn_g[["MEF_IFN-vs-MEF_noIFN--UP"]]))  ## unique in IRF35D
sh_MEF_irf3 <- setdiff(
  intersect(venn_g[["irf35d-vs-bfp--UP"]], venn_g[["MEF_IFN-vs-MEF_noIFN--UP"]]), 
  venn_g[["E14_IFN-vs-E14_noIFN--UP"]])  ## shared in MEFs and IRF35d, but not in E14
sh_E14_irf3 <- setdiff(
  intersect(venn_g[["irf35d-vs-bfp--UP"]], venn_g[["E14_IFN-vs-E14_noIFN--UP"]]), 
  venn_g[["MEF_IFN-vs-MEF_noIFN--UP"]])  ## shared in E14 and IRF35d, but not in MEFs
sh_E14_MEF <- setdiff(
  intersect(venn_g[["MEF_IFN-vs-MEF_noIFN--UP"]], venn_g[["E14_IFN-vs-E14_noIFN--UP"]]), 
  venn_g[["irf35d-vs-bfp--UP"]])  ## shared in E14 and MEFs, but not in IRF35D
### create a data frame for each context to merge
  ### initialize an empty list to hold the data frames
  df_list <- list()
  ### create a data frame for each context if not empty
    ### conserved genes
    if (length(conserved_genes) > 0) {
      df_list$conserved <- data.frame(
        context..all_trt_up = "conserved",
        gene_id = conserved_genes
      )
    }
    ### unique to E14_IFN
    if (length(u_E14) > 0) {
      df_list$unique_E14_IFN <- data.frame(
        context..all_trt_up = "unique-E14_IFN-up",
        gene_id = u_E14
      )
    }
    ### unique to MEF_IFN
    if (length(u_MEF) > 0) {
      df_list$unique_MEF_IFN <- data.frame(
        context..all_trt_up = "unique-MEF_IFN-up",
        gene_id = u_MEF
      )
    }
    ### unique to irf35d
    if (length(u_irf3) > 0) {
      df_list$unique_irf35d <- data.frame(
        context..all_trt_up = "unique-irf35d-up",
        gene_id = u_irf3
      )
    }
    ### shared by MEF and irf35d
    if (length(sh_MEF_irf3) > 0) {
      df_list$shared_MEF_irf3 <- data.frame(
        context..all_trt_up = "shared-MEF_IFN-irf35d-up",
        gene_id = sh_MEF_irf3
      )
    }
    ### shared by E14 and irf35d
    if (length(sh_E14_irf3) > 0) {
      df_list$shared_E14_irf3 <- data.frame(
        context..all_trt_up = "shared-E14_IFN-irf35d-up",
        gene_id = sh_E14_irf3
      )
    }
    ### shared by E14 and MEF
    if (length(sh_E14_MEF) > 0) {
      df_list$shared_E14_MEF <- data.frame(
        context..all_trt_up = "shared-E14_IFN-MEF_IFN-up",
        gene_id = sh_E14_MEF
      )
    }
### combine all non-empty data frames into one
gID_olap_allc_df <- do.call(rbind, df_list)
  ### remove rownames
  rownames(gID_olap_allc_df) <- NULL
### add gene annotations
gID_olap_allc_df = annotate_genes(gID_olap_allc_df, org_db)
    ### check output
    # head(gID_olap_allc_df)
    # dim(gID_olap_allc_df)
    unique(gID_olap_allc_df$context..all_trt_up)
### add significance level information for each sample
#### `irf35d-vs-bfp`
gID_olap_allc_df <- gID_olap_allc_df %>%
  dplyr::left_join(contrasts_genes$`irf35d-vs-bfp`$results %>% dplyr::select(gene_id, log2FoldChange, lfcSE, padj), by = c("gene_id" = "gene_id")) %>%
  rename(`log2fc_irf35d-vs-bfp` = log2FoldChange,
         `lfcSE_irf35d-vs-bfp` = lfcSE,
         `padj_irf35d-vs-bfp` = padj)
    #### double checked that it performed the operation correctly - yes (20.1.25)
#### `E14_IFN-vs-E14_noIFN`
gID_olap_allc_df <- gID_olap_allc_df %>%
  dplyr::left_join(contrasts_genes$`E14_IFN-vs-E14_noIFN`$results %>% dplyr::select(gene_id, log2FoldChange, lfcSE, padj), by = c("gene_id" = "gene_id")) %>%
  rename(`log2fc_E14_IFN-vs-E14_noIFN` = log2FoldChange,
         `lfcSE_E14_IFN-vs-E14_noIFN` = lfcSE,
         `padj_E14_IFN-vs-E14_noIFN` = padj)
#### `MEF_IFN-vs-MEF_noIFN`
gID_olap_allc_df <- gID_olap_allc_df %>%
  dplyr::left_join(contrasts_genes$`MEF_IFN-vs-MEF_noIFN`$results %>% dplyr::select(gene_id, log2FoldChange, lfcSE, padj), by = c("gene_id" = "gene_id")) %>%
  rename(`log2fc_MEF_IFN-vs-MEF_noIFN` = log2FoldChange,
         `lfcSE_MEF_IFN-vs-MEF_noIFN` = lfcSE,
         `padj_MEF_IFN-vs-MEF_noIFN` = padj)
### export the df to relevant results directory
readr::write_tsv(
  gID_olap_allc_df, here(outputs$results_g,
                         glue("de_genes-venn_olaps-E14_IFN-vs-MEF_IFN-vs-irf35d--up.txt")))


#################################################################################
# plotting expression and fold change for each gene set (C/IMI/MI/IR)
#################################################################################

# venn info to extract genes
head(gID_olap_allc_df)
unique(gID_olap_allc_df$context..all_trt_up)
## extract vector of genes
imi <- gID_olap_allc_df$gene_id[gID_olap_allc_df$context..all_trt_up=="shared-MEF_IFN-irf35d-up"]
mi <- gID_olap_allc_df$gene_id[gID_olap_allc_df$context..all_trt_up=="unique-MEF_IFN-up"]
ir <- gID_olap_allc_df$gene_id[gID_olap_allc_df$context..all_trt_up=="unique-irf35d-up"]
c <- gID_olap_allc_df$gene_id[gID_olap_allc_df$context..all_trt_up=="conserved"]

# load normalized counts, including summed for each condition
head(exprs)
colnames(exprs)

# filter exprs df columns
## filter exprs columns to remove summarized, and only keep sample specific
# lim_exprs <- exprs %>% 
#   dplyr::select(-gene_symbol, -gene_name, -bfp, -irf35d, -MEF_IFN, -MEF_noIFN, -mESC_IFN, -mESC_noIFN)
### doing for each group now for sig testing
lim_exprs <- exprs %>% dplyr::select(gene_id, bfp, irf35d, MEF_IFN, MEF_noIFN, mESC_IFN, mESC_noIFN)

# filter exprs df to only include the genes in each venn category (use gene_id, not symbols), and add a column for category
imi_g <- lim_exprs %>% dplyr::filter(gene_id %in% imi)
  imi_g$gene_set <- "imi"
mi_g <- lim_exprs %>% dplyr::filter(gene_id %in% mi)
  mi_g$gene_set <- "mi"
ir_g <- lim_exprs %>% dplyr::filter(gene_id %in% ir)
  ir_g$gene_set <- "ir"
c_g <- lim_exprs %>% dplyr::filter(gene_id %in% c)
  c_g$gene_set <- "c"

# combine each gene set category into one df
com_g <- dplyr::bind_rows(imi_g, mi_g, ir_g, c_g)
## check output
dim(com_g)
head(com_g)
unique(com_g$gene_set)

# prep data for plotting
### reshape the data to long format
#### convert to df (just for vis purposes)
com_g <- as.data.frame(com_g)
#### remove gene_id column (not needed for plotting)
com_g_f <- com_g %>% 
  dplyr::select(-gene_id)
#### reshape the data to long format
com_g_long <- com_g_f %>%
  tidyr::pivot_longer(!gene_set, names_to = "condition", values_to = "norm_exp")
#### check output
head(com_g_long)
unique(com_g_long$condition)
### set levels of conditions for plotting
com_g_long <- com_g_long %>%
  mutate(condition = factor(condition, levels = c("bfp", "irf35d", "mESC_noIFN", "mESC_IFN", "MEF_noIFN", "MEF_IFN")))

# plot a boxplot of expression levels in each condition, do separately for each gene set
## to check baseline expression in different conditions
# plot barplot of expression level
ge1 <- ggplot(com_g_long, aes(x = gene_set, y = norm_exp, fill = condition)) +
  geom_boxplot() +
  labs(
    title = "Expression of each overlap gene set per condition",
    x = "Gene Set",
    y = "Normalized Expression (vst)",
    fill = "Condition"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11), # Rotate x-axis labels
    text = element_text(size = 12)
  ) +
  scale_fill_brewer(palette = "Set2") # Optional: Change the color palette
ge1
  
## exporting expression barplots
ggsaver(
      "de-barplot_venn_sets-norm_exp-all_groups", ge1, path=outputs$figures, 
      height=cm_to_inch(12), width=cm_to_inch(18), units="in",
      device=c("png", "pdf"))


# statistical comparison

# filter exprs df columns
## filter exprs columns to remove summarized, and only keep summarized values for stats
mean_exprs <- exprs %>% 
  dplyr::select(gene_id, bfp, irf35d, MEF_IFN, MEF_noIFN, mESC_IFN, mESC_noIFN)
  
## converting all to long format for testing
imi_g_long <- imi_g %>%
  pivot_longer(cols = -c(gene_id, gene_set),   # Keep gene_id and gene_set as identifiers
               names_to = "condition",         # Convert column names into a 'condition' column
               values_to = "norm_exp")       # Store expression values in 'expression'
mi_g_long <- mi_g %>%
  pivot_longer(cols = -c(gene_id, gene_set),   # Keep gene_id and gene_set as identifiers
               names_to = "condition",         # Convert column names into a 'condition' column
               values_to = "norm_exp")       # Store expression values in 'expression'
ir_g_long <- ir_g %>%
  pivot_longer(cols = -c(gene_id, gene_set),   # Keep gene_id and gene_set as identifiers
               names_to = "condition",         # Convert column names into a 'condition' column
               values_to = "norm_exp")       # Store expression values in 'expression'
c_g_long <- c_g %>%
  pivot_longer(cols = -c(gene_id, gene_set),   # Keep gene_id and gene_set as identifiers
               names_to = "condition",         # Convert column names into a 'condition' column
               values_to = "norm_exp")       # Store expression values in 'expression'


# running stats testing - Wilcoxon signed-rank test, unpaired
full_g_long <- bind_rows(imi_g_long, mi_g_long, ir_g_long, c_g_long)
## run within each gene set
exp_sig_venn <- full_g_long %>%
  group_by(gene_set) %>%
  pairwise_wilcox_test(norm_exp ~ condition, 
                       paired = FALSE,
                       p.adjust.method = "BH")
## evaluate results
View(exp_sig_venn)

## export results table as tsv
readr::write_tsv(
  exp_sig_venn, here(outputs$results_g, glue("de_genes-stats_venn_sets-norm_exp-main_groups.tsv")))


#################################################################################
# plotting gene sets of interests - heatmaps, volcanoes, etc.
#################################################################################

# plotting the following
  # volcano plots of pluripotency genes
  # gsea plot of hecker IFN genes dataset (see S3A)
  # heatmap of apoptosis genes (see sent by MP) (do for all comparisons)
  #   plot with the original gene set sent by MP, and also potentially with another gene set that I find
  #   plotting two versions: one with all samples across all treatments, and one with only old datasets (old_add and irf3) to match existing version in supplementary figure drafts
  # MA plot (or volcano) with other genes highlighted, see S3B

# import/create relevant gene sets
## hecker IFN targets (sent by MP)
hecker <- read.table("./analysis/3.rna_downstream.add_for_ppr/mm39/diff_exp/data/gene_sets/HECKER_IFNB1_TARGETS.tsv", 
                     header = TRUE, sep = "\t")
  ### symbols are all capitalized, so convert them to appropriate names (convert all letters after the first to lowercase)
  hecker$SYMBOL <- str_to_title(str_to_lower(hecker$SYMBOL))
## apoptosis genes (sent by MP, only contains a list of genes)
apop <- readLines("./analysis/3.rna_downstream.add_for_ppr/mm39/diff_exp/data/gene_sets/KEGG_APOPTOSIS.v2023.1.Hs.tsv")
  ### split gene names for vectorization
  apop <- unlist(strsplit(apop, ","))
  ### they are all capitalized, so convert them to appropriate names (convert all letters after the first to lowercase)
  apop <- str_to_title(str_to_lower(apop))
## pluripotency genes (manually adding here)
plg <- c(
  "Nanog",
  "Myc",
  "Klf2",
  "Klf4",
  "Klf5",
  "Nr0b1",
  "Nr5a2",
  "Prdm14",
  "Tbx3",
  "Zfp42",
  "Dppa4",
  "Dppa2",
  "Tfcp2l1",
  "Tfap2c",
  "Tead4",
  "Tcl1",
  "Dppa3",
  "Esrrb",
  "Lifr"
)
## pluripotency genes with some IFN-up genes added (see slack for records)
plg_ifn <- c(
  "Nanog",
  "Myc",
  "Klf2",
  "Klf4",
  "Klf5",
  "Nr0b1",
  "Nr5a2",
  "Prdm14",
  "Tbx3",
  "Zfp42",
  "Dppa4",
  "Dppa2",
  "Tfcp2l1",
  "Tfap2c",
  "Tead4",
  "Tcl1",
  "Dppa3",
  "Esrrb",
  "Lifr",
  "Ifnb1",
  "Irf3",
  "Irf7",
  "Ifna4",
  "Nfkb2",
  "Rela"
)
### and for plotting points with different colors for plg_ifn, loading in the ifn separately too
p_ifn <- c(
  "Ifnb1",
  "Irf3",
  "Irf7",
  "Ifna4",
  "Nfkb2",
  "Rela"
)
## relevant genes of interest
irgen = c("ENSMUSG00000070904", "ENSMUSG00000078354", "ENSMUSG00000095270", "ENSMUSG00000095498",
          "ENSMUSG00000096682", "ENSMUSG00000100079", "ENSMUSG00000048806", "ENSMUSG00000026104",
          "ENSMUSG00000040033", "ENSMUSG00000002325", "ENSMUSG00000018899", "ENSMUSG00000003184",
          "ENSMUSG00000025498", "ENSMUSG00000002325", "ENSMUSG00000031627", "ENSMUSG00000024927",
          "ENSMUSG00000004040", "ENSMUSG00000034394", "ENSMUSG00000012396", "ENSMUSG00000021255",
          "ENSMUSG00000074637", "ENSMUSG00000024406", "ENSMUSG00000055148", "ENSMUSG00000003032",
          "ENSMUSG00000005148", "ENSMUSG00000020167", "ENSMUSG00000051176", "ENSMUSG00000018604",
          "ENSMUSG00000028163", "ENSMUSG00000025225", "ENSMUSG00000002983")


# retrieve results from contrast
## initialize empty list to store results
pl_res_o <- list()
## initialize empty list to store results dfs
pl_res <- list()
## (re-)run the analysis
for (i in names(contrasts_genes)) {
  ## report name of analysis
  print(i)
  ## retrieve results and process
  res_ma <- DESeq2::results(contrasts_genes[[i]]$contrast_args$object, 
                            contrast=contrasts_genes[[i]]$contrast_args$contrast, 
                            lfcThreshold = stats$de_l2FC_old, alpha = stats$de_alpha, 
                            altHypothesis="greaterAbs", filterFun=IHW::ihw)
  resLFC <- DESeq2::lfcShrink(contrasts_genes[[i]]$contrast_args$object, res = res_ma, type="ashr") ## logFC shrinkage via ashr
  ## convert deseq results to df and add gene_id and annotations
  resLFC_df <- as.data.frame(resLFC)
  resLFC_df$gene_id <- rownames(resLFC_df)
    ### add gene annotations
    resLFC_df = annotate_genes(resLFC_df, org_db)
  ## save the results to the list
  pl_res_o[[i]] <- resLFC
  pl_res[[i]] <- resLFC_df
}
### check the output
names(pl_res)
names(pl_res_o)
head(pl_res[[1]])


# plot GSEA of Hecker IFNb targets in IRF3 vs BFP
  ## plot GSEA running score and preranked list of GSEA result
  ## see https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
## specify and process the deseq results table with SIGNIFICANT results only for IRF3vsBFP
### specify the table
ir_res <- contrasts_genes[["irf35d-vs-bfp"]]$results
## create a ranked list, rank genes by log2 fold change (descending order)
ir_rank <- ir_res %>%
  select(gene_id, log2FoldChange) %>%  ## using gene_id to avoid duplicates
  arrange(desc(log2FoldChange)) %>%
  deframe()  ## convert to named vector for fgsea
## identify and retrieve gene_ids for the Hecker IFN targets to avoid duplicate gene_symbols during gsea analysis
### check keytypes
keytypes(org_db)
### extract ensembl gene IDs from gene symbols
heck_ens <- mapIds(org_db,
                   #multiVals = "first",  ## if multiple mappings exist, take the first (not using, not needed here (20.3.25))
                   keys = hecker$SYMBOL,
                   column = "ENSEMBL",
                   keytype = "SYMBOL")
  ### check output
  length(heck_ens)
  length(unique(heck_ens))
  length(hecker$SYMBOL)
## set up df for input into TERM2GENE parameter (requires 2 columns, 1 for gene id or symbol, 1 for gene set name)
gene_set <- data.frame(
  TERM = rep("heck_ens", length(heck_ens)),
  GENE = heck_ens
)
## run GSEA on the ranked list of DEGs from IRF3vBFP
gsea_results <- GSEA(geneList = ir_rank, 
                     TERM2GENE = gene_set, 
                     minGSSize = 10, 
                     maxGSSize = 5000, 
                     pvalueCutoff = 0.05)
  ### check the gsea results and record
  gsea_results
## plot the gsea enrichment plot for the hecker ifn targets
g <- gseaplot2(gsea_results, 
               geneSetID = "heck_ens", 
               title = "GSEA Enrichment Plot: HECKER_IFNB1_TARGETS")
### view the plot
g
### export the plot to figures directory
ggsaver(
  "gsea_enrich_plot-hecker_ifn_target_genes-irf3_vs_bfp.add_for_ppr", g, path=outputs$figures,
  height=15, width=16, units="cm")
## cleanup
rm(ir_res, ir_rank, heck_ens, gene_set, g)


# plot volcano plot of pluripotency genes (IRF3vBFP)
## set up and process full deseq results table for IRF3vsBFP for volcano plot
### specify the table
ir_res <- pl_res[["irf35d-vs-bfp"]]
### add a column ("type") specifying classification of each gene
ir_res <- ir_res %>%
  mutate(type = case_when(log2FoldChange >= 0 & padj <= 0.01 ~ "up",
                          log2FoldChange <= 0 & padj <= 0.01 ~ "down",
                          TRUE ~ "ns"))   
  ### view the output
  head(ir_res)
  tail(ir_res)
## subset the deseq results to only contain genes of interest
plg_res <- ir_res[ir_res$gene_symbol %in% plg, ]
  ### check output
  dim(plg_res)
  length(unique(plg_res$gene_symbol))
  length(plg)
## set colors for displaying the points
cols <- c("up" = "steelblue", "down" = "steelblue", "ns" = "grey")
## plot the volcano plot
g <- ggplot(data = ir_res, aes(x = log2FoldChange, y = -log10(padj))) +
  xlim(c(-6, 6)) + 
  ylim(c(0, 20)) +
  geom_point(
    aes(color = type),
    alpha = 0.5, size = 1, shape = 16) +  ## set transparency and size on those which aren't of interest (background points)
  #scale_fill_manual(values = cols) + 
  #scale_color_manual(values = cols) +
  #geom_point(colour = "grey", alpha = 0.8) + 
  geom_point(data = plg_res,  ## add new layer containing genes of interest 
             shape = 21,
             size = 2,
             #alpha = 1,
             color = "firebrick",
             fill = "firebrick") +  
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed") + 
  # geom_vline(xintercept = c(log2(0.5), log2(2)),  ## for adding a logfc threshold line
  #            linetype = "dashed") +
  labs(title = "Expression of pluripotency genes - irf3 vs bfp",
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  geom_label_repel(data = plg_res, # Add labels last to appear as the top layer  
                   aes(label = gene_symbol),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
  #                    limits = c(-10, 10)) +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.title=element_text(size=14),
        legend.position = "none") +
  theme_classic()
### view the plot
g
### save the output to figures directory
ggsaver(
  "volcano-pluripotency_genes-irf3_vs_bfp.add_for_ppr", g, path=outputs$figures,
  height=15, width=16, units="cm")
## cleanup
rm(cols, ir_res, plg_res, g)


# plot volcano plot of IFN and pathway relevant genes
## set up and process full deseq results table for IRF3vsBFP for volcano plot
### specify the table
ir_res <- pl_res[["irf35d-vs-bfp"]]
### add a column ("type") specifying classification of each gene
ir_res <- ir_res %>%
  mutate(type = case_when(log2FoldChange >= 0 & padj <= 0.01 ~ "up",
                          log2FoldChange <= 0 & padj <= 0.01 ~ "down",
                          TRUE ~ "ns"))   
### view the output
head(ir_res)
tail(ir_res)
## subset the deseq results to only contain genes of interest
irgen_res <- ir_res[ir_res$gene_id %in% irgen, ]
  ### check output
  dim(irgen_res)
  length(unique(irgen_res$gene_symbol))
  length(irgen)
## set colors for displaying the points
cols <- c("up" = "steelblue", "down" = "steelblue", "ns" = "grey")
## plot the volcano plot
g <- ggplot(data = ir_res, aes(x = log2FoldChange, y = -log10(padj))) +
  xlim(c(-6, 6)) + 
  ylim(c(0, 20)) +
  geom_point(
    aes(color = type),
    alpha = 0.5, size = 1, shape = 16) +  ## set transparency and size on those which aren't of interest (background points)
  #scale_fill_manual(values = cols) + 
  #scale_color_manual(values = cols) +
  #geom_point(colour = "grey", alpha = 0.8) + 
  geom_point(data = irgen_res,  ## add new layer containing genes of interest 
             shape = 21,
             size = 2,
             #alpha = 1,
             color = "firebrick",
             fill = "firebrick") +  
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed") + 
  # geom_vline(xintercept = c(log2(0.5), log2(2)),  ## for adding a logfc threshold line
  #            linetype = "dashed") +
  labs(title = "Expression of relevant genes - irf3 vs bfp",
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  geom_label_repel(data = irgen_res, # Add labels last to appear as the top layer  
                   aes(label = gene_symbol),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
  #                    limits = c(-10, 10)) +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.title=element_text(size=14),
        legend.position = "none") +
  theme_classic()
### view the plot
g
### save the output to figures directory
ggsaver(
  "volcano-other_relevant_genes_prev_ma-irf3_vs_bfp.add_for_ppr", g, path=outputs$figures,
  height=15, width=16, units="cm")
## cleanup
rm(cols, ir_res, irgen_res, g)


# plot heatmap of apoptosis genes - with all samples across all treatments
## subset expression matrix to include only genes of interest
### duplicate exprs and remove summarized results, keeping annotation info
exp_g <- exprs[, colnames(exprs) %in% c("gene_id", 
                                        as.character(metadata$sample), 
                                        "gene_symbol", "gene_name")]
### subset the expression matrix to only include the target genes of interest
exp_g <- exp_g[exp_g$gene_symbol %in% apop, ]
  ### check output
  length(apop)
  dim(exp_g)
  length(unique(exp_g$gene_symbol))
  anyDuplicated(exp_g$gene_symbol)  ## check for duplicated gene symbols
## remove row with duplicate value ‘Pik3r3’ (removing ENSMUSG00000111410, keeping ENSMUSG00000028698, both are named Pik3r3 with same gene description)
exp_g <- exp_g %>%
  filter(gene_id != "ENSMUSG00000111410")
## convert rows to gene_symbol for displaying gene symbols on heatmap
exp_g_df <- as.data.frame(exp_g) ## convert to a data frame (required for converting rownames)
rownames(exp_g_df) <- exp_g_df$gene_symbol ## convert rows to gene symbol
### check output
head(exp_g_df)
## calculate z-scores
### remove annotation info for z-score conversion (all column values need to be numeric)
exp_g_df <- exp_g_df %>% dplyr::select(-gene_id, -gene_symbol, -gene_name)
### calculate z-scores
z_exp_g <- t(scale(t(exp_g_df)))
### check output
head(z_exp_g)
## set up annotation bar for heatmap
metadata_hanno <- metadata[, "treatment", drop = FALSE] ## extract the `treatment` column for annotation, and maintain as df
### make named vectors for annotation bar colors
treatment_levels <- unique(as.character(metadata_hanno$treatment))
treatment_colors <- setNames(colorRampPalette(brewer.pal(2, "Set3"))(length(treatment_levels)), treatment_levels)
### design the heatmap annotation
h_ann = ComplexHeatmap::HeatmapAnnotation(df = metadata_hanno, 
                                          which="column",
                                          col = list(treatment = treatment_colors))
## plot the heatmap
g = pheatmap(
  z_exp_g,  ## specify order of rows (to match prior figure for comparison), requires no row clustering
  top_annotation = h_ann, ## specify top annotation bar
  name = "z-score", ## set name above heatmap scale
  cluster_rows = F, ## cluster rows (genes)
  cluster_cols = T, ## cluster columns (samples)
  show_rownames = T, ## show gene names
  show_colnames = T, ## show sample names
  border_color = NA, ## remove borders between cells
  #color = colorRampPalette(brewer.pal(11, "RdBu"))(11), 
  color = hCol(28),
  #color = colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(50),
  main = "Expression of Apoptosis Genes (KEGG)") ## set title
### view the plot
plot(g)
## plot the heatmap
g = pheatmap(
  z_exp_g,  ## specify order of rows (to match prior figure for comparison), requires no row clustering
  top_annotation = h_ann, ## specify top annotation bar
  name = "z-score", ## set name above heatmap scale
  cluster_rows = F, ## cluster rows (genes)
  cluster_cols = T, ## cluster columns (samples)
  show_rownames = T, ## show gene names
  show_colnames = T, ## show sample names
  border_color = NA, ## remove borders between cells
  #color = colorRampPalette(brewer.pal(11, "RdBu"))(11), 
  color = hCol(28),
  #color = colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(50),
  main = "Expression of Apoptosis Genes (KEGG)") ## set title
### view the plot
plot(g)
### export the final heatmap to figures directory
pdf(here(outputs$figures.do, "heatmap_genes-apoptosis_kegg-all_treatments.add_for_ppr.pdf"),
    width = cm_to_inch(22), height = cm_to_inch(30)
)
  plot(g)
dev.off()
png(here(outputs$figures.do, "heatmap_genes-apoptosis_kegg-all_treatments.add_for_ppr.png"),
    width = 22, height = 30, units = "cm", res = 150
)
  plot(g)
dev.off()

# plot heatmap of apoptosis genes
## subset expression matrix to include only genes of interest
### duplicate exprs and remove summarized results and new samples, keeping annotation info
exp_g <- exprs[, colnames(exprs) %in% c("gene_id", 
                                        as.character(metadata_old$sample), 
                                        "gene_symbol", "gene_name")]
### subset the expression matrix to only include the target genes of interest
exp_g <- exp_g[exp_g$gene_symbol %in% apop, ]
### check output
length(apop)
dim(exp_g)
length(unique(exp_g$gene_symbol))
anyDuplicated(exp_g$gene_symbol)  ## check for duplicated gene symbols
## remove row with duplicate value ‘Pik3r3’ (removing ENSMUSG00000111410, keeping ENSMUSG00000028698, both are named Pik3r3 with same gene description)
exp_g <- exp_g %>%
  filter(gene_id != "ENSMUSG00000111410")
## convert rows to gene_symbol for displaying gene symbols on heatmap
exp_g_df <- as.data.frame(exp_g) ## convert to a data frame (required for converting rownames)
rownames(exp_g_df) <- exp_g_df$gene_symbol ## convert rows to gene symbol
### check output
head(exp_g_df)
## calculate z-scores
### remove annotation info for z-score conversion (all column values need to be numeric)
exp_g_df <- exp_g_df %>% dplyr::select(-gene_id, -gene_symbol, -gene_name)
### calculate z-scores
z_exp_g <- t(scale(t(exp_g_df)))
### check output
head(z_exp_g)
## set up annotation bar for heatmap
metadata_hanno <- metadata_old[, "treatment", drop = FALSE] ## extract the `treatment` column for annotation, and maintain as df
### make named vectors for annotation bar colors
treatment_levels <- unique(as.character(metadata_hanno$treatment))
treatment_colors <- setNames(colorRampPalette(brewer.pal(2, "Set3"))(length(treatment_levels)), treatment_levels)
### design the heatmap annotation
h_ann = ComplexHeatmap::HeatmapAnnotation(df = metadata_hanno, 
                                          which="column",
                                          col = list(treatment = treatment_colors))
## plot the heatmap
g = pheatmap(
  z_exp_g,  ## specify order of rows (to match prior figure for comparison), requires no row clustering
  top_annotation = h_ann, ## specify top annotation bar
  name = "z-score", ## set name above heatmap scale
  cluster_rows = F, ## cluster rows (genes)
  cluster_cols = T, ## cluster columns (samples)
  show_rownames = T, ## show gene names
  show_colnames = T, ## show sample names
  border_color = NA, ## remove borders between cells
  #color = colorRampPalette(brewer.pal(11, "RdBu"))(11), 
  color = hCol(16),
  #color = colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(50),
  main = "Expression of Apoptosis Genes (KEGG)") ## set title
### view the plot
plot(g)
### export the final heatmap to figures directory
pdf(here(outputs$figures.do, "heatmap_genes-apoptosis_kegg-old_datasets.add_for_ppr.pdf"),
    width = cm_to_inch(18), height = cm_to_inch(30)
)
  plot(g)
dev.off()
png(here(outputs$figures.do, "heatmap_genes-apoptosis_kegg-old_datasets.add_for_ppr.png"),
    width = 18, height = 30, units = "cm", res = 150
)
  plot(g)
dev.off()
## cleanup
rm(exp_g, exp_g_df, z_exp_g, 
   g, metadata_hanno, h_ann, 
   treatment_levels, treatment_colors)


# plotting additional heatmaps of apoptosis genes with hierarchical clustering of both columns and rows
## plotting with all samples across all treatments
### subset expression matrix to include only genes of interest
#### duplicate exprs and remove summarized results, keeping annotation info
exp_g <- exprs[, colnames(exprs) %in% c("gene_id", 
                                        as.character(metadata$sample), 
                                        "gene_symbol", "gene_name")]
#### subset the expression matrix to only include the target genes of interest
exp_g <- exp_g[exp_g$gene_symbol %in% apop, ]
  #### check output
  length(apop)
  dim(exp_g)
  length(unique(exp_g$gene_symbol))
  anyDuplicated(exp_g$gene_symbol)  ## check for duplicated gene symbols
### remove row with duplicate value ‘Pik3r3’ (removing ENSMUSG00000111410, keeping ENSMUSG00000028698, both are named Pik3r3 with same gene description)
exp_g <- exp_g %>%
  filter(gene_id != "ENSMUSG00000111410")
### convert rows to gene_symbol for displaying gene symbols on heatmap
exp_g_df <- as.data.frame(exp_g) ## convert to a data frame (required for converting rownames)
rownames(exp_g_df) <- exp_g_df$gene_symbol ## convert rows to gene symbol
  #### check output
  head(exp_g_df)
### calculate z-scores
#### remove annotation info for z-score conversion (all column values need to be numeric)
exp_g_df <- exp_g_df %>% dplyr::select(-gene_id, -gene_symbol, -gene_name)
#### calculate z-scores
z_exp_g <- t(scale(t(exp_g_df)))
  #### check output
  head(z_exp_g)
### set up annotation bar for heatmap
metadata_hanno <- metadata[, "treatment", drop = FALSE] ## extract the `treatment` column for annotation, and maintain as df
#### make named vectors for annotation bar colors
treatment_levels <- unique(as.character(metadata_hanno$treatment))
treatment_colors <- setNames(colorRampPalette(brewer.pal(2, "Set3"))(length(treatment_levels)), treatment_levels)
#### design the heatmap annotation
h_ann = ComplexHeatmap::HeatmapAnnotation(df = metadata_hanno, 
                                          which="column",
                                          col = list(treatment = treatment_colors))
### plot the heatmap
g = pheatmap(
  z_exp_g,  ## specify order of rows (to match prior figure for comparison), requires no row clustering
  top_annotation = h_ann, ## specify top annotation bar
  name = "z-score", ## set name above heatmap scale
  cluster_rows = T, ## cluster rows (genes)
  cluster_cols = T, ## cluster columns (samples)
  show_rownames = T, ## show gene names
  show_colnames = T, ## show sample names
  border_color = NA, ## remove borders between cells
  #color = colorRampPalette(brewer.pal(11, "RdBu"))(11), 
  color = hCol(28),
  #color = colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(50),
  main = "Expression of Apoptosis Genes (KEGG)") ## set title
#### view the plot
plot(g)
#### export the final heatmap to figures directory
pdf(here(outputs$figures.do, "heatmap_genes-apoptosis_kegg-all_treatments-rows_clustered.add_for_ppr.pdf"),
    width = cm_to_inch(22), height = cm_to_inch(30)
)
plot(g)
dev.off()
png(here(outputs$figures.do, "heatmap_genes-apoptosis_kegg-all_treatments-rows_clustered.add_for_ppr.png"),
    width = 22, height = 30, units = "cm", res = 150
)
plot(g)
dev.off()

## plotting with only old dataset samples
### subset expression matrix to include only genes of interest
#### duplicate exprs and remove summarized results and new samples, keeping annotation info
exp_g <- exprs[, colnames(exprs) %in% c("gene_id", 
                                        as.character(metadata_old$sample), 
                                        "gene_symbol", "gene_name")]
#### subset the expression matrix to only include the target genes of interest
exp_g <- exp_g[exp_g$gene_symbol %in% apop, ]
  #### check output
  length(apop)
  dim(exp_g)
  length(unique(exp_g$gene_symbol))
  anyDuplicated(exp_g$gene_symbol)  ## check for duplicated gene symbols
### remove row with duplicate value ‘Pik3r3’ (removing ENSMUSG00000111410, keeping ENSMUSG00000028698, both are named Pik3r3 with same gene description)
exp_g <- exp_g %>%
  filter(gene_id != "ENSMUSG00000111410")
### convert rows to gene_symbol for displaying gene symbols on heatmap
exp_g_df <- as.data.frame(exp_g) ## convert to a data frame (required for converting rownames)
rownames(exp_g_df) <- exp_g_df$gene_symbol ## convert rows to gene symbol
  #### check output
  head(exp_g_df)
### calculate z-scores
#### remove annotation info for z-score conversion (all column values need to be numeric)
exp_g_df <- exp_g_df %>% dplyr::select(-gene_id, -gene_symbol, -gene_name)
#### calculate z-scores
z_exp_g <- t(scale(t(exp_g_df)))
  #### check output
  head(z_exp_g)
### set up annotation bar for heatmap
metadata_hanno <- metadata_old[, "treatment", drop = FALSE] ## extract the `treatment` column for annotation, and maintain as df
#### make named vectors for annotation bar colors
treatment_levels <- unique(as.character(metadata_hanno$treatment))
treatment_colors <- setNames(colorRampPalette(brewer.pal(2, "Set3"))(length(treatment_levels)), treatment_levels)
#### design the heatmap annotation
h_ann = ComplexHeatmap::HeatmapAnnotation(df = metadata_hanno, 
                                          which="column",
                                          col = list(treatment = treatment_colors))
### plot the heatmap
g = pheatmap(
  z_exp_g,  ## specify order of rows (to match prior figure for comparison), requires no row clustering
  top_annotation = h_ann, ## specify top annotation bar
  name = "z-score", ## set name above heatmap scale
  cluster_rows = T, ## cluster rows (genes)
  cluster_cols = T, ## cluster columns (samples)
  show_rownames = T, ## show gene names
  show_colnames = T, ## show sample names
  border_color = NA, ## remove borders between cells
  #color = colorRampPalette(brewer.pal(11, "RdBu"))(11), 
  color = hCol(16),
  #color = colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(50),
  main = "Expression of Apoptosis Genes (KEGG)") ## set title
#### view the plot
plot(g)
#### export the final heatmap to figures directory
pdf(here(outputs$figures.do, "heatmap_genes-apoptosis_kegg-old_datasets-rows_clustered.add_for_ppr.pdf"),
    width = cm_to_inch(18), height = cm_to_inch(30)
)
plot(g)
dev.off()
png(here(outputs$figures.do, "heatmap_genes-apoptosis_kegg-old_datasets-rows_clustered.add_for_ppr.png"),
    width = 18, height = 30, units = "cm", res = 150
)
plot(g)
dev.off()
## cleanup
rm(exp_g, exp_g_df, z_exp_g, 
   g, metadata_hanno, h_ann, 
   treatment_levels, treatment_colors)


# plot volcano plot
## set up and process full deseq results table for IRF3vsBFP for volcano plot
### specify the table
ir_res <- pl_res[["irf35d-vs-bfp"]]
### add a column ("type") specifying classification of each gene
ir_res <- ir_res %>%
  mutate(type = case_when(log2FoldChange >= 0 & padj <= 0.01 ~ "up",
                          log2FoldChange <= 0 & padj <= 0.01 ~ "down",
                          TRUE ~ "ns"))   
### view the output
head(ir_res)
tail(ir_res)
## subset the deseq results to only contain genes of interest
plg_res <- ir_res[ir_res$gene_symbol %in% plg_ifn, ]
### check output
dim(plg_res)
length(unique(plg_res$gene_symbol))
length(plg)
## set colors for displaying the points
cols <- c("up" = "steelblue", "down" = "steelblue", "ns" = "grey")
## plot the volcano plot - with axis limits
g <- ggplot(data = ir_res, aes(x = log2FoldChange, y = -log10(padj))) +
  xlim(c(-6, 6)) + 
  ylim(c(0, 20)) +
  geom_point(
    aes(color = type),
    alpha = 0.5, size = 1, shape = 16) +  ## set transparency and size on those which aren't of interest (background points)
  #scale_fill_manual(values = cols) + 
  #scale_color_manual(values = cols) +
  #geom_point(colour = "grey", alpha = 0.8) + 
  geom_point(data = plg_res,  ## add new layer containing genes of interest 
             shape = 21,
             size = 2,
             #alpha = 1,
             color = "firebrick",
             fill = "firebrick") +  
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed") + 
  # geom_vline(xintercept = c(log2(0.5), log2(2)),  ## for adding a logfc threshold line
  #            linetype = "dashed") +
  labs(title = "Expression of pluripotency genes and IFN UP genes - irf3 vs bfp",
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  geom_label_repel(data = plg_res, # Add labels last to appear as the top layer  
                   aes(label = gene_symbol),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
  #                    limits = c(-10, 10)) +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.title=element_text(size=14),
        legend.position = "none") +
  theme_classic()
### view the plot
g
### save the output to figures directory - with axis limits
ggsaver(
  "volcano-pluripotency_genes_ifn_up_genes-irf3_vs_bfp.add_for_ppr", g, path=outputs$figures,
  height=15, width=16, units="cm")
## plot the volcano plot - with NO axis limits, to see all points
g <- ggplot(data = ir_res, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(c(-6, 6)) + 
  #ylim(c(0, 20)) +
  geom_point(
    aes(color = type),
    alpha = 0.5, size = 1, shape = 16) +  ## set transparency and size on those which aren't of interest (background points)
  #scale_fill_manual(values = cols) + 
  #scale_color_manual(values = cols) +
  #geom_point(colour = "grey", alpha = 0.8) + 
  geom_point(data = plg_res,  ## add new layer containing genes of interest 
             shape = 21,
             size = 2,
             #alpha = 1,
             color = "firebrick",
             fill = "firebrick") +  
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed") + 
  # geom_vline(xintercept = c(log2(0.5), log2(2)),  ## for adding a logfc threshold line
  #            linetype = "dashed") +
  labs(title = "Expression of pluripotency genes and IFN UP genes - irf3 vs bfp",
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  geom_label_repel(data = plg_res, ## add labels last to appear as the top layer  
                   aes(label = gene_symbol),
                   force = 2,
                   nudge_y = 1,
                   max.overlaps = Inf, ## ensure all points appear
                   box.padding = 0.3, 
                   point.padding = 0.3,
                   segment.size = 0.2,
                   min.segment.length = 0) +
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
  #                    limits = c(-10, 10)) +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.title=element_text(size=14),
        legend.position = "none") +
  theme_classic()
### view the plot
g
### save the output to figures directory - with axis limits
ggsaver(
  "volcano-pluripotency_genes_ifn_up_genes-irf3_vs_bfp-no_axis_limits.add_for_ppr", g, path=outputs$figures,
  height=15, width=16, units="cm")

## plot the volcano plot - with points indicated out of bounds
### set axis limits
xlim_vals <- c(-6, 10)
ylim_vals <- c(0, 25)
### cap all -log10 padj values of the genes of interest to y-max (for plotting the points at set values)
#### add a new column with -log10 transformed pvalues
plg_res$neg_log10_pval <- -log10(plg_res$padj)
#### copy transformed results into new columns to contain capped values
plg_res$cap_log10_pval <- plg_res$neg_log10_pval
#### cap the values to y-max
plg_res$cap_log10_pval <- pmin(plg_res$cap_log10_pval, ylim_vals[2])
  #### check the output
  head(plg_res)
  min(plg_res$neg_log10_pval)
  max(plg_res$neg_log10_pval)
  min(plg_res$cap_log10_pval)
  max(plg_res$cap_log10_pval)
### split results into two dfs based on value of -log10 padj, to plot in different layers and different shapes on the volcano plot
plg_more <- plg_res[plg_res$neg_log10_pval > ylim_vals[2], ]
plg_less <- plg_res[plg_res$neg_log10_pval <= ylim_vals[2], ]
  ### check output
  head(plg_more)
  head(plg_less)
  dim(plg_more)
  dim(plg_less)
  dim(plg_res)
  max(plg_less$neg_log10_pval)
  max(plg_more$neg_log10_pval)
### plot the volcano plot - with arrows to indicate points out-of-bounds
g <- ggplot(data = ir_res, aes(x = log2FoldChange, y = -log10(padj))) +
  xlim(xlim_vals) + 
  ylim(ylim_vals) +
  geom_point(
    aes(color = type),
    alpha = 0.5, size = 1, shape = 16) +  
  #scale_fill_manual(values = cols) + 
  #scale_color_manual(values = cols) +
  #geom_point(colour = "grey", alpha = 0.8) + 
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed") + 
  # geom_vline(xintercept = c(log2(0.5), log2(2)),  ## for adding a logfc threshold line
  #            linetype = "dashed") +
  labs(title = "Expression of pluripotency genes and IFN UP genes - irf3 vs bfp",
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  geom_label_repel(data = plg_res, ## add labels for in-bounds genes
                   aes(label = gene_symbol),
                   force = 2,
                   nudge_y = 1,
                   segment.alpha = 1,
                   min.segment.length = 0.01) +
  geom_label_repel(data = plg_more, ## add labels for out-of-bounds too
                   aes(x = log2FoldChange,
                       y = cap_log10_pval,
                       label = gene_symbol),
                   force = 2,
                   nudge_y = -1,
                   nudge_x = 0.2,
                   segment.alpha = 1,
                   min.segment.length = 0.01) +
  geom_point(data = plg_less,  ## add new layer containing genes of interest (in-bounds genes)
             aes(x = log2FoldChange,
                 y = neg_log10_pval),
             shape = 21,
             size = 2,
             #alpha = 1,
             color = "firebrick",
             fill = "firebrick") +  
  geom_point(data = plg_more,  ## add new layer containing genes of interest (out-of-bounds genes)
             aes(x = log2FoldChange,
                 y = cap_log10_pval),
             shape = 24,
             size = 2.5,
             #alpha = 1,
             color = "firebrick",
             fill = "firebrick") +  
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
  #                    limits = c(-10, 10)) +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.title=element_text(size=14),
        legend.position = "none") +
  theme_classic()
### view the plot
g
### save the output to figures directory - with out of bounds arrows
ggsaver(
  "volcano-pluripotency_genes_ifn_up_genes-irf3_vs_bfp-arrows_bounds.add_for_ppr", g, path=outputs$figures,
  height=15, width=16, units="cm")

# volcano version
## subset the deseq results to only contain genes of interest
### pluripotency genes
pl_plg <- ir_res[ir_res$gene_symbol %in% plg, ]
### IFN only genes
pl_ifn <- ir_res[ir_res$gene_symbol %in% p_ifn, ]
## plot the volcano plot - with axis limits and arrows to indicate out of bounds point locations for labeled genes
#### add a new column with -log10 transformed pvalues
pl_ifn$neg_log10_pval <- -log10(pl_ifn$padj)
#### copy transformed results into new columns to contain capped values
pl_ifn$cap_log10_pval <- pl_ifn$neg_log10_pval
#### cap the values to y-max
pl_ifn$cap_log10_pval <- pmin(pl_ifn$cap_log10_pval, ylim_vals[2])
### split results into two dfs based on value of -log10 padj, to plot in different layers and different shapes on the volcano plot
plg_more <- pl_ifn[pl_ifn$neg_log10_pval > ylim_vals[2], ]
plg_less <- pl_ifn[pl_ifn$neg_log10_pval <= ylim_vals[2], ]
### plot the volcano plot - with arrows to indicate points out-of-bounds and with gene groups colored separately
g <- ggplot(data = ir_res, aes(x = log2FoldChange, y = -log10(padj))) +
  xlim(xlim_vals) + 
  ylim(ylim_vals) +
  geom_point(
    aes(color = type),
    alpha = 0.5, size = 1, shape = 16) +  ## set transparency and size on those which aren't of interest (background points)
  #scale_fill_manual(values = cols) + 
  #scale_color_manual(values = cols) +
  #geom_point(colour = "grey", alpha = 0.8) + 
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed") + 
  # geom_vline(xintercept = c(log2(0.5), log2(2)),  ## for adding a logfc threshold line
  #            linetype = "dashed") +
  labs(title = "Expression of pluripotency genes and IFN UP genes - irf3 vs bfp",
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  geom_label_repel(data = pl_plg, ## add labels pluripotency genes
                   aes(label = gene_symbol),
                   force = 2,
                   nudge_y = 1,
                   segment.alpha = 1,
                   min.segment.length = 0.01) +
  geom_label_repel(data = pl_ifn, ## add labels for in-bounds IFN genes
                   aes(label = gene_symbol),
                   force = 2,
                   nudge_y = 1,
                   segment.alpha = 1,
                   min.segment.length = 0.01) +
  geom_label_repel(data = plg_more, ## add labels for out-of-bounds IFN genes too
                   aes(x = log2FoldChange,
                       y = cap_log10_pval,
                       label = gene_symbol),
                   force = 2,
                   nudge_y = -1,
                   nudge_x = 0.2,
                   segment.alpha = 1,
                   min.segment.length = 0.01) +
  geom_point(data = pl_plg,  ## add new layer containing genes of interest - pluripotency genes
             aes(x = log2FoldChange,
                 y = -log10(padj)),
             shape = 21,
             size = 2,
             #alpha = 1,
             color = "blueviolet",
             fill = "blueviolet") +  
  geom_point(data = plg_less,  ## add new layer containing genes of interest - in-bounds IFN genes
             aes(x = log2FoldChange,
                 y = neg_log10_pval),
             shape = 21,
             size = 2,
             #alpha = 1,
             color = "firebrick",
             fill = "firebrick") +  
  geom_point(data = plg_more,  ## add new layer containing genes of interest - out-of-bounds IFN genes
             aes(x = log2FoldChange,
                 y = cap_log10_pval),
             shape = 24,
             size = 2.5,
             #alpha = 1,
             color = "firebrick",
             fill = "firebrick") +  
  scale_colour_manual(values = cols) + 
  # scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
  #                    limits = c(-10, 10)) +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.title=element_text(size=14),
        legend.position = "none") +
  theme_classic()
### view the plot
g
ggsaver(
  "volcano-pluripotency_genes_ifn_up_genes-irf3_vs_bfp-arrows_bounds-gene_groups_cols.add_for_ppr", g, path=outputs$figures,
  height=15, width=16, units="cm")
## cleanup
rm(cols, ir_res, plg_res, plg_less, plg_more, pl_plg, pl_ifn, xlim_vals, ylim_vals,
   g)

# plot MA plots of STING GOF vs BFP displaying genes (irgen) 
## prepare for plotting - set colors for displaying the points of interest
cols <- c("up" = "steelblue", "down" = "steelblue", "ns" = "grey")
## set up and process full deseq results table for stingGOF-VS-bfp for volcano plot, prepare the gene set and link to expression data
### specify the table
ir_res <- pl_res[["sting_gof-vs-bfp"]]
### specify the deseq results
ir_res_o <- pl_res_o[["sting_gof-vs-bfp"]]
  ### view the output
  head(ir_res)
  tail(ir_res)
## process the irgen list to retrieve gene symbols, not just gene ids, and extract into new vector (for ease of gene selection and filtration)
setdiff(irgen, exprs$gene_id) ## check that all are present in exprs
names <- exprs$gene_symbol[match(irgen, exprs$gene_id)]
  ### check output
  names
  length(names)
  length(irgen)

# prepare and plot MA plot - pluripotency genes removed
## process and subset the gene list to only contain genes of interest, create irgen subsets, filter the irgen results list
### remove pluripotency genes
  ### remove lif, klf4, klf5, nanog, (tbx3???), tcf3 ???, essrb, pou5f1, zfp42, klf2, sox2
#### create vector of gene names to remove
rm_list1 <- c(
  "Lif",
  "Klf4",
  "Klf5",
  "Nanog",
  "Tbx3",
  "Tcf3",
  "Esrrb",
  "Pou5f1",
  "Zfp42",
  "Klf2",
  "Sox2"
)
#### filter the irgen list to remove rows for these genes, subset the deseq results to only contain genes of interest - removing pluripotency genes
irgen_f <- names[!names %in% rm_list1]
  #### check the output
  print(irgen_f)
  length(irgen_f)
  length(irgen)
  length(rm_list1)
  length(setdiff(rm_list1, irgen_f)) ## check that none exist in the main df
## ADD IN TBK1, maybe useful, to the genes to plot
irgen_f <- c(irgen_f, "Tbk1")
## run plotMA on deseq2 results, and process to add gene information
ma_data = DESeq2::plotMA(ir_res_o, alpha=stats$de_alpha, returnData=T)
ma_data$gene_id <- ir_res$gene_id  ## add gene_id column to the MA plot data
ma_data$gene_symbol <- ir_res$gene_symbol  ## add gene_symbol column to the MA plot data
## log2 transform the mean of normalized counts in a new column (to match existing figure)
ma_data$mean[ma_data$mean == 0] <- 1e-6 
ma_data$log2_mean <- log2(ma_data$mean + 1)
## prepare y-limit
y_limit = ceiling(max(abs(quantile(ma_data$lfc, c(0.001, 0.999)))))
## plot MA plot using log2 mean of normalized counts
g <- ggplot2::ggplot(ma_data, 
                     aes(mean, lfc, color=isDE)) +
  ggplot2::geom_point(alpha=0.33, size=1, stroke=0) +
  ggplot2::geom_hline(yintercept=c(-stats$de_l2FC_old, stats$de_l2FC_old)) +
  ggplot2::scale_color_manual(values=c(`TRUE`="blue", `FALSE`="grey75")) +
  #ggplot2::coord_cartesian(ylim=c(-y_limit, y_limit)) +  ## removing upper cutoff so as not to cut off genes
  ggplot2::scale_x_continuous(trans="log10") +   ## apply log10 transformation to mean of normalized counts
  ggplot2::labs(
    x="Mean of normalized counts",
    y="Shrunken log2 Fold-Change",
    title=glue("sting_gof-vs-bfp (genes)"),
    caption=glue('DESeq results: padj 0.01, logFC 0')) +
  ggrepel::geom_label_repel(
    data = subset(ma_data, gene_symbol %in% irgen_f),
    aes(mean, lfc, label = gene_symbol),
    size = 3,
    color = "purple",
    min.segment.length = 0,  ## draw lines from labels to points
    max.overlaps = Inf,
    force = 2              ## increase force to repel labels more strongly
  ) +
  #ggplot2::geom_text(data = label_data, aes(label = label), vjust = 1.5, hjust = 0.5, color = "red") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position="none")
plot(g)
## export and save plots
ggsaver(
  glue("de-MAplot_genes-ifn_relevant_genes_no_plurip-sting_gof_vs_bfp.add_for_ppr"), g, path=outputs$figures,
  height=cm_to_inch(10), width=cm_to_inch(10), units="in",
  device=c("png", "pdf"))
## PARTIAL cleanup
rm(g, ma_data, y_limit)

# prepare and plot MA plot - filtered gene set
## process and subset the gene list to only contain genes of interest, create irgen subsets, filter the irgen results list
### remove additional uninteresting genes
  ### removing those with low counts
#### create vector of gene names to remove
rm_list2 <- c(
  "Ifnb1",
  "Ifna2",
  "Ifna1",
  "Ifnab",
  "Ifna4",
  "Ifna9",
  "Ifna5"
)
#### filter the irgen list to remove rows for these genes, subset the deseq results to only contain genes of interest - removing pluripotency genes
irgen_af <- irgen_f[!irgen_f %in% rm_list2]
  #### check the output
  print(irgen_af)
  length(irgen_af)
  length(irgen_f)
  length(rm_list2)
  length(setdiff(rm_list2, irgen_af)) ## check that none exist in the main df
## run plotMA on deseq2 results, and process to add gene information
ma_data = DESeq2::plotMA(ir_res_o, alpha=stats$de_alpha, returnData=T)
ma_data$gene_id <- ir_res$gene_id  ## add gene_id column to the MA plot data
ma_data$gene_symbol <- ir_res$gene_symbol  ## add gene_symbol column to the MA plot data
## log2 transform the mean of normalized counts in a new column (to match existing figure)
ma_data$mean[ma_data$mean == 0] <- 1e-6 ## doing this to ensure all genes are included and prevent infinite values for mean of 0 (primarily for labeling purposes, even if they're at zero)
ma_data$log2_mean <- log2(ma_data$mean + 1)
## prepare y-limit
#y_limit = ceiling(max(abs(quantile(ma_data$lfc, c(0.001, 0.999)))))
## plot MA plot using log2 mean of normalized counts
g <- ggplot2::ggplot(ma_data, 
                     aes(mean, lfc, color=isDE)) +
  #ylim(c(-1.6, 1.6)) +
  ggplot2::geom_point(alpha=0.33, size=1, stroke=0) +
  ggplot2::geom_hline(yintercept=c(-stats$de_l2FC_old, stats$de_l2FC_old)) +
  ggplot2::scale_color_manual(values=c(`TRUE`="blue", `FALSE`="grey75")) +
  #ggplot2::coord_cartesian(ylim=c(-y_limit, y_limit)) + 
  ggplot2::scale_x_continuous(trans="log10") +   ## apply log10 transformation to mean of normalized counts
  ggplot2::coord_cartesian(xlim = c(0.1, 100000), ## zoom in
                           ylim = c(-1.6, 1.6)) +
  ggplot2::labs(
    x="Mean of normalized counts",
    y="Shrunken log2 Fold-Change",
    title=glue("sting_gof-vs-bfp (genes)"),
    caption=glue('DESeq results: padj 0.01, logFC 0')) +
  ggrepel::geom_label_repel(
    data = subset(ma_data, gene_symbol %in% irgen_af),
    aes(mean, lfc, label = gene_symbol),
    size = 3,
    color = "purple",
    min.segment.length = 0,  ## draw lines from labels to points
    max.overlaps = Inf,
    force = 2              ## increase force to repel labels more strongly
  ) +
  #ggplot2::geom_text(data = label_data, aes(label = label), vjust = 1.5, hjust = 0.5, color = "red") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position="none")
plot(g)
## export and save plots
ggsaver(
  glue("de-MAplot_genes-ifn_relevant_genes_no_plurip_filtered-sting_gof_vs_bfp.add_for_ppr"), g, path=outputs$figures,
  height=cm_to_inch(10), width=cm_to_inch(10), units="in",
  device=c("png", "pdf"))
## cleanup
rm(cols, ir_res, ir_res_o, g,
   ma_data, rm_list1, rm_list2, names, 
   irgen_f, irgen_af)
