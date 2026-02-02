library(optparse)

option_list <- list(
    optparse::make_option(
        c("--config"), type="character", metavar="character",
        default = 'config.yml',
        help="the configuration .yml file containing all the settings for running the differential analysis snakemake"
    ),
    optparse::make_option(
        c("--table"), type="character", metavar="character",
        default = 'config.tsv',
        help="the .tsv file containing information about the different comparisons to run'"
    ),
    optparse::make_option(
        c("--comparison"), type="character", metavar="character",
        help="the name of the current comparison to perform the analysis on"
    ),
    optparse::make_option(
        c("--cluster"), type="character", metavar="character",
        help="the name of the current cluster to perform the analysis on"
    ),
    optparse::make_option(
        c("--test"), type="character", metavar="character",
        help='the name of the current DEG test, corresponding to a named entry inside config["tests"]'
    ),
    optparse::make_option(
        c("--outdir"), type="character", metavar="character",
        help="path to output directory."
    )
)

# parse the arguments:
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)

# --------------------------------------------------------------- #
# load libraries and configuration file
# --------------------------------------------------------------- #

library(yaml)
library(tidyverse)
library(Seurat)

# load the config file
config <- yaml::read_yaml(opt$config)

# load the config table:
config_df <- read.table(opt$table, sep='\t', header=1)

# output name
outdir <- opt$outdir

# comparison and cluster names
cur_comparison <- opt$comparison
cur_cluster <- opt$cluster
cur_test <- opt$test

# get the options for this test:
test_opt <- config$tests[[cur_test]]$options

# subset the table to get the options:
cur_config <- config_df %>% subset(
  comparison == cur_comparison & 
  cluster == cur_cluster
)

# load the seurat object based on the config file
seurat_obj <- readRDS(config$seurat_path)

# get the group variable:
group_var <- cur_config$groupby

# condition info for the DEG comparison
condition <- cur_config$condition
g1_name <- cur_config$group1 
g2_name <- cur_config$group2 

# replicate info 
replicate_col <- cur_config$replicate  
reps_remove <- strsplit(cur_config$remove_replicates, ';')[[1]]

# additional subsetting:
subset_cols <- cur_config$subset_cols
subset_groups <- cur_config$subset_groups

# additional exclusions:
exclude_cols <- cur_config$exclude_cols
exclude_groups <- cur_config$exclude_groups

# get covariates for the model:
covariates <- config$comparisons[[cur_comparison]]$covariates

# ----------------------------------------------------------------------- #
# Subset the data by cluster of interest & other variables
# ----------------------------------------------------------------------- #

# subset the seurat metadata by cluster
cur_meta <- seurat_obj@meta.data %>%
  subset(
    get(group_var) == cur_cluster & 
    get(condition) %in% c(g1_name, g2_name)
  )

# additional subsetting?
if(!any(subset_groups == "")){
  subset_cols_split <- strsplit(subset_cols, ';')[[1]]
  subset_groups_split <- strsplit(subset_groups, ';')[[1]]
  for(j in 1:length(subset_cols_split)){
    cur_subset <- subset_cols_split[j]
    cur_groups <- strsplit(subset_groups_split[j], ',')[[1]]
    cur_meta <- cur_meta %>% subset(
      get(cur_subset) %in% cur_groups
    )
  }
}

# additional exclusions?
if(!any(exclude_groups == "")){
  exclude_cols_split <- strsplit(exclude_cols, ';')[[1]]
  exclude_groups_split <- strsplit(exclude_groups, ';')[[1]]
  for(j in 1:length(exclude_cols_split)){
    cur_exclude <- exclude_cols_split[j]
    cur_groups <- strsplit(exclude_groups_split[j], ',')[[1]]
    cur_meta <- cur_meta %>% subset(
      !(get(cur_exclude) %in% cur_groups)
    )
  }
}

# filter out replicates that are underrepresented
cur_meta <- cur_meta %>% subset(!(get(replicate_col) %in% reps_remove))

# subset the seurat object based on these barcodes
seurat_test <- seurat_obj[,rownames(cur_meta)]

# ----------------------------------------------------------------------- #
# Check if we have sufficient cells to run DEG
# ----------------------------------------------------------------------- #

min_cells_group <- config$min_cells_group

counts <- table(seurat_test@meta.data[,condition])
counts_g1 <- as.numeric(counts[g1_name]); if(is.na(counts_g1)) counts_g1 <- 0
counts_g2 <- as.numeric(counts[g2_name]); if(is.na(counts_g2)) counts_g2 <- 0

can_run_deg <- (counts_g1 >= min_cells_group) & (counts_g2 >= min_cells_group)
if(!can_run_deg){
    empty_df <- data.frame(
        gene=character(),
        cluster=character(),
        group1=character(),
        group2=character(),
        de_pval=numeric(),
        de_fdr=numeric(),
        log2_de=numeric(),
        pct_group1=numeric(),
        pct_group2=numeric()
    )
        
    # write the output file:
    write.table(
        empty_df, 
        file = paste0(outdir, cur_cluster, "_DEGs.csv"), 
        quote=FALSE, row.names=FALSE, sep=','
    )
    stop()
}
# ----------------------------------------------------------------------- #
# Run the DE test with Seurat FindMarkers
# ----------------------------------------------------------------------- #

# set the Idents:
Idents(seurat_test) <- factor(
  as.character(seurat_test@meta.data[,condition]),
  levels = c(g1_name, g2_name)
)

# run the test
deg_df <- Seurat::FindMarkers(
  seurat_test, 
  ident.1 = g1_name,
  ident.2 = g2_name,
  test.use = test_opt$test_use,
  only.pos = FALSE,
  logfc.threshold = 0,
  latent.vars = covariates,
  slot = test_opt$slot,
  assay = test_opt$assay,
  min.pct = test_opt$min_pct
)

# add relevant columns:
deg_df$gene <- rownames(deg_df)
deg_df$cluster <- cur_cluster
deg_df$group1 <- g1_name
deg_df$group2 <- g2_name

# FDR correction
deg_df$de_fdr <- p.adjust(deg_df$p_val, method='fdr')

# re-name and re-order columns
deg_df <- deg_df %>% 
  dplyr::rename(c(
    de_pval = "p_val",
    log2_de = "avg_log2FC",
    de_pval_adj = "p_val_adj",
    pct_group1 = "pct.1",
    pct_group2 = "pct.2")
  ) %>% 
  dplyr::select(
    c(gene, cluster, group1, group2, de_pval, de_fdr, log2_de, pct_group1, pct_group2)
  )

# sort the results
deg_df <- deg_df %>% dplyr::arrange(desc(log2_de))

# write the output file:
write.table(
    deg_df, 
    file = paste0(outdir, cur_cluster, "_DEGs.csv"), 
    quote=FALSE, row.names=FALSE, sep=','
)
