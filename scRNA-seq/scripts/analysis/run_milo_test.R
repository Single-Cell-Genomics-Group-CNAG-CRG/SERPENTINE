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
        c("--milo"), type="character", metavar="character",
        help="the path to the milo object to run the test on."
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
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)

# load helper functions
source("scripts/misc/helper_functions.R")

# load the config file
config <- yaml::read_yaml(opt$config)

# load the config table:
config_df <- read.table(opt$table, sep='\t', header=1)

# output name
outdir <- opt$outdir

# comparison and cluster names
cur_comparison <- opt$comparison

# subset the table to get the options:
cur_config <- config_df %>% subset(
  comparison == cur_comparison 
)

# load the milo object based on the config file
milo_obj <- readRDS(opt$milo)

# get the group variable:
group_var <- unique(cur_config$groupby)

# condition info for the DEG comparison
condition <- cur_config$condition[1]
g1_name <- cur_config$group1[1]
g2_name <- cur_config$group2[1] 

# replicate info 
replicate_col <- unique(cur_config$replicate)

# settings for milo setup 
milo_setup <- config$comparisons[[cur_comparison]]$milo_setup 
k <- milo_setup$k 
prop <- milo_setup$prop
reduction <- milo_setup$reduction 
dims <- milo_setup$dims
mixed_prop <- milo_setup$mixed_prop

# get the list of covariates
covariates <- config$comparisons[[cur_comparison]]$covariates

# TODO: CHECK THAT COVARIATES ARE PRESENT in the dataset

# convert the reduction to uppercase
reduction_milo <- toupper(reduction)

# additional subsetting:
subset_cols <- cur_config$subset_cols[1]
subset_groups <- cur_config$subset_groups[1]

# additional exclusions:
exclude_cols <- cur_config$exclude_cols[1]
exclude_groups <- cur_config$exclude_groups[1]

# what columns do we need to subset?
if(!any(subset_groups == "")){
  subset_cols_split <- strsplit(subset_cols, ';')[[1]]
  subset_groups_split <- strsplit(subset_groups, ';')[[1]]
} else{
  subset_cols_split <- NULL; subset_groups <- NULL; subset_groups_split <- NULL
}

# what columns do we need to exclude?
if(!any(exclude_groups == "")){
  exclude_cols_split <- strsplit(exclude_cols, ';')[[1]]
  exclude_groups_split <- strsplit(exclude_groups, ';')[[1]]
} else{
  exclude_cols_split <- NULL; exclude_groups <- NULL; exclude_groups_split <- NULL
}

#------------------------------------------------------
# Create the design matrix
#------------------------------------------------------

meta_cols <- c(condition, subset_cols_split, exclude_cols_split, covariates)

# initialize the design matrix:
design_df <- as.data.frame(colData(milo_obj))[,c(replicate_col, meta_cols)]

# subset the design matrix if necessary
if(!is.null(subset_groups)){
  for(j in 1:length(subset_cols_split)){
    cur_subset <- subset_cols_split[j]
    cur_groups <- strsplit(subset_groups_split[j], ',')[[1]]
    design_df <- design_df %>% subset(
      as.character(get(cur_subset)) %in% cur_groups
    )
  }
}

# exclude things from the design matrix if necessary
if(!is.null(exclude_groups)){
  for(j in 1:length(exclude_cols_split)){
    cur_exclude <- exclude_cols_split[j]
    cur_groups <- strsplit(exclude_groups_split[j], ',')[[1]]
    design_df <- design_df %>% subset(
      !(as.character(get(cur_exclude)) %in% cur_groups)
    )
  }
}

design_df <- fix_metadata(
  design_df,
  cols = meta_cols
)

g1_name_fix <- fix_metadata(g1_name)
g2_name_fix <- fix_metadata(g2_name)

# create the design matrix 
design_df <- design_df[,c(replicate_col, paste0(meta_cols, '_fix'))] %>% 
  as.data.frame() %>% 
  dplyr::distinct()

# re-name the row and column names
colnames(design_df) <- c(replicate_col, meta_cols)
rownames(design_df) <- design_df[,replicate_col]

# set up the contrast of interest!
cur_contrast <- paste0(paste0(condition, g1_name_fix), '-', paste0(condition, g2_name_fix))

# define the formula
formula <- as.formula(paste0(c("~ 0", covariates, condition), collapse=" + "))

#------------------------------------------------------
# Run the differential abundance test!
#------------------------------------------------------

da_results <- miloR::testNhoods(
  milo_obj, 
  design = formula, 
  design.df = design_df, 
  model.contrasts = cur_contrast,
  fdr.weighting="graph-overlap", 
  norm.method="TMM",
  reduced.dim = reduction_milo
)

# add neighborhood annotations
da_results <- miloR::annotateNhoods(
  milo_obj, 
  da_results, 
  coldata_col = group_var
)

# determine which neighborhoods are "Mixed" with respect to the annotation based on a threshold
da_results[,group_var] <- ifelse(da_results[,paste0(group_var, '_fraction')] <= mixed_prop, "Mixed", da_results[,group_var])

# re-name:
colnames(da_results)[(ncol(da_results)-1):ncol(da_results)] <- c('cluster', 'cluster_fraction')

# add information about the specific DA test
da_results$group1 <- g1_name 
da_results$group2 <- g2_name 

# re-order columns:
da_results <- da_results %>% 
    dplyr::select(
        c(
            Nhood, cluster, group1, group2, PValue, 
            FDR, SpatialFDR, logFC, logCPM, F, cluster_fraction
        )
    )

# save the DA results table!
write.table(
  da_results,
  file = paste0(outdir, cur_comparison, '_milo_DA.csv'),
  quote=FALSE, row.names=FALSE, sep=','
)