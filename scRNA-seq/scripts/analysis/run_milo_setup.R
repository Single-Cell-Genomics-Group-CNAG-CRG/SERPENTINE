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

# load the seurat object based on the config file
seurat_obj <- readRDS(config$seurat_path)

# get the group variable:
group_var <- unique(cur_config$groupby)

# replicate info 
replicate_col <- unique(cur_config$replicate)

# settings for milo setup 
milo_setup <- config$comparisons[[cur_comparison]]$milo_setup 
k <- milo_setup$k 
prop <- milo_setup$prop
reduction <- milo_setup$reduction 
dims <- milo_setup$dims

# conver the reduction to uppercase
reduction_milo <- toupper(reduction)

# additional subsetting:
subset_cols <- cur_config$subset_cols[1]
subset_groups <- cur_config$subset_groups[1]

#------------------------------------------------------
# Create the milo object
#------------------------------------------------------

# remove outlier replicates
rep_counts <- table(seurat_obj@meta.data[,replicate_col])
reps_remove <- names(which(rep_counts < config$min_cells_replicate))
cur_meta <- seurat_obj@meta.data %>% subset(
  !(get(replicate_col) %in% reps_remove)
)

# additional subsetting?
if(milo_setup$subset){
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
}

# subset the seurat object based on these barcodes
seurat_obj <- seurat_obj[,rownames(cur_meta)]

# convert from SCE to milo object
milo_obj <- miloR::Milo(as.SingleCellExperiment(seurat_obj))

# remove the seurat object to save memory usage 
rm(seurat_obj); gc();

# compute knn
milo_obj <- miloR::buildGraph(
  milo_obj, 
  k = k, 
  d = dims, 
  reduced.dim = reduction_milo
)

# defining neighborhoods
milo_obj <- miloR::makeNhoods(
  milo_obj, 
  prop = prop, 
  k = k, 
  d = dims, 
  refined = TRUE, 
  reduced_dims = reduction_milo
)

# count cells in each neighborhood 
milo_obj <- miloR::countCells(
    milo_obj, 
    meta.data = as.data.frame(colData(milo_obj)), 
    sample=replicate_col
)

# graph of neighborhoods
milo_obj <- miloR::buildNhoodGraph(milo_obj)

# save the object!
milo_path <- paste0(outdir, cur_comparison, "_milo.rds")
saveRDS(milo_obj, milo_path)
