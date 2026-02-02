library(optparse)


# list of options needed:
# config (config yaml file for snakemake options)
# nthreads (to pass to WGCNA)

option_list <- list(
    optparse::make_option(
        c("--config"), type="character", metavar="character",
        default = 'config.yml',
        help="the configuration .yml file containing all the settings for running the differential analysis snakemake"
    ),
    # optparse::make_option(
    #     c("--table"), type="character", metavar="character",
    #     default = 'config.tsv',
    #     help="the .tsv file containing information about the different comparisons to run'"
    # ),
    # optparse::make_option(
    #     c("--comparison"), type="character", metavar="character",
    #     help="the name of the current comparison to perform the analysis on"
    # ),
    # optparse::make_option(
    #     c("--milo"), type="character", metavar="character",
    #     help="the path to the milo object to run the test on."
    # ),
    # optparse::make_option(
    #     c("--outdir"), type="character", metavar="character",
    #     help="path to output directory."
    # )
)

# parse the arguments:
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)

# --------------------------------------------------------------- #
# load libraries and configuration file
# --------------------------------------------------------------- #

library(yaml)

library(tidyverse)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(MetBrewer)

library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(Matrix)
library(igraph)
library(ggraph)
library(tidygraph)

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = opt$nthreads)
theme_set(theme_cowplot())


# load helper functions
# source("scripts/misc/helper_functions.R")

# load the config file
config <- yaml::read_yaml(opt$config)

# load the config table:
config_df <- read.table(opt$table, sep='\t', header=1)

# output name
outdir <- opt$outdir

# load the seurat object based on the config file
seurat_obj <- readRDS(config$seurat_path)


# continue with the setup 



# --------------------------------------------------------------- #
# Set up the dataset for WGCNA 
# --------------------------------------------------------------- #

# First, subset the Seurat object based on the given parameters
# this should be similar to the MiloR script 

# Feature selection for hdWGCNA 
seurat_obj <- hdWGCNA::SetupForWGCNA(
    seurat_obj,

)

# --------------------------------------------------------------- #
# load libraries and configuration file
# --------------------------------------------------------------- #
