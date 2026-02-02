library(optparse)

option_list <- list(
    optparse::make_option(
        c("--config"), type="character", metavar="character",
        default = 'config.yml',
        help="the configuration .yml file containing all the settings for running the differential analysis snakemake"
    ),
    optparse::make_option(
        c("--comparison"), type="character", metavar="character",
        help="the name of the current comparison to perform the analysis on"
    ),
    optparse::make_option(
        c("--results"), type="character", metavar="character",
        help="semicolon-delimited list of files to concatenate"
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

# load the config file
config <- yaml::read_yaml(opt$config)

# output name
outdir <- opt$outdir
comparison <- opt$comparison
test <- opt$test

# options for subsetting the data for significant results
pl_opt <- config$volcano

# get the list of results files
DEG_tests <- strsplit(opt$results, ';')[[1]]

# read the DEG results from each cluster, and combine into one table
combined <- Reduce(rbind, lapply(DEG_tests, function(file){
  read.csv(file)
}))

# write the combined file as a new csv:
write.csv(
    combined,
    file = paste0(outdir, '/', comparison, '_', test, '_DEGs.csv'),
    row.names=FALSE,
    quote=FALSE
)

# ----------------------------------------------------------------------- #
# subset for significant DEGs (use volcano options)
# ----------------------------------------------------------------------- #

signif <- combined %>% subset(
    get(pl_opt$de_signif) <= 0.05 &
    abs(get(pl_opt$de_effect)) >= pl_opt$effect_thresh
)

# write the signif file as a new csv:
write.csv(
    signif,
    file = paste0(outdir, '/', comparison, '_', test, '_DEGs_signif.csv'),
    row.names=FALSE,
    quote=FALSE
)