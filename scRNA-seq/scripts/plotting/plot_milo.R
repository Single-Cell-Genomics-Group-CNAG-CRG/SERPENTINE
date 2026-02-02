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
        c("--csv"), type="character", metavar="character",
        help="path to a .csv file containing the DEG results for a single cluster."
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
library(patchwork)
library(cowplot)
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)

# set the general plot theme
theme_set(theme_cowplot())

# load helper functions
source("scripts/misc/helper_functions.R")

# load the plotting functions
source("scripts/plotting/plotting_functions.R")

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

# load the milo test results:
da_results <- read.csv(opt$csv)

# get the group variable:
group_var <- unique(cur_config$groupby)

# condition info for the DEG comparison
condition <- cur_config$condition[1]
g1_name <- as.character(cur_config$group1[1])
g2_name <- as.character(cur_config$group2[1])

# replicate info 
replicate_col <- unique(cur_config$replicate)
reps_remove <- strsplit(cur_config$remove_replicates, ';')[[1]]

# additional subsetting:
subset_cols <- cur_config$subset_cols[1]
subset_groups <- cur_config$subset_groups[1]

# additional exclusions:
exclude_cols <- cur_config$exclude_cols[1]
exclude_groups <- cur_config$exclude_groups[1]

# settings for milo setup 
milo_setup <- config$comparisons[[cur_comparison]]$milo_setup 
k <- milo_setup$k 
prop <- milo_setup$prop
reduction <- milo_setup$reduction 
dims <- milo_setup$dims
mixed_prop <- milo_setup$mixed_prop

# get the list of covariates
covariates <- config$comparisons[[cur_comparison]]$covariates

#------------------------------------------------------
# Plot the proportions in each group as a barplot
#------------------------------------------------------

pl_opt <- config$proportion_barplot 
plot_width <- pl_opt$width
plot_height <- pl_opt$height

p <- ProportionBarPlot(
  object = milo_obj,
  group_var = group_var,
  condition = condition,
  g1_name = g1_name,
  g2_name = g2_name,
  replicate_col = replicate_col,
  reps_remove = reps_remove,
  subset_groups = subset_groups,
  subset_cols = subset_cols,
  exclude_groups = exclude_groups,
  exclude_cols = exclude_cols,
  order_groups = TRUE
) +  ggtitle(paste0(g1_name, ' vs. ', g2_name))


pdf(paste0(outdir, cur_comparison, '_milo_proportions.pdf'), width=plot_width, height=plot_height)
print(p)
dev.off()

#------------------------------------------------------
# Plot the raw number of cells per group
#------------------------------------------------------

pl_opt <- config$cell_count_barplot 
plot_width <- pl_opt$width
plot_height <- pl_opt$height

p <- ProportionBarPlot(
  object = milo_obj,
  group_var = group_var,
  condition = condition,
  g1_name = g1_name,
  g2_name = g2_name,
  replicate_col = replicate_col,
  reps_remove = reps_remove,
  subset_groups = subset_groups,
  subset_cols = subset_cols,
  exclude_groups = exclude_groups,
  exclude_cols = exclude_cols,
  plot_raw_counts = TRUE
) +  ggtitle(paste0(g1_name, ' vs. ', g2_name))

pdf(paste0(outdir, cur_comparison, '_cell_counts.pdf'), width=plot_width, height=plot_height)
print(p)
dev.off()

#------------------------------------------------------
# Plot the effect sizes as a box and whisker
#------------------------------------------------------

pl_opt <- config$milo_boxplot 
plot_width <- pl_opt$width
plot_height <- pl_opt$height
signif_col <- pl_opt$signif_col

print(head(da_results))
print('pl_opt:')
print(pl_opt)

p <- DABoxPlot(
  da_results,
  signif_col = signif_col
) + 
  ggtitle(paste0(g1_name, ' vs. ', g2_name))

pdf(paste0(outdir, cur_comparison, '_milo_distributions.pdf'), width=plot_width, height=plot_height)
print(p)
dev.off()

#------------------------------------------------------
# Plot the effect sizes on the dim reduction
#------------------------------------------------------

pl_opt <- config$milo_reduction
plot_width <- pl_opt$width
plot_height <- pl_opt$height
signif_col <- pl_opt$signif_col
reduction_plot <- toupper(pl_opt$reduction)

p <- DAReductionPlot(
    milo_obj,
    da_results,
    centroids = group_var,
    reduction = reduction_plot,
    signif_col = signif_col,
    insignif_color = "lightgrey"
) + 
  ggtitle(paste0(g1_name, ' vs. ', g2_name)) + 
  coord_fixed()

pdf(paste0(outdir, cur_comparison, '_milo_reduction.pdf'), width=plot_width, height=plot_height)
print(p)
dev.off()