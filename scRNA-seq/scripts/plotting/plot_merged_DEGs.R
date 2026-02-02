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
        c("--csv"), type="character", metavar="character",
        help="path to a .csv file containing the merged DEG results across all clusters."
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
library(patchwork)
library(cowplot)

# set the general plot theme
theme_set(theme_cowplot())

# source the plotting functions
source("scripts/plotting/plotting_functions.R")

# load the config file
config <- yaml::read_yaml(opt$config)

# load the DEG results:
deg_df <- read.csv(opt$csv)
g1_name <- unique(deg_df$group1)
g2_name <- unique(deg_df$group2)

# output name
outdir <- opt$outdir
comparison <- opt$comparison
test <- opt$test

pl_opt <- config$multivolcano
plot_width <- pl_opt$width
plot_height <- pl_opt$height

# --------------------------------------------------------------- #
# Plot the DE multi volcano plot
# --------------------------------------------------------------- #

# determine the plot height 
n_clusters <- length(unique(deg_df$cluster))
plot_h <- max(
  3,
  plot_height * (n_clusters/2) 
)

plot_title <- paste0(g1_name, ' vs. ', g2_name)

p <- MultiVolcanoPlot(
    deg_df,
    signif_col = pl_opt$de_signif,
    effect_col = pl_opt$de_effect,
    effect_thresh = pl_opt$effect_thresh,
    x_label = bquote("log"[2]~"(Fold Change), Mean"),
) + 
    ggtitle(plot_title) + 
    theme(plot.title = element_text(size=pl_opt$title_size))

# TODO
# should I dynamically change the plot size based on the number of clusters?
plot_outfile <- paste0(outdir, '/', comparison, '_', test, '_DE_multivolcano.pdf')

pdf(plot_outfile, width=plot_width, height=plot_h)
print(p)
dev.off()

# --------------------------------------------------------------- #
# Plot the DV multi volcano plot (memento only)
# --------------------------------------------------------------- #

if(any(grepl("dv", colnames(deg_df)))){

    p <- MultiVolcanoPlot(
        deg_df,
        signif_col = pl_opt$dv_signif,
        effect_col = pl_opt$dv_effect,
        effect_thresh = pl_opt$effect_thresh,
        x_label = bquote("log"[2]~"(Fold Change), Variance"),
    ) + 
        ggtitle(plot_title) + 
        theme(plot.title = element_text(size=pl_opt$title_size))

    # should I dynamically change the plot size based on the number of clusters?
    plot_outfile <- paste0(outdir, '/', comparison, '_', test, '_DV_multivolcano.pdf')

    pdf(plot_outfile, width=plot_width, height=plot_h)
    print(p)
    dev.off()

}