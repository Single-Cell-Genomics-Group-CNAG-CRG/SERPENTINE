library(optparse)

option_list <- list(
    optparse::make_option(
        c("--config"), type="character", metavar="character",
        default = 'config.yml',
        help="the configuration .yml file containing all the settings for running the differential analysis snakemake"
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

# set the general plot theme
theme_set(theme_cowplot())

# load the plotting functions
source("scripts/plotting/plotting_functions.R")

# load the config file
config <- yaml::read_yaml(opt$config)

# load the DEG results:
deg_df <- read.csv(opt$csv)
g1_name <- unique(deg_df$group1)
g2_name <- unique(deg_df$group2)

# output name
outdir <- opt$outdir

# get the cluster name from the deg_df
cur_cluster <- unique(deg_df$cluster)

# plotting options
pl_opt <- config$volcano
plot_width <- pl_opt$width
plot_height <- pl_opt$height

# load the DEG results
# deg_df <- read.csv(paste0(config$outdir, 'cluster_', cur_cluster, '_DEGs.csv'))

# --------------------------------------------------------------- #
# Plot the DE volcano plot
# --------------------------------------------------------------- #

plot_title <- paste0(cur_cluster, ' ', g1_name, ' vs. ', g2_name)

# make the DE volcano plot
p <- VolcanoPlot(
    deg_df,
    signif_col = pl_opt$de_signif,
    effect_col = pl_opt$de_effect,
    effect_thresh = pl_opt$effect_thresh,
    x_label = bquote("log"[2]~"(Fold Change), Mean")
) + 
    ggtitle(plot_title) + 
    theme(plot.title = element_text(size=pl_opt$title_size))

# define the output file name
plot_outfile <- paste0(outdir, cur_cluster, '_DE_volcano.pdf')

# write the output file
pdf(plot_outfile, width=plot_width, height=plot_height)
print(p)
dev.off()

# --------------------------------------------------------------- #
# Plot the DV volcano plot (memento only)
# --------------------------------------------------------------- #

if(any(grepl("dv", colnames(deg_df)))){

    # make the DV volcano plot
    p <- VolcanoPlot(
        deg_df,
        signif_col = pl_opt$dv_signif,
        effect_col = pl_opt$dv_effect,
        effect_thresh = pl_opt$effect_thresh,
        x_label = bquote("log"[2]~"(Fold Change), Variance"),
        y_label = bquote("-log"[10]~"(P-value)")
    ) + 
        ggtitle(plot_title) + 
        theme(plot.title = element_text(size=pl_opt$title_size))


    # define the output file name
    plot_outfile <- paste0(outdir, cur_cluster, '_DV_volcano.pdf')

    # write the output file
    pdf(plot_outfile, width=plot_width, height=plot_height)
    print(p)
    dev.off()

    # make the DE / DV comparison plot
    p <- EffectComparisonPlot(
        deg_df,
        signif_col1 = pl_opt$de_signif,
        signif_col2 = pl_opt$dv_signif,
        effect_col1 = pl_opt$de_effect,
        effect_col2 = pl_opt$dv_effect
    ) + 
        ggtitle(plot_title) + 
        theme(plot.title = element_text(size=pl_opt$title_size))

    # define the output file name
    plot_outfile <- paste0(outdir, cur_cluster, '_DE_DV_comparison.pdf')

    # write the output file
    pdf(plot_outfile, width=plot_width, height=plot_height)
    print(p)
    dev.off()

}
