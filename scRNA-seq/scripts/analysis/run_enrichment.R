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
library(enrichR)
library(simplifyEnrichment)

# set the general plot theme
theme_set(theme_cowplot())

# load the plotting functions
source("scripts/plotting/plotting_functions.R")
source("scripts/misc/helper_functions.R")

# load the config file
config <- yaml::read_yaml(opt$config)

# load the DEG results:
deg_df <- read.csv(opt$csv)
g1_name <- as.character(unique(deg_df$group1))
g2_name <- as.character(unique(deg_df$group2))

# output name
outdir <- opt$outdir

# get the cluster name from the deg_df
cur_cluster <- unique(deg_df$cluster)

# get options for filtering the DEG table from the config file
pl_opt <- config$volcano
plot_width <- config$enrichment$width
plot_height <- config$enrichment$height

# --------------------------------------------------------------- #
# Set up the gene lists for enrichR
# --------------------------------------------------------------- #

# get the genes significantly DE for this cluster
up <- deg_df %>% subset(
    get(pl_opt$de_signif) <= 0.05 & 
    get(pl_opt$de_effect) >= pl_opt$effect_thresh
)

down <- deg_df %>% subset(
    get(pl_opt$de_signif) <= 0.05 & 
    get(pl_opt$de_effect) <= -pl_opt$effect_thresh
)

input_list <- list()
input_list[[g1_name]] <- up$gene
input_list[[g2_name]] <- down$gene


# If there's Not enough DEGs, we can't run the enrichment test 
# TODO

# --------------------------------------------------------------- #
# Run enrichR for DEGs in each DB
# --------------------------------------------------------------- #

# get the list of dbs from the config file:
dbs <- config$enrichment$dbs

# use the helper function to run enrichR safely
enrich_df <- TestEnrichment(input_list, dbs)

# --------------------------------------------------------------- #
# Write the output file
# --------------------------------------------------------------- #

write.table(
    enrich_df, 
    file = paste0(outdir, "/enrichment/", cur_cluster, "_enrichR.tsv"), 
    quote=FALSE, row.names=FALSE, sep='\t'
)

# --------------------------------------------------------------- #
# plot the results with simplifyEnrichment
# --------------------------------------------------------------- #

# only plot the results from Biological Process go term:
if(any(grepl('GO_Biological', dbs))){

  cur_db <- dbs[grepl('GO_Biological', dbs)][1]

  # make the Enrichment plot for group1:
  cur_terms <- subset(enrich_df, group == g1_name & P.value <= 0.05 & db == cur_db)
  go_ids <- gsub(".*\\(GO:(\\d+)\\)", "GO:\\1", cur_terms$Term)

  plot_out <- paste0(outdir, "/figures/", cur_cluster, "_group1_simplifyEnrichment.pdf")

  # if there aren't enough GO terms (fewer than 5), make an empty plot
  if(length(go_ids) < 5){
    pdf(plot_out, height=plot_height, width=plot_width)
    dev.off()
  } else{

    # run semantic similarity analysis with simplifyEnrichment
    mat <- simplifyEnrichment::GO_similarity(go_ids, ont='BP')

    pdf(plot_out, height=plot_height, width=plot_width)
    simplifyEnrichment::simplifyGO(mat)
    dev.off()
  }

  # make the Enrichment plot for group2:
  cur_terms <- subset(enrich_df, group == g2_name & P.value <= 0.05 & db == cur_db)
  go_ids <- gsub(".*\\(GO:(\\d+)\\)", "GO:\\1", cur_terms$Term)

  plot_out <- paste0(outdir, "/figures/", cur_cluster, "_group2_simplifyEnrichment.pdf")

  # if there aren't enough GO terms (fewer than 5), make an empty plot
  if(length(go_ids) < 5){
    pdf(plot_out, height=plot_height, width=plot_width)
    dev.off()
  } else{
    # run semantic similarity analysis with simplifyEnrichment
    mat <- simplifyEnrichment::GO_similarity(go_ids, ont='BP')

    pdf(plot_out, height=plot_height, width=plot_width)
    simplifyEnrichment::simplifyGO(mat)
    dev.off()
  }

}