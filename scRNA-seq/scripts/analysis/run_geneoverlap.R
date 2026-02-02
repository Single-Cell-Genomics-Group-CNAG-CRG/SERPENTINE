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
library(GeneOverlap)
library(fgsea)
library(simplifyEnrichment)
library(Seurat)

# set the general plot theme
theme_set(theme_cowplot())

# load the plotting functions
source("scripts/plotting/plotting_functions.R")
source("scripts/misc/helper_functions.R")

# load the config file
config <- yaml::read_yaml(opt$config)

# load the seurat object
if(!file.exists(config$seurat_path)){
  stop(paste0("Path to Seurat Object specified in config file does not exist. config: ", config_file, ', seurat_path: ', config$seurat_path))
}
seurat_obj <- readRDS(config$seurat_path)

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

# get the list of db files from the config file:
dbs <- config$enrichment$db_files


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

# --------------------------------------------------------------- #
# Calculate overlap between DEGs and gene lists
# --------------------------------------------------------------- #


# get the list of dbs
dbs <- config$enrichment$db_files
cur_db <- dbs[1]

overlap_df <- data.frame()
for(cur_db in dbs){

  db_name <- str_split(cur_db, '/')[[1]]
  db_name <- tools::file_path_sans_ext(db_name[length(db_name)])

  print(db_name)

  # load the GO Biological Pathways file (downloaded from EnrichR website)
  pathways <- fgsea::gmtPathways(cur_db)

  # run the overlap test
  gom.obj <- GeneOverlap::newGOM(
    pathways, 
    input_list,
    genome.size = nrow(seurat_obj)
  )

  # need to convert this into a table!!!
  cur_overlap_df <- data.frame() 
  for(i in 1:length(gom.obj@go.nested.list)){
    for(j in 1:length(gom.obj@go.nested.list[[i]])){

      cur <- gom.obj@go.nested.list[[i]][[j]]
      cur_df <-  data.frame(
        term = names(pathways)[j],
        overlap = paste0(length(cur@intersection), '|', length(cur@listA)),
        genes = paste0(cur@intersection, collapse=','),
        pval = cur@pval,
        odds_ratio = cur@odds.ratio,
        jaccard = cur@Jaccard,
        db = db_name,
        group = names(input_list)[i],
        cluster = cur_cluster,
        ngenes = length(cur@intersection)
      )
      cur_overlap_df <- rbind(cur_overlap_df, cur_df)
    }
  }

  # remove empty entries
  cur_overlap_df <- subset(cur_overlap_df, ngenes > 0)

  # calculate adjusted p-val
  cur_overlap_df$fdr <- p.adjust(cur_overlap_df$pval, 'fdr')

  # order by significance and effect size:
  cur_overlap_df <- rbind(
    subset(cur_overlap_df, fdr <= 0.05) %>% arrange(desc(odds_ratio)),
    subset(cur_overlap_df, fdr > 0.05) %>% arrange(desc(odds_ratio))
  )

  overlap_df <- rbind(overlap_df, cur_overlap_df)

}

# --------------------------------------------------------------- #
# Write the output file
# --------------------------------------------------------------- #

write.table(
    overlap_df, 
    file = paste0(outdir, "/enrichment/", cur_cluster, "_enrichment.tsv"), 
    quote=FALSE, row.names=FALSE, sep='\t'
)

signif_df <- subset(overlap_df, fdr <= 0.05)
write.table(
    signif_df, 
    file = paste0(outdir, "/enrichment/", cur_cluster, "_enrichment_signif.tsv"), 
    quote=FALSE, row.names=FALSE, sep='\t'
)

# --------------------------------------------------------------- #
# plot the results with simplifyEnrichment
# --------------------------------------------------------------- #

dbs <- unique(overlap_df$db)

# only plot the results from Biological Process go term:
if(any(grepl('GO_Biological', dbs))){

  cur_db <- dbs[grepl('GO_Biological', dbs)][1]

  # make the Enrichment plot for group1:
  cur_terms <- subset(signif_df, group == g1_name & db == cur_db)
  go_ids <- gsub(".*\\(GO:(\\d+)\\)", "GO:\\1", cur_terms$term)

  plot_out <- paste0(outdir, "/figures/", cur_cluster, "_", g1_name, "_simplifyEnrichment.pdf")

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
  cur_terms <- subset(signif_df, group == g2_name & db == cur_db)
  go_ids <- gsub(".*\\(GO:(\\d+)\\)", "GO:\\1", cur_terms$term)

  plot_out <- paste0(outdir, "/figures/", cur_cluster, "_", g2_name, "_simplifyEnrichment.pdf")

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