library(optparse)

option_list <- list(
    optparse::make_option(
        c("--config"), type="character", metavar="character",
        default = 'config.yml',
        help="the configuration .yml file containing all the settings for running the differential analysis snakemake"
    )
)

# parse the arguments:
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)

# --------------------------------------------------------------- #
# load libraries and configuration file
# --------------------------------------------------------------- #

library(tidyverse)
library(yaml)
library(Seurat)

# load the config file
config <- yaml::read_yaml(opt$config)

# helper functions
source("scripts/misc/helper_functions.R")

# --------------------------------------------------------------- #
# Set up the configuration table
# --------------------------------------------------------------- #

# load the seurat object
if(!file.exists(config$seurat_path)){
  stop(paste0("Path to Seurat Object specified in config file does not exist. config: ", config_file, ', seurat_path: ', config$seurat_path))
}
seurat_obj <- readRDS(config$seurat_path)

# get the comparisons:
comparisons <- config$comparisons

# loop over each comparison
snake_df <- data.frame()
for(cur_name in names(config$comparisons)){

  cur_comp <- config$comparisons[[cur_name]]

  # cluster column:
  group_var <- cur_comp$groupby
  validate_meta(seurat_obj, group_var)

  # get a list of clusters in the seurat object
  if(is.null(cur_comp$clusters)){
    clusters <- unique(as.character(seurat_obj@meta.data[,group_var]))
  } else{
    clusters <- cur_comp$clusters
    validate_meta(seurat_obj, group_var, clusters)
  }

  # get condition info
  condition <- cur_comp$condition
  g1_name <- cur_comp$group1
  g2_name <- cur_comp$group2
  validate_meta(seurat_obj, condition, c(g1_name, g2_name))

  # replicate info
  replicate_col <- cur_comp$replicate
  validate_meta(seurat_obj, replicate_col)

  # subset options if needed?
  if(!is.null(cur_comp$subset)){

    subset_cols <- paste0(names(cur_comp$subset), collapse=';')

    subset_groups <- list()
    for(cur_subset in names(cur_comp$subset)){
      validate_meta(seurat_obj, cur_subset, cur_comp$subset[[cur_subset]])
      subset_groups[[cur_subset]] <- paste0(cur_comp$subset[[cur_subset]], collapse=',') 
    }
    subset_groups <- paste0(unlist(subset_groups), collapse=';')

  } else{
    subset_cols <- ""
    subset_groups <- ""
  }

  # exclude options if needed?
  if(!is.null(cur_comp$exclude)){

    exclude_cols <- paste0(names(cur_comp$exclude), collapse=';')

    exclude_groups <- list()
    for(cur_exclude in names(cur_comp$exclude)){
      validate_meta(seurat_obj, cur_exclude, cur_comp$exclude[[cur_exclude]])
      exclude_groups[[cur_exclude]] <- paste0(cur_comp$exclude[[cur_exclude]], collapse=',') 
    }
    exclude_groups <- paste0(unlist(exclude_groups), collapse=';')

  } else{
    exclude_cols <- ""
    exclude_groups <- ""
  }

  # finally make the dataframe
  df <- data.frame(
    comparison = cur_name,
    groupby = group_var,
    cluster = clusters,
    replicate = replicate_col,
    condition = condition,
    group1 = g1_name,
    group2 = g2_name,
    subset_cols = subset_cols,
    subset_groups = subset_groups,
    exclude_cols = exclude_cols,
    exclude_groups = exclude_groups,
    remove_replicates = NA
  )

  # filter out groups in cur_df that do not meet criteria:
  rows_keep <- c()
  for(i in 1:nrow(df)){
    
    cur_df <- df[i,]

    # subset the seurat metadata by cluster
    cur_meta <- seurat_obj@meta.data %>%
      subset(
        get(group_var) == cur_df$cluster & 
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
    rep_sizes <- table(cur_meta[,replicate_col])
    reps_remove <- names(which(rep_sizes < config$min_cells_replicate))
    cur_meta <- cur_meta %>% 
      subset(!(get(replicate_col) %in% reps_remove))

    # add these replicates to the table as a semicolon-delimited list 
    df[i,"remove_replicates"] <- paste0(reps_remove, collapse=';')

    # subset by the comparison groups
    if(nrow(cur_meta) >= config$min_cells_group){
      rows_keep <- c(rows_keep, i)
    } 
  }

  df <- df[rows_keep,]
  snake_df <- rbind(snake_df, df)

}

# save the file!
write.table(
  snake_df, 
  file = 'config.tsv',
  sep = '\t',
  row.names=FALSE, quote=FALSE
)