#' Validate a Metadata Column in a Seurat Object
#'
#' This function checks whether a specified metadata column exists in a Seurat object,
#' ensures that its values do not contain forbidden characters ("," or ";"),
#' and verifies that specific values (if provided) exist in the column.
#' If any of these checks fail, the function stops execution with an error message.
#'
#' @param seurat_obj A Seurat object containing single-cell metadata.
#' @param col A character string specifying the metadata column to validate.
#' @param values A character vector of expected values to check within the column. Defaults to an empty vector.
#'
#' @return NULL (Stops execution if validation fails).
#'
validate_meta <- function(
  seurat_obj,
  col,
  values = c()
){
  
  meta <- seurat_obj@meta.data
  if(!(col %in% colnames(meta))){
    if(file.exists('config.tsv')){unlink('config.tsv')}
    stop(paste0("Column ", col, " not found in colnames(seurat_obj@meta.data)."))
  }

  # check for commas or semicolons in this column:
  if(any(grepl("[,;]", seurat_obj@meta.data[,col]))){
    stop(paste0("Values in column ", col, " contain , or ; characters. Please update this column to remove these characters and try running again."))
  }

  # check that the selected values are valid in this column:
  valid <- unique(as.character(meta[,col]))
  for(value in values){
    if(!(value %in% valid)){
      if(file.exists('config.tsv')){unlink('config.tsv')}
      stop("Value ", value, " not found in seurat_obj@meta.data$", col, ".")
    }
  }
}

#' Run Enrichr with Automatic Retries on Failure
#'
#' This function performs enrichment analysis using Enrichr and automatically
#' retries in case of connection failures. If Enrichr API requests fail, the
#' function waits and re-attempts up to `max_attempts` times.
#'
#' @param gene_list A character vector of gene symbols for enrichment analysis.
#' @param dbs A character vector specifying the Enrichr databases to query.
#' @param max_attempts The maximum number of times to retry in case of failure (default: 5).
#' @param wait_time Time (in seconds) to wait between retry attempts (default: 5).
#'
#' @return A list of enrichment results, where each element corresponds to a database
#'         from `dbs`. If all attempts fail, returns `NULL`.
run_enrichr_safe <- function(gene_list, dbs, max_attempts = 5, wait_time = 5) {
    attempt <- 1
    while(attempt <= max_attempts) {
        try({
            enrich_results <- enrichR::enrichr(gene_list, dbs)
            return(enrich_results)  # Success
        }, silent = TRUE)
        
        # If it failed, wait and retry
        message(paste("Attempt", attempt, "failed. Retrying in", wait_time, "seconds..."))
        Sys.sleep(wait_time)
        attempt <- attempt + 1
    }
    
    message("Enrichr request failed after multiple attempts.")
    return(NULL)  # Return NULL if all retries fail
}


#' Perform enrichment analysis using Enrichr and return filtered results.
#'
#' This function takes a list of gene sets (`input_list`) and runs enrichment
#' analysis using the provided databases (`dbs`). It combines results across
#' clusters and filters out entries with fewer than 3 genes.
#'
#' @param input_list A named list where each element is a vector of gene symbols.
#' @param dbs A character vector of database names for Enrichr analysis.
#'
#' @return A data frame containing enrichment results, filtered for gene sets
#'         with at least 3 genes.
TestEnrichment <- function(input_list, dbs) {

  # ------------------------------------------------------------ #
  # Perform Input Checks
  # ------------------------------------------------------------ #
  
  # Ensure input_list is a named list
  if (!is.list(input_list) || is.null(names(input_list))) {
    stop("Error: 'input_list' must be a named list of gene sets.")
  }
  
  # Ensure dbs is a character vector
  if (!is.character(dbs) || length(dbs) == 0) {
    stop("Error: 'dbs' must be a non-empty character vector of database names.")
  }
  
  # Ensure input_list is not empty
  if (length(input_list) == 0) {
    warning("Warning: 'input_list' is empty. No enrichment will be performed.")
    return(data.frame())  # Return an empty dataframe
  }
  
  # ------------------------------------------------------------ #
  # Run Enrichr and Combine Outputs
  # ------------------------------------------------------------ #
  
  enrich_df <- do.call(rbind, lapply(names(input_list), function(x) {
    
    # Skip empty gene sets
    if (length(input_list[[x]]) == 0) return(data.frame())
    
    # Run enrichment safely
    cur_enrich <- run_enrichr_safe(input_list[[x]], dbs)
    
    # Check if enrichment returned results
    if (is.null(cur_enrich)) {
      warning(paste("Warning: Enrichment failed for group", x, "Skipping..."))
      return(data.frame()) 
    }
    
    # Process results for each database
    cur_df <- do.call(rbind, lapply(dbs, function(cur_db) {
      if (!is.null(cur_enrich[[cur_db]]) && nrow(cur_enrich[[cur_db]]) > 1) {
        df <- cur_enrich[[cur_db]]
        df$group <- x
        df$db <- cur_db
      } else {
        df <- data.frame()
      }
      df
    }))
    
    # Only keep entries with 3 or more genes
    if (nrow(cur_df) > 0 && "Genes" %in% colnames(cur_df)) {
      cur_df$ngenes <- sapply(strsplit(cur_df$Genes, ';'), length)
      cur_df <- subset(cur_df, ngenes >= 3)
    }
    
    return(cur_df)
  }))
  
  # Return the combined enrichment table
  return(enrich_df)
}


#' fix_metadata: Clean and Modify Metadata Column Names
#'
#' This function processes specified columns in a data frame by removing or replacing 
#' unwanted characters using regular expressions and appends a suffix to the modified column names.
#'
#' @param df A data frame containing the metadata.
#' @param cols A character vector of column names to be modified.
#' @param pattern_remove A regular expression pattern specifying characters to remove. 
#'   Default is `"[ ()-]"` (spaces, parentheses, and hyphens).
#' @param pattern_sub A string to replace matched characters with. Default is `""` (empty string).
#' @param append A string to append to modified column names. Default is `"_fix"`.
#'
#' @return A data frame with new columns containing cleaned metadata.
#'
#' @export
fix_metadata <- function(
  input, 
  cols=NULL,
  pattern_remove = "[ ()-]",
  pattern_sub = "",
  append = "_fix"
){
  if("data.frame" %in% class(input)){
    for(cur_col in cols){
      input[,paste0(cur_col, append)] <- gsub(pattern_remove, pattern_sub, input[,cur_col])
    }

  } else if("character" %in% class(input)){
    input <- gsub(pattern_remove, pattern_sub, input)
  }

  input
}
