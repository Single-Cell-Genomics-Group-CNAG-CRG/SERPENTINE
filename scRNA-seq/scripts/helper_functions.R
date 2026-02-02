library(yaml)
library(Seurat)


#-------------------------------------------------------------
# Helper: validate metadata columns and optional values
#-------------------------------------------------------------
validate_metadata_column <- function(seu, column, values = NULL, context = "",
                                     verbose = TRUE, errors = list()) {
  vcat <- function(...) if (verbose) cat(..., "\n")
  prefix <- ifelse(context != "", paste0("[", context, "] "), "")
  
  # 1. Column existence
  if (!(column %in% colnames(seu@meta.data))) {
    msg <- sprintf("❌ %sMetadata column '%s' not found.", prefix, column)
    errors <- c(errors, msg)
    return(errors)
  } else {
    vcat(green(paste0("✅ ", prefix, "Metadata column '", column, "' found.")))
  }
  
  # 2. Check values (if provided)
  if (!is.null(values)) {
    meta_vals <- unique(seu@meta.data[[column]])
    missing_vals <- setdiff(values, meta_vals)
    if (length(missing_vals) > 0) {
      msg <- sprintf("❌ %sColumn '%s' missing expected values: %s",
                     prefix, column, paste(missing_vals, collapse = ", "))
      errors <- c(errors, msg)
    } else {
      vcat(green(paste0("✅ ", prefix, "All specified values found in '", column, "'.")))
    }
  }
  return(errors)
}

#-------------------------------------------------------------
# Main: load_and_validate_config
#-------------------------------------------------------------
load_and_validate_config <- function(config_file, wgcna_name, verbose = TRUE) {
  vcat <- function(...) if (verbose) cat(..., "\n")
  errors <- list()
  
  vcat(blue$bold("🔍 Reading configuration file:"), config_file)
  
  # 1. Read config file
  if (!file.exists(config_file)) {
    stop(red$bold("❌ Config file not found: "), config_file)
  }
  config <- tryCatch(
    yaml::read_yaml(config_file),
    error = function(e) stop(red$bold("❌ Failed to parse YAML: "), e$message)
  )
  
  # 2. Check Seurat path
  if (is.null(config$seurat_path)) {
    errors <- c(errors, "❌ Missing 'seurat_path' in config.")
  } else if (!file.exists(config$seurat_path)) {
    errors <- c(errors, paste0("❌ Seurat object not found at: ", config$seurat_path))
  }
  
  # 3. Validate general section
  required_general <- c("out_dir", "fig_dir", "data_dir", "cores")
  if (is.null(config$general)) {
    errors <- c(errors, "❌ Missing 'general' section in config.")
  } else {
    missing_general <- setdiff(required_general, names(config$general))
    if (length(missing_general) > 0) {
      errors <- c(errors, paste0("❌ Missing general fields: ", paste(missing_general, collapse = ", ")))
    }
  }
  
  # Stop early if config is missing core pieces
  if (length(errors) > 0) {
    cat(red$bold("\n🚨 Critical configuration errors found. Cannot proceed.\n"))
    cat(paste(errors, collapse = "\n"))
    stop("🛑 Validation failed due to missing core config fields.")
  }
  
  # Create directories if needed
  for (d in c(config$general$out_dir, config$general$fig_dir, config$general$data_dir)) {
    if (!dir.exists(d)) {
      vcat(yellow(paste0("📁 Creating directory: ", d)))
      dir.create(d, recursive = TRUE)
    }
  }
  
  # 4. Ensure the specified network exists
  if (is.null(config$networks) || !(wgcna_name %in% names(config$networks))) {
    stop(red$bold(paste0("❌ Network '", wgcna_name, "' not found in config$networks.")))
  }
  
  net_cfg <- config$networks[[wgcna_name]]
  processed_file <- file.path(config$general$out_dir, paste0(wgcna_name, '/', config$general$data_dir, '/', wgcna_name, "_hdWGCNA.rds"))
  
  # 5. Determine recompute mode (now per-network)
  recompute <- if (!is.null(net_cfg$recompute_network)) {
    isTRUE(net_cfg$recompute_network)
  } else {
    TRUE  # default: recompute if not specified
  }
  vcat(blue$bold(paste0("⚙️ Recompute network for '", wgcna_name, "': ", recompute)))
    
  # 6. Load appropriate Seurat object
  if (recompute) {
    vcat(blue$bold("📦 Recompute requested — attempting to load raw Seurat object..."))

    # Try reading raw object
    seu <- tryCatch(
      readRDS(config$seurat_path),
      error = function(e) stop(red$bold("❌ Failed to load raw Seurat object:"), e$message)
    )

    if (!inherits(seu, "Seurat"))
      stop("❌ Object loaded is not a Seurat object.")

    vcat(green(paste0("✅ Loaded raw Seurat object with ", ncol(seu), " cells.")))

  } else {

    vcat(blue$bold("📦 Skipping recomputation — attempting to load processed hdWGCNA object..."))

    if (!file.exists(processed_file)) {

      vcat(yellow$bold("⚠️ Processed network object not found."))
      vcat(yellow("➡️ Falling back to raw Seurat object instead (recompute forced)."))

      # Fall back to raw Seurat object
      seu <- tryCatch(
        readRDS(config$seurat_path),
        error = function(e) stop(red$bold("❌ Failed to load raw Seurat object:"), e$message)
      )

      if (!inherits(seu, "Seurat"))
        stop("❌ Object loaded is not a Seurat object.")

      # Override recompute flag since we must recompute
      recompute <- TRUE

      vcat(green(paste0("✅ Loaded raw Seurat object with ", ncol(seu), " cells.")))

    } else {

      # Normal skip-recompute path
      seu <- tryCatch(
        readRDS(processed_file),
        error = function(e) stop(red$bold("❌ Failed to load processed hdWGCNA object:"), e$message)
      )

      if (!inherits(seu, "Seurat"))
        stop("❌ Object loaded is not a Seurat object.")

      vcat(green(paste0("✅ Loaded processed hdWGCNA object: ", wgcna_name)))
    }
  }

  
  # 7. Validate SetupForWGCNA
  if (!is.null(config$SetupForWGCNA)) {
    s_cfg <- config$SetupForWGCNA
    if (!is.null(s_cfg$assay) && !(s_cfg$assay %in% Assays(seu))) {
      errors <- c(errors, paste0("❌ Assay '", s_cfg$assay, "' not found."))
    }
    if (!is.null(s_cfg$group.by)) {
      errors <- validate_metadata_column(seu, s_cfg$group.by, context = "SetupForWGCNA",
                                         verbose = verbose, errors = errors)
    }
  }
  
  # 8. Validate MetacellsByGroups
  if (!is.null(config$MetacellsByGroups)) {
    m_cfg <- config$MetacellsByGroups
    for (field in c("group.by", "ident.group")) {
      if (!is.null(m_cfg[[field]])) {
        for (col in unlist(m_cfg[[field]])) {
          errors <- validate_metadata_column(seu, col, context = "MetacellsByGroups",
                                             verbose = verbose, errors = errors)
        }
      }
    }
    if (!is.null(m_cfg$reduction) && !(m_cfg$reduction %in% Reductions(seu))) {
      errors <- c(errors, paste0("❌ Reduction '", m_cfg$reduction, "' not found in Seurat object."))
    }
  }
  
  # 9. Validate the chosen network only
  vcat(blue$bold(paste0("🌐 Checking network: ", wgcna_name)))
  required_fields <- c("group.by", "replicate",  "subset")
  missing <- setdiff(required_fields, names(net_cfg))
  if (length(missing) > 0) {
    errors <- c(errors, paste0("❌ Network '", wgcna_name, "' missing fields: ",
                               paste(missing, collapse = ", ")))
  }
  for (f in c("group.by", "replicate")) {
    if (!is.null(net_cfg[[f]])) {
      errors <- validate_metadata_column(seu, net_cfg[[f]], context = wgcna_name,
                                         verbose = verbose, errors = errors)
    }
  }
  if (!is.null(net_cfg$subset)) {
    for (sub_key in names(net_cfg$subset)) {
      subset_values <- unlist(net_cfg$subset[[sub_key]])
      errors <- validate_metadata_column(seu, sub_key, values = subset_values,
                                         context = paste0(wgcna_name, " subset"),
                                         verbose = verbose, errors = errors)
    }
  }
  
  # 10. Final report
  if (length(errors) > 0) {
    cat(red$bold("\n❌ Validation failed with the following issues:\n"))
    cat(paste(errors, collapse = "\n"))
    stop("🛑 Config validation failed. Please correct the errors above.")
  } else {
    cat(green$bold("\n🎉 All validation checks passed successfully!"), "\n")
  }
  
  return(list(
    config = config,
    seurat = seu,
    wgcna_name = wgcna_name,
    processed_file = processed_file,
    recompute = recompute
  ))
}

make_markdown_summary <- function(options_list, header = "### Options:") {
  lines <- c(header, "")  # header + blank line
  
  for (key in names(options_list)) {
    value <- options_list[[key]]
    
    # Skip completely NULL values (or remove this if you want to show them)
    if (is.null(value)) {
      pretty_val <- "NULL"
    } else if (length(value) > 1) {
      # list or vector
      pretty_val <- paste(value, collapse = ", ")
    } else {
      pretty_val <- as.character(value)
    }
    
    pretty_key <- gsub("_", " ", key)
    
    lines <- c(
      lines,
      sprintf("* **%s**: %s", pretty_key, pretty_val)
    )
  }
  
  paste(lines, collapse = "\n")
}

make_markdown_summary <- function(
  options_list,
  header = "",
  format = c("table", "bullets"),
  indent = 0
) {
  format <- match.arg(format)

  # -------- Helper: format nested lists recursively ----------
  format_item <- function(key, value, indent, format) {
    indent_str <- paste0(rep("  ", indent), collapse = "")

    # Replace NULL with "default"
    if (is.null(value)) {
      display_val <- "default"
    } else if (is.list(value) && !is.null(names(value))) {
      # Named nested list
      if (format == "bullets") {
        out <- sprintf("%s* **%s**:", indent_str, key)
        children <- unlist(
          lapply(
            names(value),
            function(k) format_item(k, value[[k]], indent + 1, format)
          )
        )
        return(c(out, children))
      } else {
        child_rows <- unlist(
          lapply(
            names(value),
            function(k) format_item(k, value[[k]], indent + 1, "table")
          )
        )
        return(child_rows)
      }
    } else if (length(value) > 1) {
      display_val <- paste(value, collapse = ", ")
    } else {
      display_val <- as.character(value)
    }

    pretty_key <- gsub("_", " ", key)

    if (format == "bullets") {
      return(sprintf("%s* **%s**: %s", indent_str, pretty_key, display_val))
    } else {
      key_with_indent <- paste0(
        paste0(rep("&nbsp;&nbsp;", indent), collapse = ""),
        pretty_key
      )
      return(sprintf("| %s | %s |", key_with_indent, display_val))
    }
  }

  # -------- Build output ----------
  if (format == "bullets") {
    lines <- c(header, "")
    for (key in names(options_list)) {
      val <- options_list[[key]]
      lines <- c(lines, format_item(key, val, indent, "bullets"))
    }
    return(paste(lines, collapse = "\n"))
  }

  # TABLE FORMAT
  header_lines <- c(header, "", "| Option | Value |", "|--------|--------|")

  table_rows <- unlist(
    lapply(
      names(options_list),
      function(k) format_item(k, options_list[[k]], indent, "table")
    )
  )

  return(paste(c(header_lines, table_rows), collapse = "\n"))
}

# Plotly requires RGBA colors for opacity, so convert module hex → rgba()
hex_to_rgba <- function(hex, alpha) {
  rgb <- col2rgb(hex)
  sprintf("rgba(%d,%d,%d,%.3f)", rgb[1], rgb[2], rgb[3], alpha)
}


