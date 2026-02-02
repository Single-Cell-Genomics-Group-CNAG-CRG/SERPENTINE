#---------------------------------------------------------#
# Main functions to create the 
# Pseudobulk SummarizedExperiment object
#---------------------------------------------------------#

#' AggregatePseudobulk
#'
#' Create pseudobulk samples by aggregating single-cell (or single-nucleus) counts
#' according to replicate and group annotations, and return a SummarizedExperiment.
#'
#' @description
#' This function aggregates a gene-by-cell count matrix into gene-by-pseudobulk
#' counts. Pseudobulk groups are defined by the interaction of a replicate
#' identifier and a group identifier (for example: sample_id × cell_type). The
#' function builds a sparse design matrix that maps cells to pseudobulks,
#' multiplies the counts matrix by that mapping to obtain aggregated counts,
#' filters pseudobulks with too few contributing cells, removes genes with zero
#' variance across the kept pseudobulks, and returns a SummarizedExperiment
#' containing assay(s) and per-pseudobulk metadata (including nCells, nUMI and
#' nFeatures).
#'
#' @param X matrix or Matrix
#'   Gene-by-cell count matrix. Can be a base R matrix or a sparse Matrix
#'   (from the Matrix package). Columns must be cell identifiers that match
#'   the row names of `meta`.
#' @param meta data.frame
#'   Per-cell metadata. Row names must correspond to column names of `X`.
#'   Must contain the columns specified by `replicate_col` and `group_col`.
#' @param replicate_col character(1)
#'   Name of the column in `meta` indicating the biological replicate (for
#'   example sample or individual). Used as the first component of the
#'   interaction that defines pseudobulks.
#' @param group_col character(1)
#'   Name of the column in `meta` indicating the grouping factor (for example
#'   cell type, cluster, condition). Used as the second component of the
#'   interaction that defines pseudobulks.
#' @param min_cells integer(1), optional
#'   Minimum number of cells required for a pseudobulk to be retained. Pseudobulks
#'   with strictly greater than `min_cells` contributing cells are kept. Default
#'   value is 10. NOTE: the current implementation contains an internal
#'   assignment `min_cells <- 10` which overrides this argument — see Details.
#' @param assay_name character(1), optional
#'   Name to assign to the assay in the returned SummarizedExperiment. Default
#'   is "counts".
#'
#' @return SummarizedExperiment
#'   An object with:
#'   - assays: a named list with a single matrix-like assay (genes × pseudobulks)
#'     containing aggregated counts. The assay name equals `assay_name`.
#'   - colData: a data.frame with one row per pseudobulk (metadata built by
#'     `make_pseudobulk_metadata(meta, pb_groups)` and subset to kept pseudobulks).
#'   Additional columns added to colData:
#'     - nCells: number of cells that contributed to each pseudobulk
#'     - nUMI: sum of counts across genes for each pseudobulk
#'     - nFeatures: number of genes with non-zero counts in the pseudobulk
#'
#' @details
#' - Input checks:
#'   - `X` must be matrix-like (dense matrix or Matrix sparse object).
#'   - `meta` must be a data.frame with rownames matching `colnames(X)`.
#'   - `replicate_col` and `group_col` must exist in `meta` and contain no NAs.
#' - Grouping:
#'   - Pseudobulk groups are created by interaction(meta[[replicate_col]],
#'     meta[[group_col]], drop = TRUE). This produces factor levels representing
#'     unique replicate × group combinations.
#'   - A sparse model matrix (~0 + pb_groups) is constructed to map cells to
#'     pseudobulks. Columns of the resulting aggregated matrix are renamed by
#'     removing the "pb_groups" prefix that is created by the model matrix.
#' - Filtering:
#'   - The function computes the number of cells per pseudobulk (`n_cells`) and
#'     keeps only pseudobulks with n_cells > min_cells (strict inequality).
#'   - Genes with zero standard deviation across the retained pseudobulks are
#'     removed.
#' - Post-processing:
#'   - The returned SummarizedExperiment's colData receives nCells, nUMI and
#'     nFeatures. nUMI is computed as column sums of the aggregated counts.
#'     nFeatures is computed after thresholding counts to binary (counts > 1
#'     set to 1) and summing per column.
#'
#' @section Important note (implementation caveat):
#' The current function body contains an explicit assignment `min_cells <- 10`
#' after the per-pseudobulk cell counts are computed. This overwrites the value
#' provided by the user via the function argument, so passing a different
#' `min_cells` will have no effect unless the implementation is corrected.
#'
#' @section Edge cases and warnings:
#' - If column names of `X` are not present in rownames(meta), the function
#'   will stop and report the mismatch (it reports the missing cell ids).
#' - If all pseudobulks are filtered out by the `min_cells` threshold, the
#'   function will produce an empty SummarizedExperiment or fail in downstream
#'   steps; callers should check the returned object.
#' - The function assumes the existence of a helper `make_pseudobulk_metadata()`
#'   in the calling environment or package; this function must accept the same
#'   `meta` and `pb_groups` and return rownames corresponding to pseudobulk
#'   column order.
#' 
#' @seealso
#' make_pseudobulk_metadata, SummarizedExperiment, Matrix::sparse.model.matrix
#'
#' @export
AggregatePseudobulk <- function(
    X,
    meta,
    replicate_col, 
    group_col,
    min_cells = 10,
    assay_name = 'counts'
){

    # ------------------------------------------------------------
    # Input sanity checks
    # ------------------------------------------------------------

    # X must be a matrix-like object
    if (!inherits(X, "Matrix") && !is.matrix(X)) {
        stop("'X' must be a dense matrix or a sparse 'Matrix' object from the Matrix package.")
    }

    # meta must be a data.frame-like object
    if (!is.data.frame(meta)) {
        stop("'meta' must be a data.frame")
    }

    # Ensure that X and meta have matching cells
    if (!all(colnames(X) %in% rownames(meta))) {
        missing <- setdiff(colnames(X), rownames(meta))
        stop("Mismatch between cells in colnames(X) and rownames(meta)")
    }

    # Optionally reorder meta to match X
    meta <- meta[colnames(X), , drop = FALSE]

    cols_to_check <- c(replicate_col, group_col)

    for (col in cols_to_check) {

        # Ensure column exists
        if (!(col %in% colnames(meta))) {
            stop(sprintf("Column '%s' not found in meta.", col))
        }

        # Check for missing values
        if (any(is.na(meta[[col]]))) {
            stop(sprintf("Column '%s' contains NA values.", col))
        }
    }

    # define the pseudobulk grouping based on replicate_col
    pb_groups <- interaction(meta[,replicate_col], meta[,group_col], drop=TRUE)
    
    # calculate the number of cells per grouping
    n_cells <- table(pb_groups)

    # Create group indicator matrix: cells × pseudobulks
    G <- Matrix::sparse.model.matrix(~0 + pb_groups)

    # Pseudobulk = counts × G   (genes × cells  %*%  cells × groups)
    pb <- X %*% G
    colnames(pb) <- gsub('pb_groups', '', colnames(pb))

    # remove pseudobulk reps with insufficient cells
    pb <- pb[,as.logical(n_cells >= min_cells)]

    # calculate the standard deviation of each gene:
    good_genes <- names(which(apply(pb, 1, sd) != 0))
    pb <- pb[good_genes,]

    # create the pseudobulk meta-data table:
    pb_meta <- make_pseudobulk_metadata(meta, pb_groups)
    pb_meta <- pb_meta[colnames(pb),]

    assay_list <- list('tmp' = pb)
    names(assay_list) <- assay_name

    # Create a SummarizedExperiment object
    se <- SummarizedExperiment(
        assays = assay_list,
        colData = pb_meta
    )

    # add the number of cells
    colData(se)$nCells <- as.numeric(n_cells[colnames(se)])

    # nUMI and nFeatures detected:
    colData(se)$nUMI <- colSums(pb)
    pb[pb > 1] <- 1
    colData(se)$nFeatures <- colSums(pb)

    # return the SummarizedExperiment object:
    return(se)
}

# helper function for AggregatePseudobulk
find_replicate_columns <- function(meta, group){
  is_replicate_col <- function(col) {
    # For each group, check if all entries in that cluster-sample group are identical
    all(tapply(col, group, function(x) length(unique(x)) == 1))
  }
  
  replicate_cols <- names(meta)[sapply(meta, is_replicate_col)]
  replicate_cols
}

# helper function for AggregatePseudobulk
make_pseudobulk_metadata <- function(meta, group) {
  groups <- unique(group)

  replicate_cols <- find_replicate_columns(meta, group)
  
  # For every replicate col, extract ONE value per group
  out <- sapply(replicate_cols, function(colname) {
    tapply(meta[[colname]], group, function(x) unique(x)[1])
  })
  
  # Convert to data.frame, preserving group order
  out <- as.data.frame(out, row.names = levels(group))
  out
}

#' NormalizeCounts
#'
#' Perform per-sample normalization of a count assay stored in a
#' SummarizedExperiment and add the normalized matrix as a new assay.
#'
#' @param se A SummarizedExperiment containing a raw counts assay.
#' @param method Character; one of "CPM", "logCPM", "logNorm", "VST", "rlog".
#'   - "CPM": counts per million.
#'   - "logCPM": log2(CPM + pseudocount).
#'   - "logNorm": Seurat-style log1p(counts / size_factor * 1e4) where
#'     size_factor = colSums(counts) / median(colSums(counts)).
#'   - "VST", "rlog": variance-stabilizing transform or rlog via DESeq2.
#' @param assay_name Character scalar; name of the assay in `se` to normalize
#'   (default: "counts").
#' @param new_assay_name Character or NULL; name to assign the normalized assay.
#'   If NULL, defaults to the chosen `method`.
#' @param pseudocount Numeric scalar added to CPM before log2 in "logCPM"
#'   (default: 1).
#' @param ... Additional arguments forwarded to DESeq2::vst or DESeq2::rlog when
#'   `method` is "VST" or "rlog".
#'
#' @return A SummarizedExperiment identical to `se` but with a new assay named
#'   `new_assay_name` containing the normalized matrix.
#'
#' @details The function checks that `se` is a SummarizedExperiment and that
#'   `assay_name` exists and is a (possibly sparse) matrix. For "VST" and
#'   "rlog", DESeq2 must be installed; the function converts the assay to a
#'   dense matrix and constructs a DESeqDataSet with design ~ 1 before
#'   applying the transform. Errors are raised for invalid inputs.
#'
#' @examples
#' # Basic usage (assuming `se` is a SummarizedExperiment with a "counts" assay)
#' # se_norm <- NormalizeSE(se, method = "logCPM")
#'
#' @seealso DESeq2::vst, DESeq2::rlog
NormalizeCounts <- function(
    se,
    method = c("CPM", "logCPM", "logNorm", "VST", "rlog"),
    assay_name = "counts",
    new_assay_name = NULL,
    pseudocount = 1,
    ...
){
    method <- match.arg(method)

    # --- Check inputs ---------------------------------------------------------
    if (!inherits(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object.")
    }
    if (!assay_name %in% names(assays(se))) {
        stop(paste0("Assay '", assay_name, "' not found in se."))
    }

    X <- assay(se, assay_name)

    if (!is.matrix(X) && !inherits(X, "Matrix")) {
        stop("The assay must be a matrix or sparse Matrix.")
    }

    # default name for normalized assay
    if (is.null(new_assay_name)) {
        new_assay_name <- method
    }

    # --- Normalization methods -----------------------------------------------
    
    ## CPM ---------------------------------------------------------------------
    if (method == "CPM") {
        lib.size <- colSums(X)
        norm <- t(t(X) / lib.size) * 1e6
    }

    ## log CPM -----------------------------------------------------------------
    if (method == "logCPM") {
        lib.size <- colSums(X)
        cpm <- t(t(X) / lib.size) * 1e6
        norm <- log2(cpm + pseudocount)
    }

    ## Log-normalization (like Seurat: log1p(counts / size factor * 1e4)) ------
    if (method == "logNorm") {
        size.factor <- colSums(X) / median(colSums(X))
        scaled <- t(t(X) / size.factor) * 1e4
        norm <- log1p(scaled)
    }

    ## VST using DESeq2 --------------------------------------------------------
    if (method %in% c("VST", "rlog")) {
        if (!requireNamespace("DESeq2", quietly = TRUE)) {
            stop("DESeq2 must be installed for VST / rlog.")
        }

        # Extract the raw counts and convert to a dense matrix
        mat_dense <- as.matrix(SummarizedExperiment::assay(se, assay_name))

        # Rebuild a SE object with dense assay for DESeq2
        se_dense <- SummarizedExperiment::SummarizedExperiment(
            assays  = list(counts = mat_dense),
            colData = SummarizedExperiment::colData(se)
        )

        # Construct DESeqDataSet
        dds <- DESeq2::DESeqDataSet(se_dense, design = ~ 1)

        # Run DESeq2 normalization transform
        if (method == "VST") {
            norm <- SummarizedExperiment::assay(DESeq2::vst(dds, ...))
        } else {  # rlog
            norm <- SummarizedExperiment::assay(DESeq2::rlog(dds, ...))
        }
    }

    # --- Add normalized assay and return --------------------------------------
    assay(se, new_assay_name) <- norm
    return(se)
}




#---------------------------------------------------------#
# WGCNA-specific functions
#---------------------------------------------------------#

#' Compute module eigengenes (MEs) from a SummarizedExperiment
#'
#' Compute module eigengenes for a set of gene modules using WGCNA::moduleEigengenes
#' and return them as a new SummarizedExperiment assay.
#'
#' @param se SummarizedExperiment containing expression assays (rows = genes, cols = samples).
#' @param modules data.frame with at least columns 'gene_name' and 'module' mapping genes to modules.
#' @param assay_name character(1) name of the assay in `se` to use for expression (default: "VST").
#' @param new_assay_name character(1) name for the assay in the returned SummarizedExperiment (default: "MEs").
#' @return SummarizedExperiment with a single assay of module eigengenes (rows = modules, cols = samples)
#' and the original colData preserved.
#' @details Only genes present in the specified assay are used. Eigengenes are computed per module using
#' WGCNA::moduleEigengenes on the transposed expression matrix.
#' @export
PseudobulkModuleEigengenes <- function(
    se,
    modules,
    assay_name = 'VST',
    new_assay_name = 'MEs'
){

    # TODO: checks

    # get a list of modules to calculate
    mods <- unique(as.character(modules$module))

    # get the expression matrix:
    X <- assays(se)[[assay_name]]

    # loop through each module and calculate the ME
    eigengene_list <- lapply(mods, function(mod) {
        genes <- subset(modules, module == mod) %>% .$gene_name
        genes <- intersect(genes, rownames(X))
        colors <- rep(mod, length(genes)); names(colors) <- genes
        ME <- WGCNA::moduleEigengenes(
            t(X[genes, , drop=FALSE]), 
            colors = colors
        )$eigengenes
        return(ME[[1]])
    })

    # create the MEs matrix
    MEs <- as.data.frame(do.call(rbind, eigengene_list))
    rownames(MEs) <- mods; colnames(MEs) <- colnames(X)

    # create a new SummarizedExperiment:
    assay_list <- list('tmp' = MEs)
    names(assay_list) <- new_assay_name
    se_modules <- SummarizedExperiment(
        assay_list,
        colData = colData(se)
    )

    return(se_modules)

}

#---------------------------------------------------------#
# Principal Component analysis
#---------------------------------------------------------#

#' Run PCA on a SummarizedExperiment assay and store results in metadata
#'
#' Performs principal component analysis (via prcomp) on the specified assay of a
#' SummarizedExperiment (rows = features, cols = samples) and saves a list with
#' embedding (PC scores), loadings (rotation) and variance (proportion per PC)
#' into metadata(se)[[name]].
#'
#' @param se SummarizedExperiment object containing the assay to analyze.
#' @param features Character vector of feature names to use; if NULL all features are used.
#' @param assay_name Character; name of the assay in se containing expression values (default "VST").
#' @param name Character; metadata slot name to store PCA results (default "PCA").
#' @param n_components Integer; number of principal components to compute (default 50).
#' @param scale Logical; whether to scale features to unit variance before PCA (default TRUE).
#' @param center Logical; whether to center features before PCA (default TRUE).
#' @param overwrite Logical; if FALSE and metadata slot exists, the function errors (default FALSE).
#'
#' @return The input SummarizedExperiment with metadata[[name]] set to a list:
#'   - embedding: matrix of PC scores (samples x n_components)
#'   - loadings: matrix of variable loadings (features x n_components)
#'   - variance: numeric vector of variance proportion for each component
#'
#' @details
#' Validates input object and assay presence, checks that requested features/samples
#' are sufficient for the requested number of components, and requires no NA values
#' in the expression matrix. Uses prcomp(t(X)) internally.
#'
#' @examples
#' # se <- SummarizedExperiment(assays = list(VST = exprs))
#' # se <- PseudobulkPCA(se, features = rownames(se), assay_name = "VST", n_components = 10)
#'
#' @export
PseudobulkPCA <- function(
    se, 
    features = NULL,
    assay_name = "VST",
    reduction_name = "PCA",
    n_components = 50,
    scale = TRUE,
    center = TRUE,
    overwrite = FALSE 
){
    
    # Check if the input object is a SummarizedExperiment
    if (!is(se, "SummarizedExperiment")) {
        stop("Input object 'se' must be a SummarizedExperiment.")
    }
    
    # Check if the specified assay exists
    if (!(assay_name %in% assayNames(se))) {
        stop(paste0("Assay '", assay_name, "' not found in the SummarizedExperiment."))
    }
        
    if (reduction_name %in% names(metadata(se)) && !overwrite) {
            stop(paste0("The metadata slot '", reduction_name, 
                        "' already exists. Set overwrite = TRUE to replace it."))
        }

    # Check if the specified features exist in the assay
    if (!is.null(features)) {
        missing_features <- features[!(features %in% rownames(se))]
        if (length(missing_features) > 0) {
            stop(paste0("The following features are missing from the SummarizedExperiment: ", 
                        paste(head(missing_features, 5), collapse = ", "), 
                        if (length(missing_features) > 5) "..." else "."))
        }
    } else {
        # Default to all features if not specified
        features <- rownames(se)
    }

    # Check if the number of features is sufficient for the requested rank
    if (length(features) < n_components) {
        stop(paste0("Number of features (", length(features), 
                    ") is less than the requested number of components (", n_components, ")."))
    }

    # Check if the number of samples is sufficient for the requested rank
    if (ncol(se) < n_components) {
        stop(paste0("Number of samples (", ncol(se), 
                    ") is less than the requested number of components (", n_components, ")."))
    }
    
    # get the expression matrix 
    X <- assays(se)[[assay_name]][features,]
    if (any(is.na(X))) {
            stop("The expression matrix contains NA values. PCA requires complete data.")
    }

    # run pca
    pca_results <- prcomp(
        t(X), 
        rank. = n_components, 
        scale=scale
    )  

    # calculate percentage of variance
    pca_sd <- pca_results$sdev[1:n_components]
    pca_var <- (pca_sd^2) / sum(pca_sd^2)

    # add the results to the SummarizedExperiment object
    metadata(se)[[reduction_name]] <- list(
        embedding = pca_results$x,
        loadings = pca_results$rotation,
        variance = pca_var
    )

    return(se)

}






#' Plot PCA embedding for pseudobulk samples
#'
#' Create a 2D scatter plot of two principal components from a PCA reduction stored
#' in the metadata of a SummarizedExperiment-like object, with points colored by a
#' specified sample-level metadata column.
#'
#' @param se SummarizedExperiment-like object. Expects a PCA embedding at
#'   metadata(se)[[reduction_name]]$embedding and sample metadata in colData(se).
#' @param color_by Character. Name of the column in colData(se) used to color points.
#' @param reduction_name Character. Name of the reduction stored in metadata(se)
#'   (default: "PCA").
#' @param components Character vector of length 2. Names of the embedding columns
#'   to plot (default: c("PC1", "PC2")).
#' @param raster_dpi Resolution for rasterising the points
#' @param point_size Numeric. The size of the points. Default is 2.
#' 
#' @details
#' The function validates that the requested components exist in the embedding and
#' that row names of the embedding match row names of colData(se). It will stop with
#' a clear error when these conditions are not met.
#'
#' @return A ggplot2 object representing the PCA scatter plot (x = components[1],
#'   y = components[2], colored by color_by).
#'
#' @examples
#' # PlotPCAEmbedding(se, color_by = "sample_group")
#'
#' @importFrom ggplot2 ggplot aes geom_point theme xlab ylab labs ggtitle
#' @export
PlotPCAEmbedding <- function(
    se, 
    color_by,
    reduction_name = "PCA", 
    components = c("PC1", "PC2"),
    raster_dpi = 300,
    point_size = 2
) {
    # Get the PCA embedding matrix from the metadata slot
    pca_mat <- metadata(se)[[reduction_name]]$embedding
    pca_df <- as.data.frame(pca_mat)

    # Check if the requested components exist in the embedding
    if (!all(components %in% colnames(pca_mat))) {
        stop(paste0("Requested components (", paste(components, collapse = ", "), ") not found in the embedding matrix."))
    }
    pca_df <- pca_df[, components]
    
    # Get the pseudo-bulk metadata (colData)
    meta_df <- as.data.frame(colData(se))
    
    # Check if the color_by column exists in colData
    if (!(color_by %in% colnames(meta_df))) {
        stop(paste0("Color variable '", color_by, "' not found in colData(se)."))
    }
        
    # Merge PCA coordinates with metadata
    # Ensure rownames/sample names match before combining
    if (!all(rownames(pca_df) == rownames(meta_df))) {
        stop("Rownames of PCA embedding and colData(se) do not match.")
    }
    
    # Select the required columns and merge
    plot_df <- meta_df %>%
        dplyr::select(all_of(color_by)) %>%
        cbind(pca_df)
    
    # Build the ggplot
    p <- plot_df %>% 
        ggplot(aes(x = .data[[components[1]]], y = .data[[components[2]]])) + 
        ggrastr::rasterise(
            geom_point(aes(color = .data[[color_by]]), size = point_size),
        dpi = raster_dpi) + 
        theme(
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(linewidth = 1, color = 'black', fill = NA),
            panel.grid.major.y = element_line(linewidth = 0.25, color = 'lightgrey'),
            panel.grid.major.x = element_line(linewidth = 0.25, color = 'lightgrey'),
            plot.title = element_text(hjust = 0.5) 
        ) + 
        xlab(components[1]) + 
        ylab(components[2]) +
        labs(color = color_by) 
    
    return(p)
}

#' Plot Continuous Features on a PCA Embedding
#'
#' @description
#' Visualizes continuous variables—either from sample metadata (colData) or gene 
#' expression (assays)—on a 2D PCA coordinate system stored within a 
#' \code{SummarizedExperiment} object. 
#'
#' @param se A \code{SummarizedExperiment} object containing the PCA results in its 
#' \code{metadata} slot.
#' @param feature Character string. The name of the feature to color points by. 
#' The function first searches \code{colData(se)} for a matching column name. 
#' If not found, it searches \code{rownames(se)} for a gene/row name.
#' @param reduction_name Character string. The key in \code{metadata(se)} where 
#' the PCA list is stored. Default is "PCA".
#' @param assay_name Character string. The assay to extract values from if 
#' \code{feature} is a row name (e.g., "counts", "logcounts", "VST"). Default is "counts".
#' @param components A character vector of length 2 specifying which principal 
#' components to plot. Default is \code{c("PC1", "PC2")}.
#' @param raster_dpi Numeric. The resolution in dots per inch for rasterizing 
#' the point layer using \code{ggrastr}. Useful for large datasets. Default is 300.
#' @param point_size Numeric. The size of the points. Default is 2.
#' @param high_color Character string. The color for high values in the gradient.
#' Default is 'darkorchid4'.
#' @param low_color Character string. The color for low values in the gradient.
#' Default is 'whitesmoke'.
#' @param order Logical. If \code{TRUE} (default), points are reordered by their 
#' value so that samples with higher values are plotted on top. This prevents 
#' "overplotting" where high-signal points are hidden.
#' @param adjust_outliers Numeric or NA, quantile threshold for adjusting extreme significance 
#'    values (default: 0.999, set to NA to disable).
#' @details
#' The function enforces that the \code{feature} must be numeric. If a categorical 
#' variable from \code{colData} is passed, the function will stop and suggest 
#' using \code{PlotPCAEmbedding} instead. Point coordinates are pulled from the 
#' \code{embedding} slot of the specified reduction in \code{metadata(se)}.
#'
#' @return A \code{ggplot} object with the points rasterized according to 
#' \code{raster_dpi}.
#'
#' @export
#' @import ggplot2
#' @importFrom ggrastr rasterise
#' @importFrom SummarizedExperiment assay colData metadata
PlotPCAFeatures <- function(
    se, 
    feature,
    reduction_name = "PCA", 
    assay_name = "counts",
    components = c("PC1", "PC2"),
    raster_dpi = 300,
    point_size = 2,
    high_color = 'darkorchid4',
    low_color = 'whitesmoke',
    order = TRUE,
    adjust_outliers = 0.999
) {

    if (!inherits(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object.")
    }

   # Get the PCA embedding matrix from the metadata slot
    pca_mat <- metadata(se)[[reduction_name]]$embedding
    pca_df <- as.data.frame(pca_mat)

    # Check if the requested components exist in the embedding
    if (!all(components %in% colnames(pca_mat))) {
        stop(paste0("Requested components (", paste(components, collapse = ", "), ") not found in the embedding matrix."))
    }
    pca_df <- pca_df[, components]

    # check if the feature is in the colData or in the assay 
    if (feature %in% colnames(colData(se))){
        plot_df <- data.frame(
            value = colData(se)[[feature]]
        )
        plot_label <- feature
        if(!is.numeric(plot_df$value)){
            stop(paste0("Feature '", feature, "' in colData is not numeric. If you want to plot categorical data, please use the PlotPCAEmbedding function."))
        }
    } else if (feature %in% rownames(se)){
        if (!assay_name %in% names(assays(se))) {
            stop(paste0("Assay '", assay_name, "' not found in SummarizedExperiment object."))
        }
        plot_df <- data.frame(
            value = assay(se, assay_name)[feature,]
        )
        plot_label <- paste0(feature, "\n(", assay_name, ")")
    } else{
        stop(paste0("Feature '", feature, "' not found in colData or rownames of SummarizedExperiment object."))
    }

    # create the plot df with the data and embedding 
    plot_df <- cbind(plot_df, pca_df)

    # adjust outliers 
    if(!is.na(adjust_outliers)){
        new_max <- quantile(plot_df$value, adjust_outliers)
        plot_df$value <- ifelse(plot_df$value > new_max, new_max, plot_df$value)
    }

    # optionally order by value to plot high values on top
    if(order){
        plot_df <- plot_df[order(plot_df$value), ]
    }

    # Build the ggplot
    p <- plot_df %>% 
        ggplot(aes(x = .data[[components[1]]], y = .data[[components[2]]])) + 
        ggrastr::rasterise(
            geom_point(aes(color = value), size = point_size),
        dpi = raster_dpi) + 
        theme(
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(linewidth = 1, color = 'black', fill = NA),
            panel.grid.major.y = element_line(linewidth = 0.25, color = 'lightgrey'),
            panel.grid.major.x = element_line(linewidth = 0.25, color = 'lightgrey'),
            plot.title = element_text(hjust = 0.5) 
        ) + 
        xlab(components[1]) + 
        ylab(components[2]) +
        labs(color = plot_label) + 
        scale_color_gradient(low='whitesmoke', high='darkorchid4')
    
    return(p)
}

#' PCRegression
#'
#' Fit linear models of principal components against sample-level covariates and
#' store results in the SummarizedExperiment metadata.
#'
#' @param se A SummarizedExperiment containing PCA results in metadata(se)[[reduction_name]]$embedding.
#' @param covariates Character vector of column names in colData(se) to test as predictors.
#' @param reduction_name Character; name of the reduction stored in metadata(se) (default: "PCA").
#' @param overwrite Logical; if FALSE, do not overwrite existing regression results (default: FALSE).
#'
#' @details
#' - If a single covariate is provided, a "level" model is fit (pc ~ factor) for each PC and
#'   full coefficient tables (estimate, std. error, t, p) are returned in
#'   metadata(se)[[reduction_name]]$level_regression (with BH FDR correction).
#' - If multiple covariates are provided, separate simple regressions (pc ~ covariate) are fit
#'   for each PC and covariate; summary metrics (R2, adjusted R2, p, FDR) are stored in
#'   metadata(se)[[reduction_name]]$regression.
#' - Performs input validation: checks class of se, presence of reduction and embedding,
#'   existence and variance of covariates, and sample name alignment between embedding and colData.
#'
#' @return The input SummarizedExperiment with regression results added to metadata(se)[[reduction_name]].
#'
#' @examples
#' # PCRegression(se, covariates = c("Batch", "Age"))
#' # PCRegression(se, covariates = "Tissue", overwrite = TRUE)
#'
#' @keywords regression PCA SummarizedExperiment
PCRegression <- function(
    se,
    covariates,
    reduction_name = "PCA",
    overwrite = FALSE
){

    # Check if the input object is a SummarizedExperiment
    if (!is(se, "SummarizedExperiment")) {
        stop("Input object 'se' must be a SummarizedExperiment.")
    }
    
    # Check if the specified reduction_name (PCA) exists in the metadata
    if (!(reduction_name %in% names(metadata(se)))) {
        stop(paste0("Reduction name '", reduction_name, 
                    "' not found in metadata(se). Did you run PCA first?"))
    }
    
    # Check if the embedding matrix exists within the reduction
    if (!("embedding" %in% names(metadata(se)[[reduction_name]]))) {
        stop(paste0("Embedding matrix not found in metadata(se)$", reduction_name, 
                    ". Check the 'embedding' key."))
    }
    
    # Check if all specified covariates exist in the colData
    missing_covs <- covariates[!(covariates %in% colnames(colData(se)))]
    if (length(missing_covs) > 0) {
        stop(paste0("The following covariates are missing from colData(se): ", 
                    paste(missing_covs, collapse = ", ")))
    }

    # Determine Analysis Mode 
    level_mode <- length(covariates) == 1
    
    # Check if covariates has any usable values (no zero variance)
    # A simple check for zero variance (needed for lm)
    for (col in covariates) {
        if (length(unique(colData(se)[[col]])) <= 1) {
            stop(paste0("Covariate '", col, 
                        "' has zero or one unique value and cannot be used for regression."))
        }
    }

    # Check for overwrite conflict 
    if (!level_mode && ("regression" %in% names(metadata(se)[[reduction_name]])) && !overwrite) {
        stop(paste0("The regression slot for '", reduction_name, 
                    "' already exists. Set overwrite = TRUE to replace it."))
    }
    if (level_mode && ("level_regression" %in% names(metadata(se)[[reduction_name]])) && !overwrite) {
        stop(paste0("The level_regression slot for '", reduction_name, 
                    "' already exists. Set overwrite = TRUE to replace it."))
    }

    # get the PCA matrix
    pca_embedding <- metadata(se)[[reduction_name]]$embedding 
    metadata_df <- as.data.frame(colData(se)[, covariates, drop = FALSE]) # drop=FALSE ensures a DataFrame even with one column

    # Ensure the rows (samples) are perfectly aligned
    if (!all(rownames(pca_embedding) == rownames(metadata_df))) {
        stop("Sample names in PCA embedding and colData do not match.")
    }

    # Convert metadata columns to factors if they aren't already
    for (col in covariates) {
        if (is.character(metadata_df[[col]])) {
            metadata_df[[col]] <- factor(metadata_df[[col]])
        }
    }

    # Initialize the results table (PCs x Covariates)
    num_pcs <- ncol(pca_embedding)
    df <- data.frame()

    # Level Mode (Single Covariate)
    if (level_mode) {
        
        cov_name <- covariates[1] # The single covariate (e.g., "Tissue")
        levels_to_test <- levels(metadata_df[[cov_name]]) # Get all levels (e.g., "Liver", "Lung")
        
        # Initialize storage for coefficients and P-values
        # Output will be: PC x Model_Term (e.g., Intercept, TissueLung, TissueLiver)
        num_terms <- length(levels_to_test) 
        pc_names <- colnames(pca_embedding)
        
        # We will store the full model summary for simplicity
        results_list <- list() 

        for (i in 1:ncol(pca_embedding)) {
            pc_name <- pc_names[i]
            pc_scores <- pca_embedding[, i]
            
            # Use the factor name directly
            formula_str <- paste("pc_scores ~", cov_name)
            
            tryCatch({
                # Fit the standard LM (e.g., PC1 ~ Tissue)
                lm_fit <- lm(as.formula(formula_str), 
                             data = cbind(pc_scores = pc_scores, metadata_df))
                
                # Extract summary results
                s <- summary(lm_fit)$coefficients
                
                # Convert matrix of coefficients/stats to a data frame and store
                # This keeps all the coefficient information for the PC
                s_df <- as.data.frame(s)
                s_df$component <- i
                s_df$model_term <- rownames(s_df)
                
                results_list[[pc_name]] <- s_df
                
            }, error = function(e) {
                warning(paste("Error fitting Level Model for PC", i, "and covariate", cov_name, ":", e$message))
            })
        }
        
        # Combine all PC results into one tidy data frame
        df <- do.call(rbind, results_list)
        
        # Clean up column names (Example: Estimate is the coefficient/effect size)
        df <- df %>% 
            dplyr::rename(
                estimate = Estimate, 
                stderr = `Std. Error`,
                tval = `t value`, 
                pval = `Pr(>|t|)`
            ) %>%
            dplyr::select(component, model_term, estimate, pval)
        
        rownames(df) <- 1:nrow(df)

        # apply p-value correction
        df$fdr <- p.adjust(df$pval, method="BH")

        # Add to the SE metadata slot
        metadata(se)[[reduction_name]]$level_regression <- df
        
        return(se)
        
    } else {
        for (cov_name in covariates) {
            for (i in 1:num_pcs) {
                pc_scores <- pca_embedding[, i]
                
                # Construct the simple regression formula: PC_score ~ Covariate
                formula_str <- paste("pc_scores ~", cov_name)
                
                # Fit the Linear Model
                tryCatch({
                    lm_fit <- lm(as.formula(formula_str), 
                                data = cbind(pc_scores = pc_scores, metadata_df))

                    s <- summary(lm_fit)
                    a <- anova(lm_fit)

                    cur_df <- data.frame(
                        component = i,
                        covariate = cov_name,
                        R2 = s$r.squared,
                        adj_R2 = s$adj.r.squared,
                        pval = a$`Pr(>F)`[1]
                    )
                    df <- rbind(df, cur_df)

                }, error = function(e) {
                    warning(paste("Error fitting model for PC", i, "and covariate", cov_name, ":", e$message))
                })
            }
        }

        # apply p-value correction
        df$fdr <- p.adjust(df$pval, method="BH")

        # add to the se metadata slot:
        metadata(se)[[reduction_name]]$regression <- df  
    }

    return(se)
}


# SignatureDistPlot <- function(
#     signature_df,
#     meta_df,
#     sample_col,
#     group_by,
#     signature_col = 'module',
#     score_col = 'score',
#     shape = NA,
#     comparisons = NA,
#     comparison_method = 'wilcox',
#     raster_dpi = 300,
#     point_size = 1.5,
#     box_alpha = 0.4,
#     box_fill = 'grey',
#     box_notch = TRUE,
#     box_width = 0.5,
#     signature_cp = NA,
#     ncol = 5
# ){

#     # TODO: checks 

#     # simplify the input metadata to just the columns that we need
#     meta_df <- meta_df %>%
#         dplyr::select(all_of(c(sample_col, group_by))) %>% 
#         dplyr::distinct()


#     # merge the signature df with the metadata df
#     plot_df <- as.data.frame(left_join(
#         as.data.frame(signature_df), 
#         as.data.frame(meta_df),
#         by = sample_col
#     ))

#     # check that signature_cp is valid 
#     signatures <- as.character(unique(plot_df[,signature_col]))

#     # check that signature_col is valid

#     # check that comparisons are valid based on what is in group_by
#     valid_groups <- as.character(unique(plot_df[,group_by]))

#     # set up the plot 
#     p <- plot_df %>% ggplot(aes(y = get(score_col), x = get(group_by)))

#     # plot the points on the bottom:
#     p <- p + 
#         ggrastr::rasterise(ggbeeswarm::geom_quasirandom(
#             aes(color=get(signature_col)),
#             method = "pseudorandom",
#             size = point_size
#         ), dpi=raster_dpi)

#     # add the box 
#     p <- p + 
#         geom_boxplot(
#             outlier.shape=NA, 
#             alpha = box_alpha,
#             width = box_width,
#             fill = box_fill,
#             notch = box_notch
#         ) 

#     # color scheme for each signature?
#     if(all(signatures %in% names(signature_cp))){
#         p <- p + scale_color_manual(values = signature_cp) 
#     }

#     # compare groups?
#     if(all(unlist(comparisons) %in% valid_groups)){
#         my_comparisons <- comparisons 
#         p <- p + 
#             ggpubr::stat_compare_means(
#                 aes(label = after_stat(p.signif)),
#                 comparisons = my_comparisons, 
#                 method = comparison_method
#             ) +
#             scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
#     } 

#     # set up the theme:
#     p <- p + 
#         theme(
#             axis.line.x = element_blank(),
#             axis.line.y = element_blank(),
#             panel.border = element_rect(linewidth=1,color='black', fill=NA),
#             panel.grid.major.y = element_line(linewidth=0.25, color='lightgrey'),
#             plot.title = element_text(hjust=0.5),
#             strip.background = element_blank(),
#             strip.text = element_text(face='bold')
#         ) + Seurat::NoLegend() + Seurat::RotatedAxis() +
#         xlab('') + ylab(score_col)  

#     # facet 
#     patch <- p + 
#         facet_wrap(
#             ~get(signature_col), ncol=ncol, scales='free'
#         )

#     patch
# } 
