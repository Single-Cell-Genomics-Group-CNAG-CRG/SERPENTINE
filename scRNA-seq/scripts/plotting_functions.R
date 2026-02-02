
#' PlotEmbedding
#'
#' Create a 2D embedding plot (e.g., UMAP or PCA) of cells in a Seurat object, 
#' colored by a specified metadata variable, with options for splitting, labeling, 
#' custom colors, and other plot customizations.
#'
#' @param seurat_obj A Seurat object containing the embedding and metadata to plot.
#' @param group.by Character. Metadata column name to color the points by (e.g., cluster, condition).
#' @param reduction Character. Name of the dimensionality reduction to use (default is `'umap'`).
#' @param split.by Character (optional). Metadata column to split the plot by (e.g., sample, condition).
#' @param plot_under Logical. If `TRUE` and `split.by` is set, plots remaining cells in grey under the split panel.
#' @param label Logical. If `TRUE`, adds text labels at cluster centroids.
#' @param point_size Numeric. Size of points in the plot.
#' @param legend_point_size Numeric. Size of points in the legend.
#' @param text_size Numeric. Size of label text.
#' @param raster Logical. If `TRUE`, rasterizes the points layer for faster plotting.
#' @param order_points Character. Point plotting order: `"shuffle"` randomizes the order to reduce overplotting.
#' @param x, y Character (optional). Metadata columns to use for custom x and y coordinates instead of a reduction.
#' @param raster_dpi Numeric. Resolution (DPI) for rasterized layer.
#' @param raster_scale Numeric. Scale factor for rasterized layer.
#' @param selected Character vector (optional). Subset of groups to highlight; others will be hidden.
#' @param plot_theme A ggplot2 theme to apply to the plot (optional).
#' @param color_df Data frame (optional). Custom colors for groups. Must contain `group` and `colour` columns.
#' @param plot_ratio Logical. If `TRUE`, keeps aspect ratio square.
#'
#' @return A `ggplot` object if `split.by` is `NULL`, otherwise a list of `ggplot` objects, one per split group.
#'
#' @details  
#' - If `x` and `y` are provided, these override the selected `reduction`.  
#' - Cluster centroids are computed for labeling if `label = TRUE`.  
#' - If `split.by` is provided, separate panels are generated for each split group, 
#' optionally overlaying the remaining cells underneath.
#'
#' @examples
#' PlotEmbedding(seurat_obj = pbmc, group.by = "seurat_clusters")
#' PlotEmbedding(seurat_obj = pbmc, group.by = "seurat_clusters", split.by = "sample")
#'
#' @export
PlotEmbedding <- function(
  seurat_obj,
  group.by,
  reduction = 'umap',
  split.by = NULL,
  plot_under = FALSE,
  label = TRUE,
  point_size = 1,
  legend_point_size = 5,
  text_size = 3,
  raster = TRUE,
  order_points = "shuffle",
  x = NULL,
  y = NULL,
  raster_dpi = 400,
  raster_scale = 1,
  selected = NULL,
  plot_theme = NULL,
  color_df = NULL,
  plot_ratio=TRUE
){

  plot_df <- seurat_obj@meta.data
  if(!is.null(x) & !is.null(y)){
    plot_df$x <- plot_df[[x]]
    plot_df$y <- plot_df[[y]]
  } else if(reduction %in% names(seurat_obj@reductions)){
    plot_df$x <- seurat_obj@reductions[[reduction]]@cell.embeddings[,1]
    plot_df$y <- seurat_obj@reductions[[reduction]]@cell.embeddings[,2]
  }


  # convert to a factor:
  if(!is.factor(plot_df[[group.by]])){
    cur_groups <- unique(plot_df[[group.by]])
    cur_groups <- cur_groups[order(cur_groups)]
    plot_df[[group.by]] <- factor(
      as.character(plot_df[[group.by]]),
      levels = cur_groups
    )
  }

  # compute coordinates for cluster labels
  centroid_df <- data.frame()
  if(label){
    for(cur_cluster in unique(plot_df[[group.by]])){
      cur_meta <- plot_df[plot_df[[group.by]] == cur_cluster,]
      df <- data.frame(
        cluster = cur_cluster,
        x = mean(cur_meta$x),
        y = mean(cur_meta$y)
      )
      centroid_df <- rbind(centroid_df, df)
    }
  }

  # make a dummy ggplot to extract color scheme
  if(is.null(color_df)){

    factor_df <- data.frame(
      level = 1:length(levels(plot_df[[group.by]])),
      group_name = levels(plot_df[[group.by]])
    )

    p <- plot_df %>%
      ggplot(aes_string(x='x', y='y', color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    color_df$group <- factor_df$group_name

  }

  # only show selected groups
  if(!is.null(selected)){
    plot_df[[group.by]][!(plot_df[[group.by]] %in% selected)] <- NA

    if(is.factor(plot_df[[group.by]])){
      plot_df[[group.by]] <- droplevels(plot_df[[group.by]])
    }

    if(label){
      centroid_df <- subset(centroid_df, cluster %in% selected)
    }
  }


  # shuffle points:
  if(order_points == "shuffle"){
    plot_df <- plot_df[sample(nrow(plot_df)),]
  }

  # plot a single embedding
  if(is.null(split.by)){
    p <- .PlotSingleEmbedding(plot_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, plot_theme=plot_theme, plot_ratio=plot_ratio)
  } else{

    split_groups <- unique(plot_df[[split.by]])
    plot_list <- lapply(split_groups, function(cur_split){
      cur_df <- plot_df[plot_df[[split.by]] == cur_split,]
      split_df <- plot_df[plot_df[[split.by]] != cur_split,]
      .PlotSingleEmbedding(cur_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, split_df, cur_split, plot_theme, plot_under, plot_ratio=plot_ratio)
    })

    names(plot_list) <- split_groups

    return(plot_list)

  }

  p

}

.PlotSingleEmbedding <- function(
  plot_df,
  group.by,
  label,
  raster,
  raster_dpi,
  raster_scale,
  point_size,
  legend_point_size,
  text_size,
  color_df,
  centroid_df,
  split_df = NULL,
  cur_split = NULL,
  plot_theme = NULL,
  plot_under=FALSE,
  plot_ratio = TRUE
){

  p <- plot_df %>%
    ggplot(aes_string(x='x', y='y', color=group.by))

  # add the under layer for the split plots
  if(!is.null(split_df) & plot_under){
    # add points
    if(raster){
      p <-  p + ggrastr::rasterise(geom_point(data = split_df, size=point_size/2, color='lightgrey'), dpi=raster_dpi, scale=raster_scale)
    } else{
      p <- p + geom_point(data = split_df, size=point_size/2, color='lightgrey')
    }
  }

  # add points
  if(raster){
    p <-  p + ggrastr::rasterise(geom_point(size=point_size), dpi=raster_dpi,  scale=raster_scale)
  } else{
    p <- p + geom_point(size=point_size)
  }

  # add labels
  if(label){
    p <- p + ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=text_size)
  }

  # color scheme
  plot_colors <- color_df$colour
  names(plot_colors) <- as.character(color_df$group)
  p <- p + scale_color_manual(values=plot_colors, na.value = 'grey90')

  # add title:
  if(is.null(cur_split)){
    p <- p + ggtitle(group.by) + theme(plot.title=element_text(hjust=0.5))
  } else{
    p <- p + ggtitle(cur_split) + theme(plot.title=element_text(hjust=0.5))
  }

  # add theme:
  if(!is.null(plot_theme)){
    p <- p + plot_theme
  }

  # fixed coords
  if(plot_ratio){
    p <- p + coord_equal()
  }

  # adjust legend point size
  p <- p + guides( color = guide_legend(override.aes = list(size=legend_point_size)))

  p

}

SetupModuleNetwork <- function(
    seurat_obj,
    n_hubs = 25,
    edge_multiplier = 10,
    wgcna_name = NULL
) {

  require(igraph)
  require(dplyr)
  require(reshape2)

  # Resolve WGCNA name if NULL
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  n_edges <- n_hubs * edge_multiplier

  # ---- Modules + Colors ----
  modules <- GetModules(seurat_obj, wgcna_name = wgcna_name) %>%
    subset(module != "grey") %>%
    mutate(module = droplevels(module))

  mods <- levels(modules$module)

  mod_colors <- modules %>%
    select(module, color) %>% 
    distinct()

  mod_cp <- mod_colors$color
  names(mod_cp) <- as.character(mod_colors$module)

  # ---- Load TOM ----
  TOM <- GetTOM(seurat_obj, wgcna_name = wgcna_name)

  TOM_df <- data.frame()
  hub_df <- data.frame()

  # ---- Build TOM edges across all modules ----
  for (cur_mod in mods) {

    cur_hub_df <- GetHubGenes(seurat_obj, n_hubs = n_hubs, wgcna_name = wgcna_name) %>%
      subset(module == cur_mod)

    rownames(cur_hub_df) <- cur_hub_df$gene_name
    cur_hub_genes <- cur_hub_df$gene_name

    cur_TOM <- TOM[cur_hub_genes, cur_hub_genes]
    cur_TOM[lower.tri(cur_TOM)] <- 0

    cur_TOM_df <- reshape2::melt(cur_TOM) %>%
      rename(gene1 = Var1, gene2 = Var2, weight = value) %>%
      subset(gene1 != gene2) %>%
      arrange(desc(weight)) %>%
      slice_max(order_by = weight, n = n_edges) %>%
      mutate(module = cur_mod)

    TOM_df <- rbind(TOM_df, cur_TOM_df)
    hub_df <- rbind(hub_df, cur_hub_df)
  }

  # ---- Construct igraph ----
  graph <- igraph::graph_from_data_frame(TOM_df, directed = FALSE)

  # Add attributes
  igraph::V(graph)$kME <- hub_df[V(graph)$name, "kME"]
  igraph::V(graph)$module <- modules[V(graph)$name, "module"]

  # ---- Return everything ----
  return(list(
    graph     = graph,
    modules   = modules,
    mod_colors = mod_cp,
    TOM_df    = TOM_df,
    hub_df    = hub_df
  ))
}


ModuleNetworkPlotly <- function(
  seurat_obj,
  n_hubs = 25,
  edge_multiplier = 10,
  min_alpha = 0.2,
  max_alpha = 0.8,
  label_zoom_fraction = 0.25,
  plot_width = 1200,
  plot_height = 900,
  module_label_size = 18,
  gene_label_size = 12,
    wgcna_name = NULL
) {

  require(plotly)
  require(dplyr)
  require(scales)
  require(glue)
  require(htmlwidgets)
  require(igraph)

  # --- NEW: use helper ---
  setup <- SetupModuleNetwork(
    seurat_obj = seurat_obj,
    n_hubs = n_hubs,
    edge_multiplier = edge_multiplier,
    wgcna_name = wgcna_name
  )

  graph <- setup$graph
  modules <- setup$modules
  mod_cp <- setup$mod_colors
  TOM_df <- setup$TOM_df
  hub_df <- setup$hub_df

  # Assign node attributes
  igraph::V(graph)$kME <- hub_df[V(graph)$name, "kME"]
  igraph::V(graph)$module <- modules[V(graph)$name, "module"]

  # Graph layout
  coords <- layout_with_graphopt(graph)
  coords <- as.data.frame(coords)
  names(coords) <- c("x", "y")

  # Node dataframe
  nodes <- data.frame(
    gene = igraph::V(graph)$name,
    x = coords$x,
    y = coords$y,
    kME = igraph::V(graph)$kME,
    module = igraph::V(graph)$module
  )

  nodes$color <- mod_cp[nodes$module]


  # calculate centroid coordinates for each module
    module_centroids <- nodes %>%
    group_by(module) %>%
    summarize(
        x = median(x),
        y = median(y)
    )

  # Edge dataframe
  edges <- igraph::as_data_frame(graph, "edges")

  # Scale weights -> alpha
  edges_plot <- edges %>%
    mutate(
      x0 = nodes$x[match(from, nodes$gene)],
      y0 = nodes$y[match(from, nodes$gene)],
      x1 = nodes$x[match(to, nodes$gene)],
      y1 = nodes$y[match(to, nodes$gene)],
      raw_alpha = scales::rescale(weight),
      alpha = pmin(pmax(raw_alpha, min_alpha), max_alpha),
      module = nodes$module[match(from, nodes$gene)],
      edge_color = mod_cp[module]
    )

  # Helper: convert hex color + alpha → rgba
  hex_to_rgba <- function(hex, alpha) {
    rgb <- col2rgb(hex) / 255
    sprintf("rgba(%d,%d,%d,%.3f)",
            round(rgb[1] * 255),
            round(rgb[2] * 255),
            round(rgb[3] * 255),
            alpha)
  }

  edges_plot$rgba <- mapply(hex_to_rgba, edges_plot$edge_color, edges_plot$alpha)

# # --- FIX: customdata must be list-of-lists with one entry per node ---
# nodes$customdata <- lapply(seq_len(nrow(nodes)), function(i) {
#   list(module = nodes$module[i])
# })

# Base plot
p <- plot_ly(
  width = plot_width,
  height = plot_height
)

# Add edges grouped by module
edge_groups <- split(edges_plot, edges_plot$module)

for (m in names(edge_groups)) {
  df <- edge_groups[[m]]

  p <- add_segments(
    p,
    data = df,
    x = ~x0, y = ~y0,
    xend = ~x1, yend = ~y1,
    line = list(color = df$rgba[1], width = 0.6),
    hoverinfo = "none",
    showlegend = FALSE
  )
}

# Add nodes
p <- add_markers(
  p,
  data = nodes,
  x = ~x, y = ~y,
  color = ~module,
  colors = mod_cp,
  marker = list(
    size = ~kME * 15,
    line = list(color = "black", width = 0.8)
  ),
  hoverinfo = "text",
  text = ~paste0(
    "<b>", gene, "</b><br>",
    "Module: ", module, "<br>",
    "kME: ", round(kME, 3)
  )
)

# Label zoom threshold
x_range <- diff(range(nodes$x))
threshold_val <- x_range * label_zoom_fraction

# ---- NEW: Module-level centroids ----
module_centroids <- nodes %>%
  group_by(module) %>%
  summarize(
    x = median(x),
    y = median(y)
  )

# ---- Gene labels (hidden initially) ----
p <- add_text(
  p,
  data = nodes,
  x = ~x, y = ~y,
  text = ~gene,
  textposition = "top center",
  textfont = list(size = gene_label_size, family = "Arial", color="black", face='italic'),
  hoverinfo = "none",
  visible = FALSE,     # hidden initially
  name = "labels"
)

# ---- Module labels (visible initially) ----
p <- add_text(
  p,
  data = module_centroids,
  x = ~x, y = ~y,
  text = ~module,
  textposition = "middle center",
  textfont = list(size = module_label_size, color = "black", family = "Arial Black"),
  hoverinfo = "none",
  visible = TRUE,      # shown when zoomed out
  name = "module_labels"
)

# Base layout
p <- layout(
  p,
  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title=""),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title="")
)

# ---- JS dual-visibility zoom logic ----
js <- glue("
  function(el, x) {{
    var plot = document.getElementById(el.id);

    function compute_xrange(eventdata) {{
      if (eventdata && eventdata['xaxis.range[0]'] !== undefined)
        return eventdata['xaxis.range[1]'] - eventdata['xaxis.range[0]'];
      if (plot.layout.xaxis.range)
        return plot.layout.xaxis.range[1] - plot.layout.xaxis.range[0];
      return null;
    }}

    var threshold = {threshold_val};

    function updateVisibility(eventdata) {{
      var xr = compute_xrange(eventdata);
      if (xr === null) return;

      var showGenes = xr < threshold;
      var showModules = xr >= threshold;

      var geneIndex = plot.data.findIndex(t => t && t.name === 'labels');
      var moduleIndex = plot.data.findIndex(t => t && t.name === 'module_labels');

      if (geneIndex >= 0)
        Plotly.restyle(plot, {{visible: showGenes}}, [geneIndex]);

      if (moduleIndex >= 0)
        Plotly.restyle(plot, {{visible: showModules}}, [moduleIndex]);
    }}

    plot.on('plotly_relayout', function(e) {{ updateVisibility(e); }});
    updateVisibility(null);
  }}
")

# Attach widget JS
p <- htmlwidgets::onRender(p, js)
p

}
