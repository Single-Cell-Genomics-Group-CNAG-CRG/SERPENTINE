#' VolcanoPlot: Visualize Differential Gene Expression in a Volcano Plot
#'
#' This function creates a volcano plot to visualize differentially expressed genes based on 
#' effect size (log2 fold change) and significance (FDR-adjusted p-values). Genes are colored 
#' based on their statistical significance and fold-change thresholds, with optional gene labels.
#'
#' @param plot_df A data frame containing differential expression results, including effect 
#'    sizes and significance values.
#' @param signif_col Character, column name in `plot_df` indicating significance values 
#'    (default: 'de_fdr').
#' @param signif_thresh Numeric, significance threshold (default: 0.05). Genes with 
#'    FDR-adjusted p-values below this threshold are considered significant.
#' @param effect_col Character, column name in `plot_df` indicating log2 fold change values 
#'    (default: 'log2_de').
#' @param effect_thresh Numeric, log2 fold-change threshold (default: 1.0). Genes with absolute 
#'    effect sizes above this threshold are considered strongly differentially expressed.
#' @param raster_dpi Numeric, resolution for rasterized points to optimize performance 
#'    (default: 300).
#' @param point_size Numeric, size of plotted points (default: 1.5).
#' @param high_color Character, color for significantly upregulated genes (default: 'darkgoldenrod1').
#' @param low_color Character, color for significantly downregulated genes (default: 'dodgerblue').
#' @param insignif_color Character, color for non-significant genes (default: 'grey').
#' @param label_genes Numeric, number of most upregulated and downregulated genes to label (default: 10).
#' @param plot_lines Logical, whether to add reference lines for significance and effect-size 
#'    thresholds (default: TRUE).
#' @param line_color Character, color of threshold reference lines (default: 'lightgrey').
#' @param x_label Expression, x-axis label (default: `bquote("log"[2]~"(Fold Change)")`).
#' @param y_label Expression, y-axis label (default: `bquote("-log"[10]~"(FDR)")`).
#' @param adjust_outliers Numeric or NA, quantile threshold for adjusting extreme significance 
#'    values (default: 0.999, set to NA to disable).
#' 
#' @return A `ggplot2` object representing the volcano plot.
#'
#' @details 
#' - p-values are transformed to `-log10(FDR)`, with a significance threshold indicated by 
#'   `-log10(signif_thresh)`.
#' - Non-significant genes are shown in `insignif_color`, while significant genes are colored 
#'   according to their direction of change.
#' - A subset of genes with the highest and lowest log2 fold changes are labeled.
#' - The function adjusts infinite values in significance (`Inf` p-values are replaced with the 
#'   next highest observed value).
#' - Uses `rasterise()` from `ggrastr` to efficiently handle large datasets.
#' - Dashed reference lines indicate fold-change and significance thresholds.
#' - Annotations in the top corners display the number of significantly up- and downregulated 
#'   genes.
#' 
#' @export
VolcanoPlot <- function(
    plot_df,
    signif_col = "de_fdr",
    signif_thresh = 0.05,
    effect_col = "log2_de",
    effect_thresh = 1.0,
    raster_dpi = 300,
    point_size = 1.5,
    high_color = 'darkgoldenrod1',
    low_color = 'dodgerblue',
    insignif_color = 'grey',
    label_genes = 10,
    label_size = 3,
    label_col = 'gene',
    plot_lines = TRUE,
    line_color = 'lightgrey',
    x_label = bquote("log"[2]~"(Fold Change)"),
    y_label = bquote("-log"[10]~"(FDR)"),
    adjust_outliers = 0.999
){

    # TODO: ChatGPT add checks 

    # -log10 transform the p-vals
    log_transform <- TRUE
    if(log_transform){
        plot_df <- mutate(plot_df, p = -log10(get(signif_col)))
        signif_col <- 'p'
        signif_thresh <- -log10( signif_thresh)
    }

    # fix Inf values (where p-val was 0):
    tmp <- subset(plot_df, get(signif_col) != Inf)
    signif_max <- max(tmp[,signif_col])
    plot_df[,signif_col] <- ifelse(plot_df[,signif_col] == Inf, signif_max, plot_df[,signif_col])

    # adjust outliers 
    if(!is.na(adjust_outliers)){
        new_max <- quantile(plot_df[,signif_col], adjust_outliers)
        plot_df[,signif_col] <- ifelse(plot_df[,signif_col] > new_max, new_max, plot_df[,signif_col])
    }

    # get significant points
    signif_df <- subset(plot_df, get(signif_col) >= signif_thresh & abs(get(effect_col)) >= effect_thresh)
    signif_df$color <- ifelse(signif_df[,effect_col] > 0, high_color, low_color)

    # get non-significant points
    mid_df <- subset(plot_df, get(signif_col) >= signif_thresh & abs(get(effect_col)) < effect_thresh)
    mid_df$color <- ifelse(mid_df[,effect_col] > 0, high_color, low_color)

    # get points to label
    label_df1 <- slice_max(signif_df, n=label_genes, order_by=get(effect_col))
    label_df2 <- slice_min(signif_df, n=label_genes, order_by=get(effect_col))
    label_df <- rbind(label_df1, label_df2)
    label_df$color <- ifelse(label_df[,effect_col] > 0, high_color, low_color)

    # annotations
    up_right <- subset(signif_df, get(effect_col) > 0) %>% nrow
    up_left <- subset(signif_df, get(effect_col) < 0) %>% nrow

    # define the annotations (number of genes above the fc thresh)
    annotations <- data.frame(
        xpos = c(-Inf,Inf),
        ypos =  c(-Inf,-Inf),
        annotateText = c(as.character(up_left),as.character(up_right)),
        hjustvar = c(-1,2),
        vjustvar = c(-1,-1)) #<- adjust
    
    # define the max region for the x-axis
    x_max <- max(abs(range(plot_df[,effect_col]))) + 0.1

    # initialize plot
    p <- plot_df %>% 
        ggplot(aes(x = get(effect_col), y = get(signif_col))) 

    # plot non-signif points
    p <- p + ggrastr::rasterise(geom_point(
        data = subset(plot_df, get(signif_col) < signif_thresh),
        color = insignif_color,
        size = point_size,
        alpha = 0.2
    ), dpi=raster_dpi)

    # plot signif points that didn't meet the FC thresh
    p <- p + ggrastr::rasterise(geom_point(
            data = mid_df,
            color = mid_df$color,
            alpha = 0.2,
            size = point_size
    ), dpi=raster_dpi)

    # plot significant points
    p <- p + ggrastr::rasterise(geom_point(
            data = signif_df,
            color = signif_df$color,
            size = point_size
    ), dpi=raster_dpi) 

    # plot the points that we will label
    # this is the only layer we don't rasterize
    p <- p + geom_point(
            data = label_df,
            fill = label_df$color,
            color = 'black', 
            shape=21,
            size = point_size*2
    ) 

    # plot the labels
    p <- p + ggrepel::geom_text_repel(
            data = label_df,
            label = label_df[,label_col],
            fontface = 'italic',
            max.overlaps = Inf,
            size = label_size
    )

    # plot the dotted lines lines
    if(plot_lines){
        p <- p + 
            geom_hline(yintercept=signif_thresh, linetype='dashed', color=line_color, linewidth=0.5) +
            geom_vline(xintercept = effect_thresh, color=line_color, linewidth=0.5, linetype='dashed') + 
            geom_vline(xintercept = -effect_thresh, color=line_color, linewidth=0.5, linetype='dashed') 
    }

    # plot the annotations
    p <- p + geom_text(
        inherit.aes=FALSE, 
        data=annotations,
        aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)
    ) 

    # plot theme & additional settings
    p <- p + 
        xlim(c(-x_max, x_max)) + 
        xlab(x_label) +
        ylab(y_label) +
        theme(
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            plot.title = element_text(hjust=0.5)
        ) 

    # return the final plot 
    p

}

#' MultiVolcanoPlot: Visualize Differential Gene Expression Effect Sizes Across Clusters
#'
#' This function creates a multi-group volcano-style plot to visualize the effect sizes 
#' (log2 fold changes) of differentially expressed genes across multiple clusters. 
#' Each cluster is displayed as a horizontal strip, with points representing individual genes.
#'
#' @param plot_df A data frame containing differential expression results, including effect 
#'    sizes and significance values.
#' @param signif_col Character, column name in `plot_df` indicating significance values 
#' @param signif_thresh Numeric, significance threshold (default: 0.05). Genes with 
#'    significance values below this threshold are considered differentially expressed.
#' @param effect_col Character, column name in `plot_df` indicating log2 fold change values 
#' @param effect_thresh Numeric, log2 fold-change threshold (default: 1). Genes with absolute 
#'    effect sizes above this threshold are highlighted.
#'    (default: 'de_fdr').
#'    (default: 'log2_de').
#' @param raster_dpi Numeric, resolution for rasterized points (default: 300).
#' @param point_size Numeric, size of plotted points (default: 2).
#' @param high_color Character, color for high effect sizes (default: 'darkorchid3').
#' @param low_color Character, color for low effect sizes (default: 'seagreen').
#' @param mid_color Character, color for neutral values (default: 'whitesmoke').
#' @param plot_lines Logical, whether to add reference lines for significance and effect-size 
#'    thresholds (default: TRUE).
#' @param line_color Character, color of threshold reference lines (default: 'lightgrey').
#' @param x_label Expression, x-axis label (default: `bquote("log"[2]~"(Fold Change)")`).
#' @param adjust_outliers adjust extreme effect sizes to the 0.0001 and 
#'    0.9999 quantiles. Set to NULL if you don't want to use it.
#' 
#' @return A `ggplot2` object representing the multi-group volcano plot.
#'
#' @details 
#' - Significantly upregulated and downregulated genes are colored based on their effect size.
#' - Non-significant genes are plotted in a neutral color (`mid_color`).
#' - The function automatically orders clusters by the number of differentially expressed genes.
#' - Dashed vertical lines indicate the fold-change threshold (`effect_thresh`).
#' - Text annotations show the number of up- and down-regulated genes per cluster.
#' - Uses `rasterise()` from `ggrastr` to improve performance with large datasets.
#' 
#' @export
MultiVolcanoPlot <- function(
    plot_df,
    signif_col = 'de_fdr',
    signif_thresh = 0.05,
    effect_col = 'log2_de',
    effect_thresh = 1,
    raster_dpi = 300,
    point_size = 2,
    high_color = 'darkgoldenrod1',
    low_color = 'dodgerblue',
    mid_color = 'whitesmoke',
    plot_lines = TRUE,
    line_color = 'lightgrey',
    x_label = bquote("log"[2]~"(Fold Change)"),
    adjust_outliers = 0.9999
){
    
    # # calculate the number of up & down DEGs per cluster
    signif_df <- subset(plot_df, get(signif_col) < signif_thresh & abs(get(effect_col)) > effect_thresh)
    n_up <- subset(signif_df, get(effect_col) > 0) %>% .$cluster %>% table
    n_down <- subset(signif_df, get(effect_col) < 0) %>% .$cluster %>% table
    n_de  <- table(signif_df$cluster)

    # order groups by n_DEGs
    ordered_clusters <- names(n_de)[order(n_de)]
    plot_df$cluster <- factor(
        as.character(plot_df$cluster), 
        levels = ordered_clusters
    )

    # adjust outliers
    if(!is.null(adjust_outliers)){
        adjust_top <- adjust_outliers; adjust_bot <- 1 - adjust_outliers
        new_max <- quantile(plot_df[,effect_col], adjust_top)
        new_min <- quantile(plot_df[,effect_col], adjust_bot)
        plot_df[,effect_col] <- ifelse(plot_df[,effect_col] > new_max, new_max, plot_df[,effect_col])
        plot_df[,effect_col] <- ifelse(plot_df[,effect_col] < new_min, new_min, plot_df[,effect_col])
    }

    # create a dataframe to annotate the plot
    n_df <- data.frame(
        cluster = factor(ordered_clusters, levels=ordered_clusters),
        n_up = as.numeric(n_up[ordered_clusters]),
        n_down = as.numeric(n_down[ordered_clusters])
    ) %>% dplyr::arrange(cluster)
    text_max <- max(plot_df[,effect_col]) + 0.25
    text_min <- min(plot_df[,effect_col]) - 0.25

    # initialize the plot
    p <- plot_df %>% subset(get(signif_col) < signif_thresh) %>%
        ggplot(aes(y = cluster, x=get(effect_col), color=get(effect_col)))

    # add the non-significant points
    p <- p +
        ggrastr::rasterise(
            geom_jitter(
                data = subset(plot_df, get(signif_col) > signif_thresh),
                inherit.aes = FALSE,
                aes(y = cluster, x=get(effect_col)),
                alpha = 0.2,
                color=mid_color,
                size=point_size/3, 
                width=0, height=0.25
            ), dpi=raster_dpi
        ) 

    # add the significant points:
    p <- p + 
        ggrastr::rasterise(
            geom_jitter(size=point_size, width=0, height=0.25), 
            dpi=raster_dpi
        )
    
    # add vertical lines
    if(plot_lines){
        p <- p +
            geom_vline(xintercept = effect_thresh, color=line_color, linewidth=0.5, linetype='dashed') + 
            geom_vline(xintercept = -effect_thresh, color=line_color, linewidth=0.5, linetype='dashed')
    }
        
    # add text annotations 
    p <- p +
        geom_text(
            data = n_df,
            inherit.aes=FALSE,
            aes(y=cluster, x=text_max, label=n_up)
        ) +
        geom_text(
            data = n_df,
            inherit.aes=FALSE,
            aes(y=cluster, x=text_min, label=n_down)
        )

    # add thematic elements:
    p <- p +
        scale_color_gradient2(high=high_color, low=low_color, mid=mid_color) +
        ylab('') + 
        xlab(x_label) 

        # add the forest plot and final thematic elements
    p <- p + 
        ggforestplot::geom_stripes(
            aes(y=cluster), inherit.aes=FALSE, 
            data=group_by(plot_df, cluster) %>% slice_max(order_by=get(effect_col), n=1)
        ) +
        theme(
            plot.title = element_text(hjust=0.5)
        ) + labs(color="") 

    p
}

# this function should not do any subsetting etc. 

#' EffectComparisonPlot: Compare Effect Sizes from Two Differential Expression Tests
#'
#' This function plots effect sizes from two differential expression (DE) analyses against each other. 
#' The x-axis represents the effect size from one test, while the y-axis represents the effect size from 
#' another test. This allows for visualizing consistency, directionality, and significance of gene expression 
#' changes across different conditions, methods, or statistical tests.
#'
#' @param plot_df A data frame containing gene-level DE results, including effect sizes (`log2 Fold Change`) 
#'   and statistical significance values for both tests.
#' @param effect_col1 A character string specifying the column name for the x-axis effect size (e.g., `"log2_de"`). 
#'   Default is `"log2_de"`.
#' @param effect_col2 A character string specifying the column name for the y-axis effect size (e.g., `"log2_dv"`). 
#'   Default is `"log2_dv"`.
#' @param signif_col1 A character string specifying the column name for the significance values of test 1 
#'   (e.g., `"de_fdr"`). Default is `"de_fdr"`.
#' @param signif_col2 A character string specifying the column name for the significance values of test 2 
#'   (e.g., `"dv_fdr"`). Default is `"dv_fdr"`.
#' @param signif_thresh A numeric value specifying the significance threshold. Genes with p-values below this 
#'   threshold are considered significant. Default is `0.05`.
#' @param effect_thresh A numeric value specifying the minimum absolute effect size (log2 fold-change) 
#'   for a gene to be considered significant in either test. Default is `1`.
#' @param plot_dpi A numeric value specifying the resolution (dots per inch) for rasterized points. 
#'   Default is `300`.
#' @param point_size A numeric value controlling the size of plotted points. Default is `2`.
#' @param consistent_color A character string specifying the color for genes that are significant in both tests 
#'   with the same direction of effect. Default is `"darkgoldenrod1"`.
#' @param inconsistent_color A character string specifying the color for genes that are significant in both tests 
#'   but with opposite directions of effect. Default is `"dodgerblue"`.
#' @param insignif_color A character string specifying the color for genes that are not significant in either test. 
#'   Default is `"grey"`.
#' @param plot_lines A logical value indicating whether to draw threshold lines on the plot at `± effect_thresh`. 
#'   Default is `TRUE`.
#' @param line_color A character string specifying the color of the threshold lines. Default is `"lightgrey"`.
#' @param x_label A character string or `bquote` expression for labeling the x-axis. Default is 
#'   `bquote("log"[2]~"(Fold Change), Mean")`.
#' @param y_label A character string or `bquote` expression for labeling the y-axis. Default is 
#'   `bquote("log"[2]~"(Fold Change), Variance")`.
#'
#' @return A `ggplot2` object showing a scatter plot of effect sizes from two DE tests, 
#'   highlighting significant genes and directional consistency.
#'
#' @details
#' - The function categorizes genes into four groups: 
#'   - **Both significant:** Genes significant in both tests (`signif_col1` & `signif_col2`).
#'   - **Test 1 only:** Genes significant only in `effect_col1`.
#'   - **Test 2 only:** Genes significant only in `effect_col2`.
#'   - **Neither significant:** Genes that do not pass significance thresholds.
#' - Points are colored based on consistency: 
#'   - **Consistent** (`same direction`) → `consistent_color` (default: `"darkgoldenrod1"`).
#'   - **Inconsistent** (`opposite directions`) → `inconsistent_color` (default: `"dodgerblue"`).
#'   - **Not significant** → `insignif_color` (default: `"grey"`).
#' - Quadrant counts are displayed in the four corners of the plot to indicate the number of genes 
#'   in each effect-size direction.
#' - The `ggrastr::rasterise()` function improves rendering performance when plotting large datasets.
#' - The plot is generated using `ggplot2`, and the axes are scaled symmetrically.
#'
#' @import ggplot2 dplyr ggrastr
#' @export
EffectComparisonPlot <- function(
    plot_df,
    effect_col1 = 'log2_de',
    effect_col2 = 'log2_dv',
    signif_col1 = 'de_fdr',
    signif_col2 = 'dv_fdr',
    signif_thresh = 0.05,
    effect_thresh = 1,
    plot_dpi = 300,
    point_size = 2,
    consistent_color = 'darkgoldenrod1',
    inconsistent_color = 'dodgerblue',
    insignif_color = 'grey',
    plot_lines = TRUE,
    line_color = 'lightgrey',
    x_label = bquote("log"[2]~"(Fold Change), Mean"),
    y_label = bquote("log"[2]~"(Fold Change), Variance")
) {

    # subset the dataframe to get genes that pass the p-val thresh and the
    signif_df <- subset(
        plot_df, abs(get(effect_col1)) > effect_thresh & abs(get(effect_col2)) > effect_thresh & 
        get(signif_col1) < signif_thresh & get(signif_col2) < signif_thresh
    )

    # genes that are signif in group1 or group2
    signif_g1 <- subset(plot_df, get(signif_col1) < signif_thresh & abs(get(effect_col1)) >= effect_thresh) %>% .$gene
    signif_g2 <- subset(plot_df, get(signif_col2) < signif_thresh & abs(get(effect_col2)) >= effect_thresh) %>% .$gene 
    
    # genes that are signif in both groups, either group, neither groups, etc.
    signif_shared <- intersect(signif_g1, signif_g2)
    signif_g1_only <- setdiff(signif_g1, signif_shared)
    signif_g2_only <- setdiff(signif_g2, signif_shared)
    signif_neither <- subset(plot_df,! gene %in% unique(c(signif_g1, signif_g2)) ) %>% .$gene
    signif_either <- union(signif_g1, signif_g2)

    # add the group labels:
    plot_df$gene_group <- ifelse(plot_df$gene %in% signif_shared, 'both', NA)
    plot_df$gene_group <- ifelse(plot_df$gene %in% signif_g1_only, 'g1', plot_df$gene_group)
    plot_df$gene_group <- ifelse(plot_df$gene %in% signif_g2_only, 'g2', plot_df$gene_group)
    plot_df$gene_group <- ifelse(plot_df$gene %in% signif_neither, 'neither', plot_df$gene_group)

    # are the directions consistent or inconsistent?
    plot_df$consistent <- ifelse(sign(plot_df[,effect_col1]) == sign(plot_df[,effect_col2]), 'consistent', 'inconsistent')
    plot_df$consistent <- ifelse(plot_df$gene %in% signif_either, plot_df$consistent, 'not signif')

    # define shapes and colors
    plot_shapes <- c(23, 24, 25, 21); names(plot_shapes) <- c('both', 'g1', 'g2', 'neither')
    plot_colors <- c(insignif_color, consistent_color, inconsistent_color) 
    names(plot_colors) <- c('not signif', 'consistent', 'inconsistent')

    # annotations
    up_right <- subset(signif_df, get(effect_col1) > 0 & get(effect_col2) > 0) %>% nrow
    down_right <- subset(signif_df, get(effect_col1) > 0 & get(effect_col2) < 0) %>% nrow
    up_left <- subset(signif_df, get(effect_col1) < 0 & get(effect_col2) > 0) %>% nrow
    down_left <- subset(signif_df, get(effect_col1) < 0 & get(effect_col2) < 0) %>% nrow

    annotations <- data.frame(
        xpos = c(-Inf,-Inf,Inf,Inf),
        ypos =  c(-Inf, Inf,-Inf,Inf),
        group = c('Consistent', 'Inconsistent', 'Inconsistent', 'Consistent'),
        annotateText = c(as.character(down_left),as.character(up_left), as.character(down_right),as.character(up_right)),
        hjustvar = c(-1,-1,2,2),
        vjustvar = c(-1,2,-1,2)) #<- adjust

    # get plotting limits
    plot_lim <- max(c(abs(plot_df[,effect_col1]), abs(plot_df[,effect_col2])))

    # initialize the plot
    p <- plot_df %>% 
        ggplot(aes(x=get(effect_col1), y=get(effect_col2), fill = consistent, color = consistent))
    
    # plot the insignificant points
    p <- p + ggrastr::rasterise(
            geom_point(
                data = subset(plot_df, gene_group != 'both'),
                aes(shape=gene_group), alpha=0.5, color=insignif_color, fill=insignif_color
            ), dpi=plot_dpi
        ) 

    # plot the significant points
    p <- p + ggrastr::rasterise(
            geom_point(
                data = subset(plot_df, gene_group == 'both'),
                aes(shape=gene_group), alpha=1
            ), dpi=plot_dpi
        ) 

    # plot the lines 
    if(plot_lines){
        p <- p + 
            geom_vline(xintercept = -effect_thresh, color=line_color, linewidth=0.5, linetype='dashed') + 
            geom_vline(xintercept = effect_thresh, color=line_color, linewidth=0.5, linetype='dashed') + 
            geom_hline(yintercept = -effect_thresh, color=line_color, linewidth=0.5, linetype='dashed') +
            geom_hline(yintercept = effect_thresh, color=line_color, linewidth=0.5, linetype='dashed')  
    }
    
    # add text annotations 
    p <- p + geom_text(
        inherit.aes=FALSE, 
        data=annotations,
        aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), 
        color='black'
    ) 

    # add thematic elements 
    p <- p +
        scale_shape_manual(values=plot_shapes, guide="none") + 
        scale_color_manual(values = plot_colors, guide="none") + 
        scale_fill_manual(values = plot_colors, guide="none") + 
        xlim(c(-plot_lim, plot_lim)) + ylim(c(-plot_lim, plot_lim)) +
        xlab(x_label) +
        ylab(y_label) +
        coord_fixed() + 
        theme(
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            plot.title = element_text(hjust=0.5)
        ) 
 
    p

}

# MiloR differential abundance plot 

#' DABoxPlot: Visualize Differential Abundance Results from MiloR
#'
#' This function generates a box plot with overlaid swarm plots to visualize 
#' differential cell type abundance (DA) results from MiloR. Each point represents 
#' a neighborhood, with colors indicating statistical significance and log-fold change.
#'
#' @param plot_df A data frame containing DA results, including log fold-change (`logFC`) 
#'   and significance values (e.g., `SpatialFDR`).
#' @param group_by A character string specifying the column name used to group neighborhoods 
#'   (e.g., cell types or clusters).
#' @param signif_col A character string specifying the column name containing significance values. 
#'   Default is `"SpatialFDR"`.
#' @param signif_thresh A numeric threshold for significance. Points below this threshold are 
#'   considered differentially abundant. Default is `0.05`.
#' @param raster_dpi A numeric value indicating the resolution (dots per inch) for rasterized points. 
#'   Default is `300`.
#' @param point_size A numeric value specifying the size of points in the swarm plot. Default is `1`.
#' @param high_color A character string specifying the color for high positive log fold-change values. 
#'   Default is `"darkgoldenrod1"`.
#' @param low_color A character string specifying the color for high negative log fold-change values. 
#'   Default is `"dodgerblue"`.
#' @param insignif_color A character string specifying the color for non-significant points. 
#'   Default is `"lightgrey"`.
#'
#' @return A `ggplot2` object displaying a box plot with overlaid swarm plots, showing 
#'   log fold-change distributions across groups.
#'
#' @details
#' - Groups are ordered by mean `logFC`, with text annotations showing their values.
#' - Significant points are colored using `scale_color_gradient2()`, with `insignif_color` for 
#'   non-significant values.
#' - The `ggrastr::rasterise()` function is used to improve plotting performance with large datasets.
#' - The x-axis represents log2 fold-change (`logFC`), while the y-axis represents the `group_by` variable.
#'
#' @import ggplot2 dplyr ggrastr ggbeeswarm
#' @export
DABoxPlot <- function(
    plot_df, 
    group_by = 'cluster',
    signif_col = 'SpatialFDR',
    signif_thresh = 0.05,
    raster_dpi = 300,
    point_size = 1,
    order_groups = TRUE,
    high_color = 'darkgoldenrod1',
    low_color = 'dodgerblue',
    insignif_color = 'lightgrey',
    exclude_mixed = TRUE
){  

    if(exclude_mixed){
        plot_df <- plot_df %>%
            subset(get(group_by) != 'Mixed') 
    }

    # order each group by the mean logFC
    order_df <- plot_df %>% 
        group_by(get(group_by)) %>% 
        summarise(mean = mean(logFC)) %>% 
        as.data.frame()
    colnames(order_df)[1] <- group_by

    if(order_groups){
        order_df <- order_df %>%
            dplyr::arrange(desc(mean)) 

        order_df[,group_by] <- factor(
            as.character(order_df[,group_by]),
            levels = as.character(order_df[,group_by])
        )
    }

    # format the number for cleaner plotting 
    order_df$lab <- format(order_df$mean, digits=1)

    plot_df[,group_by] <- factor(
        as.character(plot_df[,group_by]),
        levels = levels(order_df[,group_by])
    )

    plot_range <- max(abs(plot_df$logFC)) + 0.25

    plot_df <- plot_df %>%
        mutate(logFC_color = ifelse(plot_df[,signif_col] <= signif_thresh, logFC, NA)) %>%
        mutate(Nhood=factor(Nhood, levels=unique(Nhood))) 

    # put the NA points on the top so they get plotted first (on bottom)
    plot_df <- rbind(
        subset(plot_df, is.na(logFC_color)),
        subset(plot_df, !is.na(logFC_color))
    )

    # initialize the plot
    p <- plot_df %>%
        ggplot(aes(x = logFC, y = get(group_by))) 

    # plot the points
    p <- p + 
        ggrastr::rasterise(
            ggbeeswarm::geom_quasirandom(
                method = "pseudorandom",
                aes(color=logFC_color),
                size=point_size
            ), 
        dpi=raster_dpi)
    
    # plot the box plot 
    p <- p + geom_boxplot(fill=NA, outlier.shape=NA)

    # modify the plot limits
    p <- p + xlim(-plot_range, plot_range)

    # color scheme:
    p <- p + 
        scale_color_gradient2(
            high=high_color,
            mid=insignif_color, low=low_color, 
            midpoint=0, 
            na.value=insignif_color
        ) 
    
    # set the theme
    p <- p + theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(linewidth=1,color='black', fill=NA),
      panel.grid.major.x = element_line(linewidth=0.5, color='lightgrey'),
      panel.grid.major.y = element_line(linewidth=0.25, color='lightgrey'),
      plot.title = element_text(hjust=0.5),
    ) + guides(color="none") + 
        ylab('') + 
        xlab(bquote("log"[2]~"(Fold Change)"))

    # add text annotations 
    p <- p +
        geom_text(
            data = order_df,
            inherit.aes=FALSE,
            aes(y=get(group_by), x=plot_range - 0.25, 
            label=lab)
        )
    p

}

#' DAReductionPlot: Visualize Differential Abundance on a Dimensionality Reduction Plot
#'
#' This function overlays differential cell type abundance (DA) results onto a 
#' 2D embedding (e.g., UMAP, t-SNE, or PCA) derived from a `MiloR` object. Each 
#' point represents a neighborhood, with colors indicating log fold-change (logFC) 
#' and sizes corresponding to neighborhood size.
#'
#' @param milo_obj A `Milo` object containing the neighborhood graph and dimensionality reduction coordinates.
#' @param da_results A data frame containing DA results, including log fold-change (`logFC`) 
#'   and significance values (e.g., `SpatialFDR`).
#' @param reduction A character string specifying the name of the dimensionality reduction 
#'   to use for plotting (e.g., `"X_UMAP"` for UMAP, `"X_TSNE"` for t-SNE, `"X_PCA"` for PCA). 
#'   Default is `"X_UMAP"`.
#' @param signif_col A character string specifying the column name containing significance values. 
#'   Default is `"SpatialFDR"`.
#' @param signif_thresh A numeric threshold for significance. Points below this threshold are 
#'   considered differentially abundant. Default is `0.05`.
#' @param raster_dpi A numeric value indicating the resolution (dots per inch) for rasterized points. 
#'   Default is `300`.
#' @param point_sizes A numeric vector of length 2 specifying the range of point sizes 
#'   for neighborhoods. Default is `c(0.1, 10)`.
#' @param high_color A character string specifying the color for high positive log fold-change values. 
#'   Default is `"darkgoldenrod1"`.
#' @param low_color A character string specifying the color for high negative log fold-change values. 
#'   Default is `"dodgerblue"`.
#' @param insignif_color A character string specifying the color for non-significant points. 
#'   Default is `"whitesmoke"`.
#' @param centroids Either `FALSE` (default) or a character string specifying a column in `milo_obj`
#'   metadata used to compute centroids for visualization.
#' @param adjust_outliers A numeric value between 0 and 1 specifying the quantile threshold 
#'   for trimming extreme log fold-change values. Default is `0.999`.
#'
#' @return A `ggplot2` object displaying a scatter plot of neighborhoods in the chosen 
#'   dimensionality reduction space, colored by log fold-change and sized by neighborhood size.
#'
#' @details
#' - The function extracts the 2D embedding from `milo_obj` using `reducedDim()`.
#' - The `nhoodGraph()` function is used to extract neighborhood graph metadata.
#' - Points are sized by the number of cells in each neighborhood (`nhood_size`).
#' - Significant points (`SpatialFDR <= signif_thresh`) are placed on top to emphasize them.
#' - The `ggrastr::rasterise()` function improves rendering performance for large datasets.
#' - Axis labels and ticks are removed to focus on the distribution of neighborhoods.
#' - If `centroids` is specified, cluster centroids are computed and displayed with labels.
#' - If `adjust_outliers` is set, log fold-change values are clipped at the specified quantile.
#'
#' @import ggplot2 dplyr ggrastr ggrepel
#' @export
DAReductionPlot <- function(
    milo_obj,
    da_results,
    reduction = 'X_UMAP',
    signif_col = 'SpatialFDR',
    signif_thresh = 0.05,
    raster_dpi = 300,
    point_sizes = c(0.1, 10),
    high_color = 'darkgoldenrod1',
    low_color = 'dodgerblue',
    insignif_color = 'whitesmoke',
    centroids = FALSE,
    adjust_outliers = 0.999
){

    # get the neighborhood graph from the milo object 
    nh_graph <- miloR::nhoodGraph(milo_obj)

    # get the 2D coordinates
    layout <- reducedDim(milo_obj, reduction)[as.numeric(igraph::vertex_attr(nh_graph)$name),] %>% as.data.frame()
    colnames(layout) <- c('x_coord', 'y_coord')
    rownames(layout) <- 1:nrow(layout)

    # combine with DA results:
    plot_df <- cbind(layout, da_results)
    plot_df$nhood_size <- as.numeric(igraph::vertex_attr(nh_graph)$size)

    print(head(plot_df))
    print(sum(is.na(da_results)))

    # adjust outliers?
    if(!is.na(adjust_outliers)){
        new_max <- quantile(plot_df$logFC, adjust_outliers)
        new_min <-  quantile(plot_df$logFC, 1-adjust_outliers)
        plot_df$logFC <- ifelse(plot_df$logFC> new_max, new_max, plot_df$logFC)
        plot_df$logFC <- ifelse(plot_df$logFC < new_min, new_min, plot_df$logFC)
    }

    # arrange the data by effect size, and put the insignificant points on the bottom
    plot_df <- plot_df %>% dplyr::arrange(abs(logFC))
    plot_df <- rbind(
        subset(plot_df, get(signif_col) > signif_thresh),
        subset(plot_df, get(signif_col) <= signif_thresh)
    )

    # initialize the plot, add the points
    p <- plot_df %>% 
        ggplot(aes(x=x_coord, y=y_coord, color=logFC, size=nhood_size)) +
        ggrastr::rasterise(geom_point(
            data = subset(plot_df, get(signif_col) > signif_thresh),
            color = insignif_color
        ), dpi=raster_dpi) +
        ggrastr::rasterise(geom_point(
            data = subset(plot_df, get(signif_col) <= signif_thresh),
        ), dpi=raster_dpi) 
    
    # add the color scale and size scale
    p <- p + 
        scale_color_gradient2(
            high = high_color,
            mid = insignif_color,
            low = low_color
        ) + 
        scale_size(range = point_sizes) 

    if(!(centroids == FALSE)){
        
        reduction_df <- as.data.frame(reducedDims(milo_obj)[[reduction]])
        colnames(reduction_df) <- c('x_coord', 'y_coord')
        reduction_df$cluster <- as.character(colData(milo_obj)[,centroids])

        centroid_df <- reduction_df %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarise(x_coord = mean(x_coord), y_coord = mean(y_coord))
            
        # add the centroids:
        p <- p + geom_point(
            inherit.aes=FALSE,
            aes(x=x_coord, y=y_coord),
            data = centroid_df,
            color = 'black',
            shape = 18, size = 3
        ) + ggrepel::geom_text_repel(
            inherit.aes=FALSE,
            aes(x=x_coord, y=y_coord, label=cluster),
            data = centroid_df,
            color = 'black'
        )

    }

    # add the plot theme
    p <- p + 
        labs(
            color = bquote("log"[2]~"(FC)"),
            size = "Neighborhood\nsize"
        ) + 
        theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust=0.5)
        )

    # return the plot 
    p

}

#' Generate a Proportion Bar Plot for Group Comparisons
#'
#' This function creates a bar plot showing the scaled proportions of a grouping variable 
#' (e.g., clusters) across two conditions. It supports input from Seurat and Milo objects, 
#' allows for optional subsetting, and normalizes proportions based on the total count of 
#' each condition.
#'
#' @param object A Seurat or Milo object containing single-cell metadata.
#' @param group_var Character. The column name specifying the grouping variable (e.g., clusters).
#' @param condition Character. The column name specifying the condition variable (e.g., treatment).
#' @param g1_name Character. The first condition to compare.
#' @param g2_name Character. The second condition to compare.
#' @param replicate_col Character. Column specifying biological replicates (default: "Replicate").
#' @param reps_remove Character vector. Replicates to exclude (default: NA).
#' @param subset_groups Character. Semicolon-separated list of groups to subset (default: "").
#' @param subset_cols Character. Semicolon-separated list of corresponding columns for subsetting (default: "").
#' @param order_groups Logical. Whether to order groups based on g1_name proportions (default: TRUE).
#' @param plot_raw_counts Logical. Do you want to plot the proportioh bar chart, or the raw cell counts per group? (default: TRUE)
#'
#' @return A ggplot2 object displaying the scaled proportions of each group across conditions.
#' 
#' @import ggplot2 dplyr
#' @export
ProportionBarPlot <- function(
    object,
    group_var,
    condition, 
    g1_name,
    g2_name,
    replicate_col = "Replicate",
    reps_remove = NA,
    subset_groups = "",
    subset_cols = "",
    exclude_groups = "",
    exclude_cols = "",
    order_groups = TRUE,
    plot_raw_counts = FALSE
){

    # get the meta-data from the object:
    if(class(object) == "Milo"){
        cur_meta <- as.data.frame(colData(object))
    } else if(class(object) == "Seurat"){
        cur_meta <- as.data.frame(object@meta.data)
    }

    # subset the seurat metadata by cluster
    cur_meta <- cur_meta %>%
        subset(get(condition) %in% c(g1_name, g2_name))

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
    cur_meta <- cur_meta %>% subset(!(get(replicate_col) %in% reps_remove))

    # set up for calculating proportions
    conditions <- c(g1_name, g2_name)
    clusters  <- cur_meta[,group_var] %>% unique %>% as.character
    scale_vec <- table(cur_meta[,condition])

    # loop over each cluster and calculate the proportion for each condition
    proportion_df <- data.frame()
    for(i in 1:length(clusters)){
        cluster_meta <- subset(cur_meta, get(group_var) == clusters[i])
        #cur_df <- as.data.frame(table(cluster_meta[,condition])) %>% dplyr::rename(Count = Freq)

        cur_df <- data.frame(
          Var1 = c(g1_name, g2_name),
          Count = c(
            sum(cluster_meta[,condition] == g1_name),
            sum(cluster_meta[,condition] == g2_name)
          )
        )

        # compute the non-normalized % of cells in each group
        cur_df$Freq <- cur_df$Count / sum(cur_df$Count)

        # scale frequency to the total number of clusters in each Tissue
        cur_df$Scaled <- cur_df$Count / scale_vec
        cur_df$Scaled <- cur_df$Scaled / sum(cur_df$Scaled)

        # add to ongoing proportion df
        cur_df$cluster <- clusters[i]
        proportion_df <- rbind(proportion_df, cur_df)
    }

    proportion_df <- proportion_df %>% dplyr::rename(condition = Var1)

    # order the groups 
    if(order_groups){
        order_clusters <- proportion_df %>% 
            subset(condition == g1_name) %>%
            dplyr::arrange(desc(Scaled)) %>% .$cluster
        proportion_df$cluster <- factor(as.character(proportion_df$cluster), levels=order_clusters)
        proportion_df <- proportion_df %>% dplyr::arrange(cluster)
    }

    proportion_df$hjust_val <- ifelse(proportion_df$Count > 500, 1.1, -0.1)

    # set the color scheme
    cp <- c('darkgoldenrod1', 'dodgerblue')
    names(cp) <- conditions

    if(!plot_raw_counts){

        # initialize the plot
        p <- ggplot(proportion_df, aes(y=Scaled, x=cluster, fill=condition)) +
            geom_bar(stat='identity') +
            geom_hline(yintercept = 0.5, linetype='dashed', color='black', linewidth=0.5)

        # add thematic elements 
        p <- p + scale_y_continuous(expand = c(0,0)) +
            scale_fill_manual(values=cp) +
            ylab("Scaled proportion") +
            theme(
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                axis.text.x = element_text(angle=45, hjust=1),
                plot.title = element_text(hjust=0.5),
                axis.title.x = element_blank(),
                legend.title = element_blank(),
                axis.line.y = element_line(colour = "black"),
                axis.line.x = element_blank()
            )
    } else{

        # order by highest counts:
        total_counts <- proportion_df %>% group_by(cluster) %>% summarise(
            n_cells = sum(Count)
        ) 

        cluster_order <- total_counts %>% arrange(n_cells) %>% .$cluster

        proportion_df$cluster <- factor(
            as.character(proportion_df$cluster),
            levels = cluster_order
        )

        # initialize the plot
        p <- ggplot(proportion_df, aes(x=Count, y=cluster, fill=condition)) +
            geom_bar(position='dodge', stat='identity') +
            geom_text(
            aes(label = Count, hjust=hjust_val),
            position = position_dodge(width=1),
        )

        # add thematic elements 
        p <- p + scale_x_continuous(expand = c(0,0)) +
            scale_fill_manual(values=cp) +
            xlab("Cell counts") + ylab('') +
            theme(
                axis.line.x = element_line(colour='black'),
                axis.line.y = element_line(colour='black'),
                panel.grid.major.x = element_line(linewidth=0.5, color='lightgrey'),
                plot.title = element_text(hjust=0.5),
                legend.position = 'bottom',
                legend.title=element_blank()
            )

    }

    # return the plot 
    p

}

ProportionBarPlotMulti <- function(
    object,
    group_var,           # e.g., "seurat_clusters"
    condition,           # e.g., "Treatment"
    replicate_col = "Replicate",
    reps_remove = NA,
    subset_groups = "",
    subset_cols = "",
    exclude_groups = "",
    exclude_cols = "",
    order_groups = TRUE,
    plot_raw_counts = FALSE,
    custom_colors = NULL  # Optional: pass a vector of colors
){
    library(dplyr)
    library(ggplot2)

    # 1. Extract Meta-data
    if(inherits(object, "Milo")){
        cur_meta <- as.data.frame(colData(object))
    } else if(inherits(object, "Seurat")){
        cur_meta <- as.data.frame(object@meta.data)
    }

    # 2. Subsetting & Exclusions
    if(!any(subset_groups == "")){
        subset_cols_split <- strsplit(subset_cols, ';')[[1]]
        subset_groups_split <- strsplit(subset_groups, ';')[[1]]
        for(j in seq_along(subset_cols_split)){
            cur_groups <- strsplit(subset_groups_split[j], ',')[[1]]
            cur_meta <- cur_meta %>% filter(!!sym(subset_cols_split[j]) %in% cur_groups)
        }
    }

    if(!any(exclude_groups == "")){
        exclude_cols_split <- strsplit(exclude_cols, ';')[[1]]
        exclude_groups_split <- strsplit(exclude_groups, ';')[[1]]
        for(j in seq_along(exclude_cols_split)){
            cur_groups <- strsplit(exclude_groups_split[j], ',')[[1]]
            cur_meta <- cur_meta %>% filter(!(!!sym(exclude_cols_split[j]) %in% cur_groups))
        }
    }

    # Filter replicates
    cur_meta <- cur_meta %>% filter(!(!!sym(replicate_col) %in% reps_remove))

    # 3. Dynamic setup for N conditions
    # Determine the unique conditions present
    conditions <- unique(as.character(cur_meta[[condition]]))
    clusters <- unique(as.character(cur_meta[[group_var]]))
    
    # Scale vec: Total number of cells per condition (for normalization)
    scale_vec <- table(cur_meta[[condition]])

    # 4. Calculate Proportions
    # We use a more 'tidyverse' approach here to handle N groups efficiently
    proportion_df <- cur_meta %>%
        group_by(!!sym(group_var), !!sym(condition)) %>%
        summarise(Count = n(), .groups = 'drop') %>%
        dplyr::rename(cluster = !!sym(group_var), condition = !!sym(condition))

    # Ensure all combinations exist (fill missing with 0)
    proportion_df <- proportion_df %>%
        tidyr::complete(cluster, condition, fill = list(Count = 0))

    # Compute Normalized Metrics
    proportion_df <- proportion_df %>%
        group_by(cluster) %>%
        mutate(
            Freq = Count / sum(Count),
            # Scaled = (Count / total cells in that condition) then re-normalized to 1
            Scaled = Count / as.numeric(scale_vec[condition])
        ) %>%
        mutate(Scaled = Scaled / sum(Scaled)) %>%
        ungroup()

    # 5. Ordering
    if(order_groups){
        # Order by the proportion of the FIRST condition (as a reference)
        ref_cond <- conditions[1]
        cluster_order <- proportion_df %>%
            filter(condition == ref_cond) %>%
            arrange(desc(Scaled)) %>%
            pull(cluster)
    } else{
        cluster_order <- levels(cur_meta[,group_var])
    }
    proportion_df$cluster <- factor(proportion_df$cluster, levels = cluster_order)
    

    # 6. Plotting
    # Set colors: If not provided, use a standard palette
    if(is.null(custom_colors)){
        cp <- scales::hue_pal()(length(conditions))
    } else {
        cp <- custom_colors
    }
    # names(cp) <- conditions

    if(!plot_raw_counts){
        p <- ggplot(proportion_df, aes(y=Scaled, x=cluster, fill=condition)) +
            geom_bar(stat='identity', position='stack') +
            scale_y_continuous(expand = c(0,0)) +
            scale_fill_manual(values=cp) +
            labs(y = "Scaled Proportion", x = NULL) +
            theme_classic() +
            theme(axis.text.x = element_text(angle=45, hjust=1), legend.title = element_blank())
            
        # Add the 50% line only if there are exactly 2 groups (as it makes less sense for 3+)
        if(length(conditions) == 2){
            p <- p + geom_hline(yintercept = 0.5, linetype='dashed')
        }

    } else {
        p <- ggplot(proportion_df, aes(x=Count, y=cluster, fill=condition)) +
            geom_bar(position='dodge', stat='identity') +
            scale_fill_manual(values=cp) +
            labs(x = "Cell Counts", y = NULL) +
            theme_bw() +
            theme(legend.position = 'bottom', legend.title = element_blank())
    }

    return(p)
}

#-----------------------------------------------------------------
# Plot distributions of scores as a functiohn
#-----------------------------------------------------------------

# update this to take one DF for the signatures and one for 
# the meta-data?

# TODO: Create docstring



SignatureDistPlot <- function(
    signature_df,
    meta_df,
    sample_col,
    group_by,
    signature_col = 'module',
    score_col = 'score',
    shape = NA,
    comparisons = NA,
    comparison_method = 'wilcox',
    raster_dpi = 300,
    point_size = 1.5,
    box_alpha = 0.4,
    box_fill = 'grey',
    box_notch = TRUE,
    box_width = 0.5,
    signature_cp = NA,
    ncol = 5
){

    # TODO: checks 

    # simplify the input metadata to just the columns that we need
    meta_df <- meta_df %>%
        dplyr::select(all_of(c(sample_col, group_by))) %>% 
        dplyr::distinct()


    # merge the signature df with the metadata df
    plot_df <- as.data.frame(left_join(
        as.data.frame(signature_df), 
        as.data.frame(meta_df),
        by = sample_col
    ))

    # check that signature_cp is valid 
    signatures <- as.character(unique(plot_df[,signature_col]))

    # check that signature_col is valid

    # check that comparisons are valid based on what is in group_by
    valid_groups <- as.character(unique(plot_df[,group_by]))

    # set up the plot 
    p <- plot_df %>% ggplot(aes(y = get(score_col), x = get(group_by)))

    # plot the points on the bottom:
    p <- p + 
        ggrastr::rasterise(ggbeeswarm::geom_quasirandom(
            aes(color=get(signature_col)),
            method = "pseudorandom",
            size = point_size
        ), dpi=raster_dpi)

    # add the box 
    p <- p + 
        geom_boxplot(
            outlier.shape=NA, 
            alpha = box_alpha,
            width = box_width,
            fill = box_fill,
            notch = box_notch
        ) 

    # color scheme for each signature?
    if(all(signatures %in% names(signature_cp))){
        p <- p + scale_color_manual(values = signature_cp) 
    }

    # compare groups?
    if(all(unlist(comparisons) %in% valid_groups)){
        my_comparisons <- comparisons 
        p <- p + 
            ggpubr::stat_compare_means(
                aes(label = after_stat(p.signif)),
                comparisons = my_comparisons, 
                method = comparison_method
            ) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
    } 

    # set up the theme:
    p <- p + 
        theme(
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(linewidth=1,color='black', fill=NA),
            panel.grid.major.y = element_line(linewidth=0.25, color='lightgrey'),
            plot.title = element_text(hjust=0.5),
            strip.background = element_blank(),
            strip.text = element_text(face='bold')
        ) + Seurat::NoLegend() + Seurat::RotatedAxis() +
        xlab('') + ylab(score_col)  

    # facet 
    patch <- p + 
        facet_wrap(
            ~get(signature_col), ncol=ncol, scales='free'
        )

    patch
} 

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

  print(levels(plot_df[[group.by]]))

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
      print(cur_split)
      cur_df <- plot_df[plot_df[[split.by]] == cur_split,]
      split_df <- plot_df[plot_df[[split.by]] != cur_split,]
      .PlotSingleEmbedding(cur_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, split_df, cur_split, plot_theme, plot_under, plot_ratio=plot_ratio)
    })

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


FeatureEmbedding <- function(
  seurat_obj,
  features,
  raster = TRUE,
  facet_by = NULL,
  slot = 'data',
  assay = NULL,
  reduction = 'umap',
  dpi = 200,
  dpi_scale=0.5,
  ncol = 4,
  combine=TRUE,
  order_points=TRUE, # 'shuffle'
  point_size = 0.5,
  plot_max = 'q100',
  plot_min = 'q0',
  same_range = FALSE,
  colfunc = viridis::inferno,
  rev_colors = TRUE
){

  input_plot_max <- plot_max
  input_plot_min <- plot_min

  plot_list <- list()

  if(same_range){
    tmp <- GetAssayData(seurat_obj, slot=slot, assay=assay)[features,]
    plot_range <- range(tmp)
    print(plot_range)

    if(is.null(plot_max)){plot_max <- plot_range[2]}
    if(is.null(plot_min)){plot_min <- plot_range[1]}
    print(plot_max)
    print(plot_min)
  }


  for(feature in features){

    plot_df <- seurat_obj@meta.data
    plot_df$plot_x_coord <-  seurat_obj@reductions[[reduction]]@cell.embeddings[,1]
    plot_df$plot_y_coord <-  seurat_obj@reductions[[reduction]]@cell.embeddings[,2]


    # check if the feature is in the meta-data
    if(feature %in% rownames(seurat_obj)){
      if(is.null(assay)){assay <- seurat_obj@active.assay}
      plot_df$plotfeature <- GetAssayData(seurat_obj, slot=slot, assay=assay)[feature,]
    } else if(feature %in% colnames(plot_df)){
      if(!is.numeric(plot_df[[feature]])){
        stop("Specified feature is not numeric. Try plotting with VisDimPlot?")
      }
      plot_df$plotfeature <- plot_df[[feature]]
    } else{
      stop("feature not found in rownames(seurat_obj) or in colnames(seurat_obj@meta.data).")
    }

    if(!same_range){
      plot_max <- input_plot_max
      plot_min <- input_plot_min

      plot_range <- range(plot_df$plotfeature)
      if(!is.null(plot_max)){
        if(is.character(plot_max)){
          quant <- as.numeric(gsub('q', '', plot_max)) / 100
          plot_max <- as.numeric(quantile(plot_df$plotfeature, quant))
        }
        plot_range[2] <- plot_max
        print(plot_max)
        plot_df$plotfeature <- ifelse(
          plot_df$plotfeature > plot_max,
          plot_max,
          plot_df$plotfeature
        )
      }

      if(is.character(plot_min)){
        quant <- as.numeric(gsub('q', '', plot_min)) / 100
        plot_min <- as.numeric(quantile(plot_df$plotfeature, quant))
      }
      plot_range[1] <- plot_min
    }
    #
    # shuffle points:
    if(order_points == TRUE){
      plot_df <- plot_df %>% dplyr::arrange(plotfeature)
    } else if(order_points == "shuffle"){
      plot_df <- plot_df[sample(nrow(plot_df)),]
    }


    # initialize gpgplot
    p <- plot_df %>% subset(plotfeature > plot_min) %>%
      ggplot(aes_(x=~plot_x_coord, y=~plot_y_coord,color=~plotfeature))

    # add all the grey dots with low/zero expression
    if(raster){
      p <- p +
        ggrastr::rasterise(
          geom_point(inherit.aes=FALSE, data = subset(plot_df, plotfeature <= plot_min), aes(x=plot_x_coord, y=plot_y_coord),color='lightgrey', size=point_size),
          dpi=dpi, scale=dpi_scale
        )  +
        ggrastr::rasterise(
          geom_point(size=point_size),
          dpi=dpi, scale=dpi_scale
        )
    } else{
      p <- p +
        geom_point(inherit.aes=FALSE, data = subset(plot_df, plotfeature <= plot_min), aes(x=plot_x_coord, y=plot_y_coord),color='lightgrey', size=point_size) +
        geom_point(size=point_size)
    }

    # add extras to plot:
    colors <- colfunc(256)
    if(rev_colors){colors <- rev(colors)}
    p <- p +
      labs(color = feature) +
      scale_color_gradientn(colors=colors, limits = plot_range) +
      coord_equal() +
      theme(
        plot.title = element_text(hjust=0.5, face='plain'),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
      )

    plot_list[[feature]] <- p

  }

  # if there's only one feature plotted
  if(length(plot_list) == 1){
    p <- plot_list[[1]]

    # do we want to facet?
    if(!is.null(facet_by)){
      p <- p + facet_wrap( ~ get(facet_by), ncol=ncol)
    }

    return(p)
  }

  if(combine){
    patch <- wrap_plots(plot_list, ncol=ncol)
    return(patch)
  } else{
    return(plot_list)
  }


}
