#' @title Plotly Log2FC Network Plot
#'
#' @description This function returns a list with interactive 
#' networkplot based on log2 fold change data.
#' 
#' @param metalyzer_se A MetAlyzer Object
#' @param q_value A numeric value specifying the cutoff for q-value
#' @param metabolite_node_size The text size of the metabolite Nodes
#' @param connection_width The line width of connections between metabolites
#' @param pathway_text_size The text size of pathway annotations
#' @param pathway_width The line width of pathway-specific connection coloring
#' @param plot_height The height of the Plot in pixel [px]
#' @param plot_column_name Column name in the Log2FC dataframe to plot; 
#' for multiple metabolites per node, the mean is used.
#' @return plotly object
#' 
#' @import dplyr
#' @import viridis
#' @import viridisLite
#' @import SummarizedExperiment
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom plotly plot_ly add_annotations layout
#' @export
#' 
#' @examples
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_mutation_data_xl())
#' metalyzer_se <- filterMetabolites(
#'   metalyzer_se,
#'   drop_metabolites = "Metabolism Indicators"
#' )
#' metalyzer_se <- renameMetaData(
#'   metalyzer_se,
#'   Mutant_Control = "Sample Description"
#' )
#' 
#' metalyzer_se <- calculate_log2FC(
#'   metalyzer_se,
#'   categorical = "Mutant_Control",
#'   impute_perc_of_min = 0.2,
#'   impute_NA = FALSE
#' )
#' 
#' p_network <- plotly_network(metalyzer_se, q_value = 0.05)
plotly_network <- function(metalyzer_se,
    q_value=0.05,
    metabolite_node_size=11,
    connection_width=1.25,
    pathway_text_size=20,
    pathway_width=10,
    plot_height=800,
    plot_column_name="log2FC") {
    log2FC_df <- metalyzer_se@metadata$log2FC 
    
    pathway_file <- MetAlyzer::pathway()
    ## Read network nodes, edges and annotations
    pathways <- read_named_region(pathway_file, "Pathways_Header")
    invalid_annotations <- which(
        is.na(pathways$Label) |
        duplicated(pathways$Label) |
        is.na(pathways$x) |
        is.na(pathways$y) |
        is.na(pathways$Color)
    )
    if (length(invalid_annotations) > 0) {
        # print warning and remove
        cat("Warning: Removing", length(invalid_annotations), "invalid pathways.\n")
        pathways <- pathways[-invalid_annotations, ]
    }
    rownames(pathways) <- pathways$Label
    
    nodes <- read_named_region(pathway_file, "Metabolites_Header")
    nodes$Pathway[is.na(nodes$Pathway)] <- ""
    invalid_nodes <- which(
        is.na(nodes$Label) |
        duplicated(nodes$Label) |
        is.na(nodes$x) |
        is.na(nodes$y) |
        !nodes$Pathway %in% c(rownames(pathways), "")
    )
    if (length(invalid_nodes) > 0) {
        # print warning and remove
        cat("Warning: Removing", length(invalid_nodes), "invalid nodes.\n")
        nodes <- nodes[-invalid_nodes, ]
    }
    rownames(nodes) <- nodes$Label
    # Remove #1 at the end
    nodes$Label <- gsub("#[0-9]+", "", nodes$Label)
    
    edges <- read_named_region(pathway_file, "Connections_Header")
    invalid_edges <- which(
        !edges$Node1 %in% rownames(nodes) |
        !edges$Node2 %in% rownames(nodes) |
        edges$Node1 == edges$Node2
    )
    if (length(invalid_edges) > 0) {
        # print warning and remove
        cat("Warning: Removing", length(invalid_edges), "invalid connections.\n")
        edges <- edges[-invalid_edges, ]
    }
    
    edges$x_start <- nodes[edges$Node1, "x"]
    edges$y_start <- nodes[edges$Node1, "y"]
    edges$x_end <- nodes[edges$Node2, "x"]
    edges$y_end <- nodes[edges$Node2, "y"]
    edges$Color <- sapply(rownames(edges), function(rowname) {
        from <- edges[rowname, "Node1"]
        to <- edges[rowname, "Node2"]
        from_pathway <- nodes[from, "Pathway"]
        to_pathway <- nodes[to, "Pathway"]
        color <- NA
        if (from_pathway == to_pathway & from_pathway != "") {
        color <- pathways[from_pathway, "Color"]
        }
        return(color)
    })
    
    ## Add log2FC to nodes_df
    signif_df <- filter(log2FC_df,
                        !is.na(.data$log2FC),
                        !is.na(.data$qval),
                        !is.na(.data$pval),
                        .data$qval <= q_value)

    nodes$log2FC <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
        tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
        if (nrow(tmp_df) > 0) {
        # Alteast 1 significantly changed
        l2fc <- sum(tmp_df$log2FC) / nrow(tmp_df)
        } else if (nrow(tmp_df) == 0 && any(m_vec %in% log2FC_df$Metabolite)) {
        # Not significantly changed but measured
        l2fc <- sum(log2FC_df$log2FC[which(log2FC_df$Metabolite %in% m_vec)]) / length(m_vec)
        } else {
        # Not measured
        l2fc <- NA
        }
        return(l2fc)
    })
    
    ## Add p-value to nodes_df
    nodes$qval <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
        tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
        if (nrow(tmp_df) > 0) {
        # Alteast 1 significantly changed
        qval <- sum(tmp_df$qval) / nrow(tmp_df)
        } else if (nrow(tmp_df) == 0 && any(m_vec %in% log2FC_df$Metabolite)) {
        # Not significantly changed but measured
        qval <- sum(log2FC_df$qval[which(log2FC_df$Metabolite %in% m_vec)]) / length(m_vec)
        } else {
        # Not measured
        qval <- NA
        }
        return(qval)
    })

    ## Add p-value to nodes_df
    nodes$pval <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
        tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
        if (nrow(tmp_df) > 0) {
        # Alteast 1 significantly changed
        pval <- sum(tmp_df$pval) / nrow(tmp_df)
        } else if (nrow(tmp_df) == 0 && any(m_vec %in% log2FC_df$Metabolite)) {
        # Not significantly changed but measured
        pval <- sum(log2FC_df$pval[which(log2FC_df$Metabolite %in% m_vec)]) / length(m_vec)
        } else {
        # Not measured
        pval <- NA
        }
        return(pval)
    })

    nodes$add_col <- sapply(strsplit(nodes$Metabolites, ";"), function(m_vec) {
        tmp_df <- filter(log2FC_df, .data$Metabolite %in% m_vec)
        if (nrow(tmp_df) > 0) {
        # Alteast 1 value per Node
        value <- sum(tmp_df[plot_column_name]) / nrow(tmp_df)
        } else {
        # Not measured
        value <- NA
        }
        return(value)
    })

    # In case the plotted column is log2FC or qval, as they are calculated differently
    if (plot_column_name == "log2FC") {
        nodes$add_col = nodes$log2FC
    } else if (plot_column_name == "qval") {
        nodes$add_col = nodes$qval
    } else if (plot_column_name == "pval") {
        nodes$add_col = nodes$pval
    }

    ## Draw network
    # Create a plot of the network using ggplotly

    # Preparing Hexcodes for Annotation Colors
    nodes$color <- sapply(nodes$add_col, function(value) {
        if (is.na(value)) {
        return("grey")
        } else {
        # Using the viridis color scale, adjust 'option' based on your preference
        color_scale <- viridis(10, option = "D")
        nodes_range <- na.omit(nodes$add_col)
        
        color_index <- findInterval(value, seq(min(nodes_range), max(nodes_range)+0.1, length.out = length(color_scale) + 1))
        return(color_scale[color_index])
        }
    })

    # Prepare the Edges List
    area_shapes <- list()
    for (i in seq_along(edges$Color)) {
        if (is.na(edges$Color[i])) {
        edges$Color[i] <- "white"   # Assign edges without pathway white
        }
    }
    # Create background of edges
    for (i in 1:nrow(edges)) {
        area_shape <- list(
        type = "line",
        x0 = edges$x_start[i],
        y0 = edges$y_start[i],
        x1 = edges$x_end[i],
        y1 = edges$y_end[i],
        line = list(
            color = edges$Color[i],
            width = pathway_width
        )
        )
        area_shapes[[i]] <- area_shape
    }
    # Create the edges
    edge_shapes <- list()
    for (i in 1:nrow(edges)) {
        edge_shape <- list(
        type = "line",
        x0 = edges$x_start[i],
        y0 = edges$y_start[i],
        x1 = edges$x_end[i],
        y1 = edges$y_end[i],
        line = list(
            color = "grey",
            width = connection_width
        )
        )
        edge_shapes[[i]] <- edge_shape
    }
    edges_area_combined <- c(area_shapes, edge_shapes)
    
    # Create the nodes
    network <- plot_ly(nodes,
                        x = nodes$x,
                        y = nodes$y,
                        type = "scatter",
                        mode = "markers",
                        marker = list(color = nodes$add_col, 
                                    colorbar = list(title = ""),
                                    colorscale='Viridis',
                                    showscale = TRUE),
                        height = plot_height)

    # Add the edges
    p_network <- layout(
        network,
        title = 'Log2(FC) with FDR Correction',
        shapes = edges_area_combined,
        xaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
        yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
        hovermode = FALSE) 
    
    # Add annotations over the nodes
    for (i in 1:nrow(nodes)) {
        p_network <- p_network %>% add_annotations(
        text = nodes$Label[i],
        x = nodes$x[i],
        y = nodes$y[i],
        arrowhead = 0,
        font = list(size = metabolite_node_size, color = "white"),
        ax = 0,
        ay = 0,
        bgcolor = nodes$color[i],
        opacity = 1,
        hovertext = paste0("log2 Fold Change: ", round(nodes$log2FC[i], 5),
                         "\nPathway: ", nodes$Pathway[i],
                         "\nadj. p-value: ", nodes$qval[i],
                         "\np-value: ", round(nodes$pval[i], 5),
                         "\n", plot_column_name, ":", nodes$add_col[i])
        )
    }
    
    # Add annotations for the pathways
    for (i in 1:nrow(pathways)) {
        p_network <- p_network %>% add_annotations(
        text = pathways$Pathway[i],
        x = pathways$x[i],
        y = pathways$y[i],
        xref = "x",
        yref = "y",
        showarrow = FALSE,
        font = list(size = pathway_text_size, color = pathways$Color[i])
        )
    }
    
    return(p_network)
}