#' Plot a taxon-pathway interaction Sankey diagram
#'
#' Renders a Sankey diagram of confounder-adjusted taxon-pathway
#' correlations for a single disease, as computed by
#' \code{\link{compute_residualCLR_correlations}}. Positive and negative
#' correlations are rendered as separate diagrams (matching the two-panel
#' "Figure A/B" convention used throughout this package's companion analyses)
#' to preserve clear biological interpretation.
#'
#' @param interactions A \code{micara_interactions} object from
#'   \code{\link{compute_residualCLR_correlations}}.
#' @param disease Name of the disease to plot (must be present in
#'   \code{names(interactions)}).
#' @param direction One of \code{"both"}, \code{"positive"}, or
#'   \code{"negative"} (default \code{"both"}). When \code{"both"}, returns
#'   and/or saves separate diagrams for positive and negative links.
#' @param top_n_per_pathway Optional integer; if supplied, keeps only
#'   the \code{top_n_per_pathway} strongest links (by \code{abs(rho)})
#'   per pathway node, for readability in densely-connected diagrams.
#'   Default \code{NULL} (no trimming — all links passing the cutoffs in
#'   \code{compute_residualCLR_correlations()} are shown).
#' @param node_colors Named character vector mapping node group to
#'   colour. Default \code{c(up = "#2196F3", down = "#F44336", unchanged =
#'   "grey70", stable = "gray70")} — "up"/"down" refer to the ANCOM-BC2 direction
#'   of that taxon/pathway (enriched/upregulated vs depleted/downregulated),
#'   consistent with \code{node_directions$group} from
#'   \code{compute_residualCLR_correlations()}.
#' @param font_size,node_width Passed to
#'   \code{networkD3::sankeyNetwork()}.
#' @param save_path Optional file path (with or without extension) to save the
#'   diagram(s). If supplied, saves interactive \code{.html} (via
#'   \code{htmlwidgets::saveWidget()}) and static \code{.png} (via
#'   \code{webshot2::webshot()} or \code{webshot::webshot()}) files. When
#'   \code{direction = "both"}, suffixes \code{_positive} and \code{_negative}
#'   are automatically appended to the base path. Default \code{NULL}.
#'
#' @return A \code{networkD3} htmlwidget object (for \code{"positive"} or
#'   \code{"negative"}), or a named list of two \code{networkD3} htmlwidget
#'   objects (\code{$positive} and \code{$negative}) when \code{direction = "both"}
#'   (returned invisibly if \code{save_path} is supplied).
#'
#' @export
#' 
#' 
plot_interaction_sankey <- function(
    interactions,
    disease,
    direction          = c("both", "positive", "negative"),
    top_n_per_pathway  = NULL,
    node_colors        = c(up = "#2196F3", down = "#F44336", unchanged = "grey70", stable = "gray70"),
    font_size          = 12,
    node_width         = 35,
    save_path          = NULL
) {
  
  direction <- match.arg(direction)
  
  if (!inherits(interactions, "micara_interactions")) {
    stop("'interactions' must be the output of compute_residualCLR_correlations().")
  }
  
  if (!is.character(disease) || length(disease) != 1L || is.na(disease) || trimws(disease) == "") {
    stop("'disease' must be a single non-empty character string.", call. = FALSE)
  }
  
  if (is.null(names(node_colors)) || any(names(node_colors) == "")) {
    stop("'node_colors' must be a named vector.", call. = FALSE)
  }
  
  # Fallback check for sanitised disease labels
  if (!disease %in% names(interactions)) {
    disease_clean <- make.names(disease)
    if (disease_clean %in% names(interactions)) {
      disease <- disease_clean
    }
  }
  
  if (!disease %in% names(interactions)) {
    stop(
      "'", disease, "' not found in interactions. Available diseases: ",
      paste(names(interactions), collapse = ", ")
    )
  }
  
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    stop("Package 'networkD3' is required for plot_interaction_sankey() but is not installed.")
  }
  
  disease_data <- interactions[[disease]]
  all_links    <- disease_data$links
  
  if (is.null(all_links) || nrow(all_links) == 0L) {
    stop(disease, " has no links passing the significance/strength cutoffs; nothing to plot.")
  }
  
  ## ------------------------------------------------------------------
  ## Inner helper function to construct and optionally save a single widget
  ## ------------------------------------------------------------------
  build_single_sankey <- function(links_df, link_label, file_suffix) {
    if (is.null(links_df) || nrow(links_df) == 0L) {
      warning(disease, " has no ", link_label, " links passing the cutoffs; skipping.")
      return(NULL)
    }
    
    # Trim top-N links per pathway if requested
    if (!is.null(top_n_per_pathway)) {
      if (!is.numeric(top_n_per_pathway) || length(top_n_per_pathway) != 1L || top_n_per_pathway < 1) {
        stop("'top_n_per_pathway' must be a positive integer.", call. = FALSE)
      }
      links_df$abs_rho <- abs(links_df$rho)
      links_split    <- split(links_df, links_df$pathway)
      links_split    <- lapply(links_split, function(df) {
        df <- df[order(-df$abs_rho), , drop = FALSE]
        utils::head(df, top_n_per_pathway)
      })
      links_df         <- do.call(rbind, links_split)
      links_df$abs_rho <- NULL
      rownames(links_df) <- NULL
    }
    
    # Build node data frame (colored by expression/enrichment status)
    node_names <- unique(c(links_df$taxon, links_df$pathway))
    
    group_lookup <- if (!is.null(disease_data$node_directions) &&
                        all(c("name", "group") %in% names(disease_data$node_directions))) {
      stats::setNames(
        disease_data$node_directions$group,
        disease_data$node_directions$name
      )
    } else {
      character(0)
    }
    
    nodes <- data.frame(
      name  = node_names,
      group = ifelse(node_names %in% names(group_lookup), group_lookup[node_names], "stable"),
      stringsAsFactors = FALSE
    )
    
    links_df$IDsource <- match(links_df$taxon,   nodes$name) - 1L
    links_df$IDtarget <- match(links_df$pathway, nodes$name) - 1L
    links_df$value    <- abs(links_df$rho)
    
    color_scale <- htmlwidgets::JS(
      sprintf(
        "d3.scaleOrdinal().domain(%s).range(%s)",
        .to_js_array(names(node_colors)),
        .to_js_array(unname(node_colors))
      )
    )
    
    # Render Sankey widget (Link strips remain default grey)
    widget <- networkD3::sankeyNetwork(
      Links       = as.data.frame(links_df),
      Nodes       = as.data.frame(nodes),
      Source      = "IDsource",
      Target      = "IDtarget",
      Value       = "value",
      NodeID      = "name",
      NodeGroup   = "group",
      colourScale = color_scale,
      fontSize    = font_size,
      nodeWidth   = node_width,
      iterations  = 0
    )
    
    # Save output files if path provided
    if (!is.null(save_path)) {
      if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
        stop("Package 'htmlwidgets' is required to use save_path but is not installed.", call. = FALSE)
      }
      
      base_path  <- sub("\\.(html|png)$", "", save_path, ignore.case = TRUE)
      out_prefix <- if (file_suffix != "") paste0(base_path, "_", file_suffix) else base_path
      
      html_path <- paste0(out_prefix, ".html")
      htmlwidgets::saveWidget(widget, html_path, selfcontained = TRUE)
      
      has_webshot2 <- requireNamespace("webshot2", quietly = TRUE)
      has_webshot  <- requireNamespace("webshot",  quietly = TRUE)
      
      if (has_webshot2) {
        webshot2::webshot(html_path, paste0(out_prefix, ".png"))
      } else if (has_webshot) {
        webshot::webshot(html_path, paste0(out_prefix, ".png"), delay = 2)
      } else {
        warning(
          "Neither 'webshot2' nor 'webshot' is installed; only the .html file was saved.",
          call. = FALSE
        )
      }
    }
    
    return(widget)
  }
  
  ## ------------------------------------------------------------------
  ## Process directions and return plot widget(s)
  ## ------------------------------------------------------------------
  if (direction == "positive") {
    return(build_single_sankey(all_links[all_links$rho > 0, , drop = FALSE], "positive", ""))
    
  } else if (direction == "negative") {
    return(build_single_sankey(all_links[all_links$rho < 0, , drop = FALSE], "negative", ""))
    
  } else { # direction == "both"
    pos_widget <- build_single_sankey(all_links[all_links$rho > 0, , drop = FALSE], "positive", "positive")
    neg_widget <- build_single_sankey(all_links[all_links$rho < 0, , drop = FALSE], "negative", "negative")
    
    plots <- list(positive = pos_widget, negative = neg_widget)
    plots <- plots[!sapply(plots, is.null)]
    
    if (length(plots) == 0L) {
      stop(disease, " has no links passing the cutoffs; nothing to plot.")
    }
    
    if (!is.null(save_path)) {
      return(invisible(plots))
    }
    return(plots)
  }
}

#' @keywords internal
#' @noRd
.to_js_array <- function(x) {
  paste0("[", paste0("\"", x, "\"", collapse = ", "), "]")
}