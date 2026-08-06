#' Plot a taxon-pathway interaction Sankey diagram
#'
#' Renders a Sankey diagram of confounder-adjusted taxon-pathway
#' correlations for a single disease, as computed by
#' \code{\link{compute_residualCLR_correlations}}. Positive and negative
#' correlations are rendered as separate diagrams by default (matching
#' the two-panel "Figure A/B" convention used throughout this package's
#' companion analyses), since mixing directions in one diagram makes the
#' biological interpretation ambiguous.
#'
#' @param interactions A \code{micara_interactions} object from
#'   \code{\link{compute_residualCLR_correlations}}.
#' @param disease Name of the disease to plot (must be present in
#'   \code{names(interactions)}).
#' @param direction One of \code{"positive"}, \code{"negative"}, or
#'   \code{"both"} (default \code{"positive"}). \code{"both"} renders a
#'   single combined diagram — use with caution, as edge sign is then
#'   only visible via hover text, not colour.
#' @param top_n_per_pathway Optional integer; if supplied, keeps only
#'   the \code{top_n_per_pathway} strongest links (by \code{abs(rho)})
#'   per pathway node, for readability in densely-connected diagrams.
#'   Default \code{NULL} (no trimming — all links passing the cutoffs in
#'   \code{compute_residualCLR_correlations()} are shown).
#' @param node_colors Named character vector mapping node group to
#'   colour. Default \code{c(up = "#2196F3", down = "#F44336", unchanged =
#'   "grey70")} — "up"/"down" refer to the ANCOM-BC2 direction of that
#'   taxon/pathway (enriched/upregulated vs depleted/downregulated),
#'   consistent with \code{node_directions$group} from
#'   \code{compute_residualCLR_correlations()}.
#' @param font_size,node_width Passed to
#'   \code{networkD3::sankeyNetwork()}.
#' @param save_path Optional file path (without extension) to save the
#'   diagram. If supplied, saves both an interactive \code{.html} (via
#'   \code{htmlwidgets::saveWidget()}) and a static \code{.png} (via
#'   \code{webshot2::webshot()} or \code{webshot::webshot()}, whichever
#'   is installed). Default \code{NULL} (no file saved; the widget is
#'   still returned and will display in an interactive session).
#'
#' @return A \code{networkD3} A htmlwidget object (invisibly if
#'   \code{save_path} is supplied).
#'
#' @export
plot_interaction_sankey <- function(
    interactions,
    disease,
    direction          = c("positive", "negative", "both"),
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
  
  #Fallback check for santised disease labels
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
  links        <- disease_data$links

  if (is.null(links) || nrow(links) == 0L) {
    stop(disease, " has no links passing the significance/strength cutoffs; nothing to plot.")
  }

  ## ------------------------------------------------------------------
  ## 1. Filter by direction (explicit sign-based — no ambiguous
  ##    magnitude threshold reused from a prior filtering step)
  ## ------------------------------------------------------------------

  links <- switch(
    direction,
    positive = links[links$rho > 0, , drop = FALSE],
    negative = links[links$rho < 0, , drop = FALSE],
    both     = links
  )

  if (nrow(links) == 0) {
    stop(disease, " has no ", direction, " links passing the cutoffs; nothing to plot.")
  }

  ## ------------------------------------------------------------------
  ## 2. Optional top-N-per-pathway trimming (visual decluttering only —
  ##    does not affect the underlying statistics stored in `interactions`)
  ## ------------------------------------------------------------------

  if (!is.null(top_n_per_pathway)) {
    if (!is.numeric(top_n_per_pathway) || length(top_n_per_pathway) != 1L || top_n_per_pathway < 1) {
      stop("'top_n_per_pathway' must be a positive integer.", call. = FALSE)
    }
    links$abs_rho <- abs(links$rho)
    links_split   <- split(links, links$pathway)
    links_split   <- lapply(links_split, function(df) {
      df <- df[order(-df$abs_rho), , drop = FALSE]
      utils::head(df, top_n_per_pathway)
    })
    links <- do.call(rbind, links_split)
    links$abs_rho <- NULL
    rownames(links) <- NULL
  }

  ## ------------------------------------------------------------------
  ## 3. Build node table (only nodes actually appearing in the filtered
  ##    links — orphan nodes make no sense in a Sankey), joined to
  ##    node_directions for colour grouping
  ## ------------------------------------------------------------------

  node_names <- unique(c(links$taxon, links$pathway))
  
  # Defensive lookup check if node_directions is missing or malformed
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
  
  links$IDsource <- match(links$taxon,   nodes$name) - 1L
  links$IDtarget <- match(links$pathway, nodes$name) - 1L
  links$value    <- abs(links$rho)

  ## ------------------------------------------------------------------
  ## 4. Render
  ## ------------------------------------------------------------------

  color_scale <- htmlwidgets::JS(
    sprintf(
      "d3.scaleOrdinal().domain(%s).range(%s)",
      .to_js_array(names(node_colors)),
      .to_js_array(unname(node_colors))
    )
  )

  widget <- networkD3::sankeyNetwork(
    Links       = as.data.frame(links),
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

  ## ------------------------------------------------------------------
  ## 5. Save output if path specified
  ## ------------------------------------------------------------------
  
  if (!is.null(save_path)) {
    
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
      stop("Package 'htmlwidgets' is required to use save_path but is not installed.", call. = FALSE)
    }
    
    # Strip trailing file extensions if provided by user
    save_path <- sub("\\.(html|png)$", "", save_path, ignore.case = TRUE)
    
    html_path <- paste0(save_path, ".html")
    htmlwidgets::saveWidget(widget, html_path, selfcontained = TRUE)
    
    has_webshot2 <- requireNamespace("webshot2", quietly = TRUE)
    has_webshot  <- requireNamespace("webshot",  quietly = TRUE)
    
    if (has_webshot2) {
      webshot2::webshot(html_path, paste0(save_path, ".png"))
    } else if (has_webshot) {
      webshot::webshot(html_path, paste0(save_path, ".png"), delay = 2)
    } else {
      warning(
        "Neither 'webshot2' nor 'webshot' is installed; only the .html file was saved. ",
        "Install one of them to also export a static .png.",
        call. = FALSE
      )
    }
    
    return(invisible(widget))
  }
  
  widget
}

#' @keywords internal
#' @noRd
.to_js_array <- function(x) {
  paste0("[", paste0("\"", x, "\"", collapse = ", "), "]")
}