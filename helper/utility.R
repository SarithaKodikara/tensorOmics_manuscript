#' Generalized Loading Plot Function for Longitudinal Data
#'
#' Creates a combined loading and trajectory plot for tensor-based analysis results.
#' Works with various data structures by specifying column mappings.
#'
#' @param loadings Matrix of loadings from tensor analysis (features × components)
#' @param merged_df Data frame containing the original data with group, time, ID, and feature columns
#' @param comp Integer indicating which component to plot (default: 1)
#' @param top Integer specifying number of top features to display (default: 8)
#' @param negGroup Character vector of group names with negative loadings (for highlighting)
#' @param posGroup Character vector of group names with positive loadings (for highlighting)
#' @param group_col Name of the group column in merged_df (default: "group")
#' @param time_col Name of the time column in merged_df (default: "time")
#' @param id_col Name of the subject/replicate ID column in merged_df (default: "id")
#' @param feature_start_col Index of first feature column in merged_df (default: 7)
#' @param feature_label Name for features in plot (default: "Feature")
#' @param time_label Label to show at time=0 on x-axis (default: "0")
#' @param group_colors Named vector of colors for groups (optional)
#' @param use_letters Logical, if TRUE use letters (A, B, C...) for y-axis labels, 
#'                    if FALSE use full feature names (default: TRUE)
#' @param gap_spacing Numeric value for spacing between features (default: 1.1)
#'
#' @return A ggplot object
#'
#' @examples
#' # Example 1: Microbiome data
#' loadingPlot(loadings = tplsda_res$x_loadings, 
#'             merged_df = merged_df,
#'             comp = 1, 
#'             top = 10,
#'             negGroup = c("PBX"),
#'             posGroup = c("CTR", "FMT"),
#'             group_col = "rGroup",
#'             time_col = "rDay",
#'             id_col = "Participant",
#'             feature_start_col = 7,
#'             feature_label = "OTU",
#'             group_colors = c("CTR" = "blue", "FMT" = "orange", "PBX" = "green"))
#'
#' # Example 2: Metabolite data
#' loadingPlot(loadings = result$x_loadings,
#'             merged_df = metabolite_df,
#'             comp = 1,
#'             top = 8,
#'             negGroup = c("Chloride"),
#'             posGroup = c("Control", "Carbonate"),
#'             group_col = "group",
#'             time_col = "number_of_days",
#'             id_col = "replicate",
#'             feature_start_col = 7,
#'             feature_label = "Metabolite",
#'             time_label = "6",
#'             group_colors = c("Control" = "blue", "Carbonate" = "orange", 
#'                            "Chloride" = "green", "Phosphate" = "hotpink"))
#'
#' # Example 3: Clinical data
#' loadingPlot(loadings = analysis$x_loadings,
#'             merged_df = clinical_df,
#'             comp = 1,
#'             top = 8,
#'             negGroup = c("Placebo"),
#'             posGroup = c("FMT"),
#'             group_col = "Group",
#'             time_col = "Week_mod",
#'             id_col = "SampleID",
#'             feature_start_col = 6,
#'             feature_label = "Feature",
#'             time_label = "6",
#'             group_colors = c("FMT" = "orange", "Placebo" = "blue"))

loadingPlot <- function(loadings, 
                        merged_df,
                        comp = 1,
                        top = 8,
                        negGroup = NULL,
                        posGroup = NULL,
                        group_col = "group",
                        time_col = "time",
                        id_col = "id",
                        feature_start_col = 7,
                        feature_label = "Feature",
                        time_label = "0",
                        group_colors = NULL,
                        use_letters = TRUE,
                        gap_spacing = 1.1) {
  
  # Load required packages
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(scales)
  require(ggtext)
  require(ggnewscale)
  
  # Extract feature columns
  feature_cols <- colnames(merged_df)[feature_start_col:ncol(merged_df)]
  
  # Extract loadings for specified component
  loadings_cmp <- loadings[, comp]
  names(loadings_cmp) <- feature_cols
  
  # Sort by absolute loading values
  sorted_indices <- order(abs(loadings_cmp), decreasing = TRUE)
  sorted_features <- feature_cols[sorted_indices[1:top]]
  sorted_loadings <- loadings_cmp[sorted_features]
  max_loading <- max(abs(sorted_loadings))
  
  # Create loading dataframe
  tmp_df <- data.frame(
    feature = sorted_features,
    loading = sorted_loadings
  )
  tmp_df$feature <- factor(tmp_df$feature, levels = sorted_features)
  
  tmp_df <- tmp_df %>% 
    mutate(
      sign = ifelse(loading >= 0, "positive", "negative"),
      ypos = (as.numeric(feature) - 1) * gap_spacing + 1
    )
  
  # Select relevant columns and calculate group means
  merged_df_sel <- merged_df[, c(group_col, time_col, id_col, sorted_features)]
  
  # Filter out time = 0 if present (for metabolite case)
  if (time_col == "number_of_days" && 0 %in% merged_df_sel[[time_col]]) {
    merged_df_sel <- merged_df_sel %>% filter(.data[[time_col]] != 0)
  }
  
  merged_df_mean <- merged_df_sel %>% 
    summarise(
      .by = c(all_of(group_col), all_of(time_col)),
      across(
        all_of(sorted_features),
        ~ mean(.x, na.rm = TRUE),
        .names = "{.col}"
      )
    )
  
  # Convert to long format
  long_dat <- pivot_longer(
    merged_df_mean,
    cols = -c(all_of(group_col), all_of(time_col)),
    names_to = "feature",
    values_to = "mean_value"
  )
  long_dat$feature <- factor(long_dat$feature, levels = sorted_features)
  
  # Add y-position
  long_dat <- long_dat %>% 
    mutate(ypos = (as.numeric(feature) - 1) * gap_spacing + 1)
  
  y_breaks <- (0:(top - 1)) * gap_spacing + 1
  
  # Add loading information and scale values
  long_dat <- long_dat %>%
    left_join(tmp_df %>% select(feature, loading), by = "feature") %>%
    mutate(
      sign_type = ifelse(loading >= 0, "positive", "negative"),
      alpha_label = case_when(
        sign_type == "positive" & .data[[group_col]] %in% posGroup ~ "High alpha",
        sign_type == "negative" & .data[[group_col]] %in% negGroup ~ "High alpha",
        TRUE ~ "Low alpha"
      )
    ) %>%
    group_by(feature) %>%
    mutate(
      value_scaled = scales::rescale(mean_value, to = c(unique(ypos) - 0.4, unique(ypos) + 0.4)),
      time = rescale(.data[[time_col]], to = c(0, max_loading))
    ) %>%
    ungroup()
  
  # Set default colors if not provided
  if (is.null(group_colors)) {
    unique_groups <- unique(merged_df[[group_col]])
    group_colors <- setNames(
      scales::hue_pal()(length(unique_groups)),
      unique_groups
    )
  }
  
  # Create labels for caption
  if (use_letters) {
    y_axis_labels <- LETTERS[1:top]
    caption_labels <- paste0(LETTERS[1:top], ": ", sorted_features)
  } else {
    y_axis_labels <- sorted_features
    caption_labels <- NULL
  }
  
  long_labels <- paste(caption_labels, collapse = "\n")
  
  # Create plot
  P <- ggplot(long_dat, aes(x = time)) +
    # Loading bars
    geom_rect(
      data = tmp_df,
      aes(
        xmin = pmin(0, -abs(loading) / 3), 
        xmax = pmax(0, -abs(loading)),
        ymin = ypos - 0.4, 
        ymax = ypos + 0.4, 
        fill = sign
      ),
      alpha = 0.8, 
      inherit.aes = FALSE
    ) +
    scale_y_continuous(breaks = y_breaks, labels = y_axis_labels) + 
    scale_fill_manual(
      values = c(positive = "#e39894", negative = "#8fafd9"),
      name = "Loading"
    ) +
    guides(fill = guide_legend(reverse = TRUE)) +
    # New fill scale for ribbons
    ggnewscale::new_scale_fill() +
    geom_ribbon(
      aes(
        ymin = ypos - 0.4, 
        ymax = ypos + 0.4,
        group = interaction(.data[[group_col]], feature), 
        fill = sign_type
      ), 
      alpha = 0.02, 
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c(positive = "red", negative = "blue"),
      name = "Loading"
    ) +
    # Trajectory lines
    geom_line(
      aes(
        y = value_scaled,
        group = interaction(.data[[group_col]], feature),  
        colour = .data[[group_col]],
        linetype = alpha_label,
        alpha = alpha_label
      ),
      linewidth = 1
    ) +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_alpha_manual(
      values = c("High alpha" = 1, "Low alpha" = 0.35), 
      guide = "none"
    ) +
    scale_linetype_manual(
      values = c("High alpha" = 1, "Low alpha" = 1), 
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = waiver(),
      labels = function(b) ifelse(b == 0, time_label, "")
    ) +
    labs(y = "", x = NULL, title = "") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      axis.line.x = element_line(),
      legend.position = "bottom",
      plot.margin = margin(20, 10, 60, 10),
      axis.text = element_text(size = 12)
    ) +
    geom_vline(xintercept = 0, colour = "black", linetype = "dashed")
  
  # Add caption with feature names if using letters
  if (use_letters && !is.null(caption_labels)) {
    P <- P + 
      labs(caption = long_labels) +
      theme(
        plot.caption = element_textbox_simple(
          halign = 0.5,
          hjust = 0,
          size = 10,
          margin = margin(t = 5),
          width = grid::unit(1, "npc"),
          face = "italic"
        )
      )
  }
  
  return(P)
}


# ============================================================================
# USAGE EXAMPLES FOR EACH CASE
# ============================================================================

#' Example 1: Microbiome Data (Case Study 1)
#' 
#' Data structure:
#' - Group column: "rGroup" (CTR, FMT, PBX)
#' - Time column: "rDay"
#' - ID column: "Participant"
#' - Features start at column 7
#'
example_1_microbiome <- function() {
  loadingPlot(
    loadings = tplsda_res$x_loadings, 
    merged_df = merged_df,
    comp = 1, 
    top = 10,
    negGroup = c("PBX"),
    posGroup = c("CTR", "FMT"),
    group_col = "rGroup",
    time_col = "rDay",
    id_col = "Participant",
    feature_start_col = 7,
    feature_label = "OTU",
    time_label = "0",
    use_letters = TRUE,
    group_colors = c("CTR" = "blue", "FMT" = "orange", "PBX" = "green")
  )
}

#' Example 2: Metabolite/Bioreactor Data
#' 
#' Data structure:
#' - Group column: "group" (Control, Carbonate, Chloride, Phosphate)
#' - Time column: "number_of_days"
#' - ID column: "replicate"
#' - Features start at column 7
#' - Note: Filters out time = 0
#'
example_2_metabolite <- function() {
  loadingPlot(
    loadings = result$x_loadings,
    merged_df = metabolite_df,
    comp = 1,
    top = 8,
    negGroup = c("Chloride"),
    posGroup = c("Control", "Carbonate"),
    group_col = "group",
    time_col = "number_of_days",
    id_col = "replicate",
    feature_start_col = 7,
    feature_label = "Metabolite",
    time_label = "6",
    use_letters = FALSE,  # Show full metabolite names
    group_colors = c(
      "Control" = "blue", 
      "Carbonate" = "orange", 
      "Chloride" = "green", 
      "Phosphate" = "hotpink"
    )
  )
}

#' Example 3: Clinical/FMT Trial Data
#' 
#' Data structure:
#' - Group column: "Group" (FMT, Placebo)
#' - Time column: "Week_mod"
#' - ID column: "SampleID"
#' - Features start at column 6
#'
example_3_clinical <- function() {
  loadingPlot(
    loadings = analysis$x_loadings,
    merged_df = clinical_df,
    comp = 1,
    top = 8,
    negGroup = c("Placebo"),
    posGroup = c("FMT"),
    group_col = "Group",
    time_col = "Week_mod",
    id_col = "SampleID",
    feature_start_col = 6,
    feature_label = "Feature",
    time_label = "6",
    use_letters = FALSE,  # Show full feature names
    group_colors = c("FMT" = "orange", "Placebo" = "blue")
  )
}

