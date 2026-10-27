# --- Yossef et al. functions ---
# Sourced after setup/setup.R. Only functions that differ from functions.R, renamed with the "_yossef" suffix.

plot_reactive_top_yossef <- function(
    object,
    pivot_cols = NULL,
    top_n = NULL,
    sample = NULL
) {
  df_long <- object@meta.data %>%
    select(reactive, all_of(pivot_cols), category, cdr3, TCRb, patient) %>%
    filter(patient %in% sample) %>%
    pivot_longer(
      cols = all_of(pivot_cols),
      names_to = "feature",
      values_to = "value"
    )
  
  max_per_clonotype <- df_long %>%
    group_by(feature, cdr3) %>%
    summarise(
      max_value = max(value, na.rm = TRUE),
      reactive_clonotype = any(category == "Reactive", na.rm = TRUE),
      .groups = "drop"
    )
  
  top_info <- max_per_clonotype %>%
    group_by(feature) %>%
    arrange(desc(max_value), cdr3) %>%
    slice_head(n = top_n) %>%
    mutate(is_top = TRUE) %>%
    ungroup()
  
  thresholds <- top_info %>%
    group_by(feature) %>%
    summarise(thresh_value = min(max_value), .groups = "drop") %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  
  df_long <- df_long %>%
    left_join(top_info %>% select(feature, cdr3, is_top),
              by = c("feature", "cdr3")) %>%
    mutate(is_top = ifelse(is.na(is_top), FALSE, TRUE)) %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  
  top_clonotypes_with_key <- df_long %>%
    filter(is_top) %>%
    distinct(feature, cdr3) %>%
    group_by(feature) %>%
    summarise(top_clonotypes = list(cdr3), .groups = "drop") %>%
    { setNames(.$top_clonotypes, .$feature) }
  
  top_clonotype_table <- df_long %>%
    filter(is_top) %>%
    group_by(feature, cdr3, category) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    arrange(feature, desc(n_cells))
  
  medians <- df_long %>%
    group_by(feature, category) %>%
    summarise(median_value = median(value, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(category)) %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  
  max_possible_reactive <- df_long %>%
    filter(category == "Reactive") %>%
    group_by(feature) %>%
    summarise(max_reactive = n_distinct(cdr3), .groups = "drop") %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  
  id_top_status <- df_long %>%
    group_by(feature, cdr3) %>%
    summarise(
      id_is_reactive = any(category == "Reactive"),
      id_is_top = any(is_top),
      .groups = "drop"
    )
  
  reactive_in_top <- id_top_status %>%
    filter(id_is_top & id_is_reactive) %>%
    group_by(feature) %>%
    summarise(reactive_top = n_distinct(cdr3), .groups = "drop") %>%
    right_join(id_top_status %>% distinct(feature), by = "feature") %>%
    replace_na(list(reactive_top = 0)) %>%
    arrange(feature)
  
  n_labels <- reactive_in_top %>%
    full_join(max_possible_reactive, by = "feature") %>%
    mutate(
      reactive_top = ifelse(is.na(reactive_top), 0, reactive_top),
      max_reactive = ifelse(is.na(max_reactive), 0, max_reactive),
      percent = ifelse(max_reactive == 0, 0,
                       round(100 * reactive_top / max_reactive, 1)),
      label = paste0(reactive_top, "/", max_reactive, "\n(", percent, "%)"),
      feature = factor(feature, levels = pivot_cols)
    )
  
  df_long <- df_long %>%
    mutate(has_match = category == "Reactive") %>%
    mutate(color = case_when(
      is_top & has_match ~ "darkred",
      is_top & !has_match ~ "black",
      !is_top & has_match ~ "gray70",
      TRUE ~ NA_character_
    ))
  
  n_levels <- n_distinct(df_long$category)
  shapes <- c(16, 17, 15, 3, 7, 8, 9, 10, 1, 2)[1:n_levels]
  
  p <- ggplot(df_long, aes(x = feature, y = value)) +
    geom_violin(fill = "white", color = "black", alpha = 0.2, linewidth = 0.25) +
    geom_jitter(
      data = df_long %>% filter(category == "Reactive"),
      aes(color = color),
      width = 0.2,
      size = 1,
      alpha = 0.75
    ) +
    geom_jitter(
      data = medians,
      aes(x = feature, y = median_value, shape = category),
      color = "black",
      size = 1,
      alpha = 0.75,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = thresholds,
      aes(
        x = as.numeric(feature) - 0.3,
        xend = as.numeric(feature) + 0.3,
        y = thresh_value,
        yend = thresh_value
      ),
      color = "black",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    geom_text(
      data = n_labels,
      aes(x = feature, y = 0.8, label = label)
    ) +
    scale_color_identity() +
    scale_shape_manual(values = shapes) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    ylim(-0.1, 1) +
    ggtitle(paste0(
      paste(toupper(sample), collapse = "_"),
      "_top_", top_clonotypes))
  
  reactive_id_table <- df_long %>%
    # unimos el threshold por feature
    left_join(thresholds, by = "feature") %>%
    filter(category == "Reactive") %>%
    group_by(feature, cdr3) %>%
    summarise(
      n_cells_top = sum(value >= thresh_value, na.rm = TRUE),
      n_cells_total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      percent_top = round(100 * n_cells_top / n_cells_total, 1)
    ) %>%
    arrange(feature, desc(n_cells_top)) %>%
    filter(n_cells_top > 0)
  
  return(list(
    plot = p,
    top_clonotypes = top_clonotypes_with_key,
    top_clonotype_table = top_clonotype_table,
    reactive_id_table = reactive_id_table
  ))
}

plot_reactive_tcrs_vs_tested_yossef <- function(object,
                                         step = 250,
                                         max_n = 50000,
                                         pivot_cols = NULL,
                                         sample_id = NULL
) {
  library(dplyr)
  library(ggplot2)
  library(tibble)
  
  meta <- object@meta.data
  meta <- meta %>% subset(meta$patient == sample_id)
  
  required_cols <- c("cdr3", "reactive", "category", pivot_cols)
  if (!all(required_cols %in% colnames(meta))) {
    stop("meta.data must contain: ", paste(required_cols, collapse = ", "))
  }
  meta <- meta %>% select(cdr3, category, pivot_cols)
  
  n_clonotypes <- n_distinct(meta$cdr3)
  max_n <- min(max_n, n_clonotypes)
  x_vals <- seq(0, max_n, by = step)
  
  all_clonotypes <- meta %>%
    distinct(cdr3, category) %>%
    slice_sample(n = n_clonotypes)
  
  # Baseline aleatoria
  baseline <- lapply(x_vals, function(n) {
    if (n == 0) {
      return(tibble(
        TCR_tested = 0,
        Reactive_identified = 0,
        Method = "Random"
      ))
    }
    top <- all_clonotypes %>% slice_head(n = n)
    # contamos reactive reactivos
    reactive_count <- meta %>%
      filter(cdr3 %in% top$cdr3, category == "Reactive") %>%
      summarise(n = n_distinct(cdr3)) %>%
      pull(n)
    
    tibble(
      TCR_tested = n,
      Reactive_identified = reactive_count,
      Method = "Random"
    )
  }) %>% bind_rows()
  
  # Métodos con score
  scored_results <- lapply(pivot_cols, function(score_col) {
    df <- meta
    if (!"cell_id" %in% colnames(df)) {
      df <- tibble::rownames_to_column(df, "cell_id")
    }
    
    df <- df %>%
      select(cell_id, cdr3, category, !!sym(score_col)) %>%
      arrange(desc(.data[[score_col]]))
    
    seen_keys <- character(0)
    counted_keys <- integer(0)
    results <- list()
    
    results[[1]] <- tibble(
      TCR_tested = 0,
      Reactive_identified = 0,
      Method = score_col
    )
    
    for (i in seq_len(nrow(df))) {
      this_key <- df$cdr3[i]
      if (!(this_key %in% seen_keys)) {
        seen_keys <- c(seen_keys, this_key)
      }
      
      n_unique_keys <- length(seen_keys)
      
      if (n_unique_keys %in% x_vals && !(n_unique_keys %in% counted_keys)) {
        counted_keys <- c(counted_keys, n_unique_keys)
        
        reactive_count <- df %>%
          filter(cdr3 %in% seen_keys, category == "Reactive") %>%
          summarise(n = n_distinct(cdr3)) %>%
          pull(n)
        
        results[[length(results) + 1]] <- tibble(
          TCR_tested = n_unique_keys,
          Reactive_identified = reactive_count,
          Method = score_col
        )
      }
      
      if (n_unique_keys >= max_n) break
    }
    
    bind_rows(results)
  }) %>% bind_rows()
  
  all_results <- bind_rows(baseline, scored_results)
  
  y_max <- all_results %>%
    summarise(max_reactive = max(Reactive_identified, na.rm = TRUE)) %>%
    pull(max_reactive)
  
  y_max <- ceiling(y_max / 5) * 5
  
  p <- ggplot(all_results, aes(x = TCR_tested, y = Reactive_identified, color = Method)) +
    geom_line(size = 0.5) +
    labs(
      title = sample_id,
      x = "Interrogated TCRs",
      y = "Reactive IDs"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    scale_color_manual(values = c(
      "Random" = "gray40",
      "TRTPBL_UCell" = "purple",
      "TRTPBL_NeoTCRPBL_UCell" = "purple4",
      "Yossef_NeoTCR_PBL_UCell" = "darkblue",
      "iTAG_pre_UCell" = "darkred",
      "iTAG_post1_UCell" = "brown2",
      "iTAG_post2_UCell" = "orange",
      "iTAG_all_UCell" = "green2"
    )) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, by = 5),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(
      limits = c(0, max_n)
    )
  
  return(p)
}

plot_ROC_curves_yossef <- function(object,
                            pivot_cols = NULL,
                            sample_id = NULL
) {
  library(pROC)
  library(dplyr)
  library(ggplot2)
  
  meta <- object@meta.data
  meta <- meta %>% subset(meta$patient == sample_id)
  
  required_cols <- c("cdr3", "reactive", "category", pivot_cols)
  if (!all(required_cols %in% colnames(meta))) {
    stop("Missing required columns in meta.data")
  }
  
  # --- Label each clonotype as reactive if *any* of its reactive are reactive ---
  clonotype_df <- meta %>%
    group_by(cdr3) %>%
    summarise(
      # a clonotype is reactive if any of its reactive are reactive
      is_reactive = as.integer(any(category == "Reactive")),
      # average score across its cells
      across(all_of(pivot_cols), ~mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # --- Compute ROC curves ---
  roc_objs <- lapply(pivot_cols, function(score_col) {
    roc(clonotype_df$is_reactive, clonotype_df[[score_col]], quiet = TRUE)
  })
  names(roc_objs) <- pivot_cols
  
  # --- Build ROC data for plotting ---
  roc_results <- lapply(pivot_cols, function(score_col) {
    coords_df <- coords(roc_objs[[score_col]], "all", ret=c("specificity", "sensitivity"), transpose = FALSE)
    coords_df$score <- score_col
    coords_df$auc <- auc(roc_objs[[score_col]])
    coords_df
  }) %>%
    bind_rows()
  
  # --- Compute AUCs and DeLong test if applicable ---
  auc_table <- roc_results %>%
    group_by(score) %>%
    summarise(AUC = unique(round(auc, 3)))
  
  p_label <- ""
  if (length(pivot_cols) == 2) {
    delong <- roc.test(roc_objs[[pivot_cols[1]]], roc_objs[[pivot_cols[2]]], method = "delong")
    p_label <- paste0("DeLong p = ", signif(delong$p.value, 3))
  }
  
  auc_text <- paste0(auc_table$score, ": AUC = ", auc_table$AUC, collapse = "\n")
  
  # --- Plot ---
  p <- ggplot(roc_results, aes(x = 1 - specificity, y = sensitivity, color = score)) +
    geom_line(size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    annotate("text", x = 0.3, y = 0.1, label = auc_text, hjust = 0, size = 2.8) +
    labs(
      title = sample_id,
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Signature"
    ) +
    theme_minimal(base_size = 13) +
    theme(aspect.ratio = 1, legend.position = "none") +
    scale_color_manual(values = c(
      "Random" = "gray40",
      "TRTPBL_UCell" = "purple",
      "TRTPBL_NeoTCRPBL_UCell" = "purple4",
      "Yossef_NeoTCR_PBL_UCell" = "darkblue",
      "iTAG_pre_UCell" = "darkred",
      "iTAG_post1_UCell" = "brown2",
      "iTAG_post2_UCell" = "orange",
      "iTAG_all_UCell" = "green2"
    ))
  
  return(p)
}

