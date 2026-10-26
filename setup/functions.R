# --- Custom functions ---
plot_umap_shadow <- function(object, group_var, color_vector, na_value = "none", ncol = 4) {
  plot_data <- FetchData(object, vars = c("umap_1", "umap_2", group_var)) %>%
    mutate(
      !!group_var := as.character(.data[[group_var]]),
      plot_group = ifelse(is.na(.data[[group_var]]) | .data[[group_var]] == na_value, na_value, .data[[group_var]])
    )
  bg_layer <- plot_data %>% filter(plot_group == na_value)
  group_levels <- sort(unique(plot_data$plot_group))
  group_levels <- group_levels[group_levels != na_value]
  color_vector[[na_value]] <- color_vector[[na_value]] %||% "gray85"
  plot_list <- lapply(group_levels, function(g) {
    fg_layer <- plot_data %>% filter(plot_group == g)
    ggplot() +
      geom_point(data = bg_layer, aes(x = umap_1, y = umap_2),
        color = "gray80", size = 0.1, alpha = 0.5) +
      geom_point(data = fg_layer, aes(x = umap_1, y = umap_2, color = plot_group),
        size = 0.5, alpha = 1) +
      scale_color_manual(values = color_vector, guide = "none") +
      ggtitle(g) +
      theme_void() +
      theme(aspect.ratio = 1)
  })
  wrap_plots(plot_list, ncol = ncol)
}

plot_cluster_percentages <- function(object,
  clusters = NULL,
  category = NULL,
  title = NULL,
  colors = NULL,
  horizontal = FALSE,
  limits = NULL,
  exclude = NULL,
  remove = TRUE) {
  data <- FetchData(object, vars = c(clusters, category))
  if (remove) {
    data <- data %>%
      filter(
        !is.na(.data[[clusters]]),
        !is.na(.data[[category]]),
        if (!is.null(exclude)) !.data[[category]] %in% exclude else TRUE
      )
  }
  if (!is.null(colors)) {
    valid_levels <- intersect(names(colors), unique(data[[category]]))
    data[[category]] <- factor(data[[category]], levels = valid_levels)
    colors <- colors[valid_levels]
  } else {
    data[[category]] <- factor(data[[category]])
  }
  counts <- data %>%
    group_by(across(all_of(c(clusters, category)))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(clusters))) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup()
  p <- ggplot(counts, aes_string(x = clusters, y = "percentage", fill = category)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    labs(x = clusters, y = "Percentage", fill = category, title = title) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank()
    )
  if (!is.null(limits)) {
    p <- p + ylim(limits)
  }
  if (horizontal) {
    p <- p + coord_flip()
  }
  q <- ggplot(counts, aes(x = "", y = percentage, fill = .data[[category]])) +
    geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors) +
    labs(title = title) +
    theme_void() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 10)
    )
  if (!is.null(clusters) && clusters %in% colnames(counts)) {
    q <- q + facet_wrap(as.formula(paste("~", clusters)), nrow = 2)
  }
  wrap_plots(p, q, widths = c(2, 1))
}

plot_box_with_stats <- function(seurat_object,
  feature,
  group_by,
  comparisons,
  category_colors
) {
  df <- FetchData(seurat_object, vars = c(feature, group_by))
  colnames(df) <- c("feature_value", "group")
  df$group <- as.factor(df$group)
  df <- df %>%
    group_by(group) %>%
    filter(n() > 1) %>%
    ungroup()
  group_medians <- df %>%
    group_by(group) %>%
    summarize(median_value = median(feature_value, na.rm = TRUE))
  ordered_groups <- group_medians$group[order(group_medians$median_value)]
  df$group <- factor(df$group, levels = ordered_groups)
  valid_comparisons <- lapply(comparisons, function(x) {
    x_in_data <- x[x %in% levels(df$group)]
    if (length(x_in_data) == 2) return(x_in_data) else return(NULL)
  })
  valid_comparisons <- valid_comparisons[!sapply(valid_comparisons, is.null)]
  if (length(valid_comparisons) == 0) valid_comparisons <- NULL
  y_max <- max(df$feature_value, na.rm = TRUE)
  y_positions <- if (!is.null(valid_comparisons)) seq(y_max * 1.05, y_max * 1.55, length.out = length(valid_comparisons)) else NULL
  ggplot(df, aes(x = group, y = feature_value, fill = group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.75) +
    scale_fill_manual(values = scales::alpha(category_colors, 0.75)) +
    stat_compare_means(
      comparisons = valid_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      label.y = y_positions
    ) +
    labs(
      title = paste(feature, "by", group_by),
      x = group_by,
      y = feature
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

run_GSEA <- function(filt.object,
  tags,
  timepoints,
  categories = c("H", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
) {
  data <- subset(object, time %in% timepoints)
  data <- subset(data, category %in% tags)
  markers.genes <- wilcoxauc(X = data, annotations)
  markers.genes$log10p <- -log10(markers.genes$padj)
  markers.genes <- markers.genes %>%
    filter(padj <= 0.05, group == "Reactive") %>%
    arrange(desc(logFC)) %>%
    select(feature, logFC, log10p)
  ranks <- deframe(markers.genes)
  barplot(sort(ranks, decreasing = TRUE), main = "Gene Ranking")
  GSEA <- list()
  for (cat in categories) {
    m_df <- msigdbr(species = "Homo sapiens", category = cat)
    fgsea_sets <- split(m_df$gene_symbol, m_df$gs_name)
    fgseaRes <- fgsea(fgsea_sets, ranks, maxSize = 500)
    GSEA[[cat]] <- fgseaRes[order(fgseaRes$padj, abs(fgseaRes$NES)), ]
  }
  return(GSEA)
}

extract_significant_fgsea_sets <- function(GSEA_results,
  categories = names(GSEA_results),
  padj_cutoff = 0.05
) {
  all_signatures <- list()
  for (cat in categories) {
    gsea_df <- GSEA_results[[cat]]
    sig_paths <- gsea_df$pathway[gsea_df$padj < padj_cutoff]
    if (length(sig_paths) == 0) next
    m_df <- msigdbr::msigdbr(species = "Homo sapiens", category = cat)
    fgsea_sets <- split(m_df$gene_symbol, m_df$gs_name)
    sig_sets <- fgsea_sets[sig_paths]
    all_signatures <- c(all_signatures, sig_sets)
  }
  return(all_signatures)
}

plot_reactive_top <- function(
  object,
  pivot_cols = NULL,
  top_n = NULL,
  time_filter = NULL,
  enrichment_filter = NULL
) {
  df_long <- object@meta.data %>%
    select(reactivity, all_of(pivot_cols), category, key, matched_ID, TCRb, time, enrichment) %>%
    filter(time %in% time_filter & enrichment %in% enrichment_filter) %>%
    pivot_longer(
      cols = all_of(pivot_cols),
      names_to = "feature",
      values_to = "value"
    )
  max_per_clonotype <- df_long %>%
    group_by(feature, key) %>%
    summarise(
      max_value = max(value, na.rm = TRUE),
      reactive_clonotype = any(category == "Reactive", na.rm = TRUE),
      .groups = "drop"
    )
  top_info <- max_per_clonotype %>%
    group_by(feature) %>%
    arrange(desc(max_value), key) %>%
    slice_head(n = top_n) %>%
    mutate(is_top = TRUE) %>%
    ungroup()
  thresholds <- top_info %>%
    group_by(feature) %>%
    summarise(thresh_value = min(max_value), .groups = "drop") %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  df_long <- df_long %>%
    left_join(
      top_info %>% select(feature, key, is_top),
      by = c("feature", "key")
    ) %>%
    mutate(is_top = ifelse(is.na(is_top), FALSE, TRUE)) %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  top_clonotypes_with_key <- df_long %>%
    filter(is_top) %>%
    distinct(feature, key) %>%
    group_by(feature) %>%
    summarise(top_clonotypes = list(key), .groups = "drop") %>%
    { setNames(.$top_clonotypes, .$feature) }
  top_clonotype_table <- df_long %>%
    filter(is_top) %>%
    group_by(feature, key, category) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    arrange(feature, desc(n_cells))
  medians <- df_long %>%
    group_by(feature, category) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(category)) %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  max_possible_reactive <- df_long %>%
    filter(
      category == "Reactive",
      !is.na(matched_ID),
      matched_ID != ""
    ) %>%
    group_by(feature) %>%
    summarise(
      max_reactive = n_distinct(matched_ID),
      .groups = "drop"
    ) %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  matched_ID_top_status <- df_long %>%
    filter(
      !is.na(matched_ID),
      matched_ID != ""
    ) %>%
    group_by(feature, matched_ID) %>%
    summarise(
      matched_ID_is_reactive = any(category == "Reactive"),
      matched_ID_is_top = any(is_top),
      .groups = "drop"
    )
  reactive_in_top <- matched_ID_top_status %>%
    filter(matched_ID_is_top & matched_ID_is_reactive) %>%
    group_by(feature) %>%
    summarise(
      reactive_top = n_distinct(matched_ID),
      .groups = "drop"
    ) %>%
    right_join(
      matched_ID_top_status %>% distinct(feature),
      by = "feature"
    ) %>%
    replace_na(list(reactive_top = 0)) %>%
    arrange(feature)
  n_labels <- reactive_in_top %>%
    full_join(max_possible_reactive, by = "feature") %>%
    mutate(
      reactive_top = ifelse(is.na(reactive_top), 0, reactive_top),
      max_reactive = ifelse(is.na(max_reactive), 0, max_reactive),
      percent = ifelse(
        max_reactive == 0,
        0,
        round(100 * reactive_top / max_reactive, 1)
      ),
      label = paste0(
        reactive_top,
        "/",
        max_reactive,
        "\n(",
        percent,
        "%)"
      ),
      feature = factor(feature, levels = pivot_cols)
    )
  df_long <- df_long %>%
    mutate(has_match = category == "Reactive") %>%
    mutate(
      color = case_when(
        is_top & has_match ~ "darkred",
        is_top & !has_match ~ "black",
        !is_top & has_match ~ "gray70",
        TRUE ~ NA_character_
      )
    )
  n_levels <- n_distinct(df_long$category)
  shapes <- c(16, 17, 15, 3, 7, 8, 9, 10, 1, 2)[1:n_levels]
  p <- ggplot(df_long, aes(x = feature, y = value)) +
    geom_violin(
      fill = "white",
      color = "black",
      alpha = 0.2,
      linewidth = 0.25
    ) +
    geom_jitter(
      data = df_long %>% filter(category == "Reactive"),
      aes(color = color),
      width = 0.2,
      size = 1,
      alpha = 0.75
    ) +
    geom_jitter(
      data = medians,
      aes(
        x = feature,
        y = median_value,
        shape = category
      ),
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
      aes(
        x = feature,
        y = 0.8,
        label = label
      )
    ) +
    scale_color_identity() +
    scale_shape_manual(values = shapes) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
    ) +
    ylim(-0.01, 1) +
    ggtitle(
      paste0(
        paste(toupper(enrichment_filter), collapse = "_"),
        "_",
        paste(toupper(time_filter), collapse = "_"),
        "_top_",
        top_n
      )
    )
  reactive_id_table <- df_long %>%
    left_join(thresholds, by = "feature") %>%
    filter(
      category == "Reactive",
      !is.na(matched_ID),
      matched_ID != ""
    ) %>%
    group_by(feature, matched_ID) %>%
    summarise(
      n_cells_top = sum(value >= thresh_value, na.rm = TRUE),
      n_cells_total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      percent_top = round(
        100 * n_cells_top / n_cells_total,
        1
      )
    ) %>%
    arrange(feature, desc(n_cells_top)) %>%
    filter(n_cells_top > 0)
  return(
    list(
      plot = p,
      top_clonotypes = top_clonotypes_with_key,
      top_clonotype_table = top_clonotype_table,
      reactive_id_table = reactive_id_table
    )
  )
}

plot_reactive_tcrs_vs_tested <- function(object,
  step = 25,
  max_n = 250,
  pivot_cols = NULL,
  title = "Reactivity Enrichment by Signature") {
  meta <- object@meta.data
  required_cols <- c(
    "key",
    "matched_ID",
    "category",
    pivot_cols
  )
  if (!all(required_cols %in% colnames(meta))) {
    stop(
      "meta.data must contain: ",
      paste(required_cols, collapse = ", ")
    )
  }
  n_clonotypes <- n_distinct(meta$key)
  max_n <- min(max_n, n_clonotypes)
  x_vals <- seq(0, max_n, by = step)
  all_clonotypes <- meta %>%
    distinct(key, category) %>%
    slice_sample(n = n_clonotypes)
  baseline <- lapply(x_vals, function(n) {
    if (n == 0) {
      return(
        tibble(
          TCR_tested = 0,
          Reactive_identified = 0,
          Method = "Random"
        )
      )
    }
    top <- all_clonotypes %>%
      slice_head(n = n)
    reactive_count <- meta %>%
      filter(
        key %in% top$key,
        category == "Reactive",
        !is.na(matched_ID),
        matched_ID != ""
      ) %>%
      summarise(
        n = n_distinct(matched_ID)
      ) %>%
      pull(n)
    tibble(
      TCR_tested = n,
      Reactive_identified = reactive_count,
      Method = "Random"
    )
  }) %>%
    bind_rows()
  scored_results <- lapply(pivot_cols, function(score_col) {
    df <- meta
    if (!"cell_id" %in% colnames(df)) {
      df <- tibble::rownames_to_column(
        df,
        "cell_id"
      )
    }
    df <- df %>%
      select(
        cell_id,
        key,
        matched_ID,
        category,
        !!sym(score_col)
      ) %>%
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
      this_key <- df$key[i]
      if (!(this_key %in% seen_keys)) {
        seen_keys <- c(
          seen_keys,
          this_key
        )
      }
      n_unique_keys <- length(seen_keys)
      if (
        n_unique_keys %in% x_vals &&
        !(n_unique_keys %in% counted_keys)
      ) {
        counted_keys <- c(
          counted_keys,
          n_unique_keys
        )
        reactive_count <- df %>%
          filter(
            key %in% seen_keys,
            category == "Reactive",
            !is.na(matched_ID),
            matched_ID != ""
          ) %>%
          summarise(
            n = n_distinct(matched_ID)
          ) %>%
          pull(n)
        results[[length(results) + 1]] <- tibble(
          TCR_tested = n_unique_keys,
          Reactive_identified = reactive_count,
          Method = score_col
        )
      }
      if (n_unique_keys >= max_n) {
        break
      }
    }
    bind_rows(results)
  }) %>%
    bind_rows()
  all_results <- bind_rows(
    baseline,
    scored_results
  )
  y_max <- all_results %>%
    summarise(
      max_reactive = max(
        Reactive_identified,
        na.rm = TRUE
      )
    ) %>%
    pull(max_reactive)
  y_max <- max(5, ceiling(y_max / 5) * 5)
  p <- ggplot(
    all_results,
    aes(
      x = TCR_tested,
      y = Reactive_identified,
      color = Method
    )
  ) +
    geom_line(size = 0.5) +
    labs(
      title = title,
      x = "Interrogated TCRs",
      y = "Reactive matched IDs"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      aspect.ratio = 1,
      legend.position = "right",
      legend.title = element_blank()
    ) +
    scale_color_manual(
      values = signature_colors
    ) +
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

plot_ROC_curves <- function(object,
  pivot_cols = c("Yossef_NeoTCR_PBL_UCell", "iTAG_UCell"),
  title) {
  meta <- object@meta.data
  required_cols <- c(
    "key",
    "matched_ID",
    "category",
    pivot_cols
  )
  if (!all(required_cols %in% colnames(meta))) {
    stop(
      "Missing required columns in meta.data: ",
      paste(
        required_cols[
          !required_cols %in% colnames(meta)
        ],
        collapse = ", "
      )
    )
  }
  empty_roc_plot <- function(message) {
    ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.58,
        label = "ROC unavailable",
        size = 5
      ) +
      annotate(
        "text",
        x = 0.5,
        y = 0.42,
        label = message,
        size = 3.5
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      labs(
        title = title,
        x = NULL,
        y = NULL
      ) +
      theme_void() +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)
      )
  }
  clonotype_df <- meta %>%
    filter(
      !is.na(key),
      key != ""
    ) %>%
    group_by(key) %>%
    summarise(
      is_reactive = as.integer(
        any(
          category == "Reactive" &
          !is.na(matched_ID) &
          matched_ID != ""
        )
      ),
      across(
        all_of(pivot_cols),
        ~mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  response_counts <- clonotype_df %>%
    count(
      is_reactive,
      name = "n_clonotypes"
    )
  if (nrow(clonotype_df) == 0) {
    return(
      empty_roc_plot(
        "No clonotypes are available in this subset"
      )
    )
  }
  if (n_distinct(clonotype_df$is_reactive) < 2) {
    response_label <- paste0(
      ifelse(
        response_counts$is_reactive == 1,
        "Reactive",
        "Non-reactive"
      ),
      ": ",
      response_counts$n_clonotypes,
      collapse = "\n"
    )
    return(
      empty_roc_plot(
        paste0(
          "Only one response class is present\n",
          response_label
        )
      )
    )
  }
  roc_objs <- lapply(pivot_cols, function(score_col) {
    score_df <- clonotype_df %>%
      transmute(
        is_reactive,
        score = .data[[score_col]]
      ) %>%
      filter(
        !is.na(score),
        is.finite(score)
      )
    if (
      nrow(score_df) == 0 ||
      n_distinct(score_df$is_reactive) < 2
    ) {
      return(NULL)
    }
    tryCatch(
      roc(
        score_df$is_reactive,
        score_df$score,
        quiet = TRUE
      ),
      error = function(e) NULL
    )
  })
  names(roc_objs) <- pivot_cols
  roc_objs <- roc_objs[
    !vapply(
      roc_objs,
      is.null,
      logical(1)
    )
  ]
  if (length(roc_objs) == 0) {
    return(
      empty_roc_plot(
        "No signature has valid scores for both response classes"
      )
    )
  }
  roc_results <- lapply(names(roc_objs), function(score_col) {
    coords_df <- coords(
      roc_objs[[score_col]],
      x = "all",
      transpose = FALSE
    )
    coords_df$score <- score_col
    coords_df$auc <- as.numeric(
      auc(roc_objs[[score_col]])
    )
    coords_df
  }) %>%
    bind_rows()
  auc_table <- roc_results %>%
    group_by(score) %>%
    summarise(
      AUC = unique(round(auc, 3)),
      .groups = "drop"
    )
  auc_text <- paste0(
    auc_table$score,
    ": AUC = ",
    auc_table$AUC,
    collapse = "\n"
  )
  p <- ggplot(
    roc_results,
    aes(
      x = 1 - specificity,
      y = sensitivity,
      color = score
    )
  ) +
    geom_line(linewidth = 0.5) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray60"
    ) +
    annotate(
      "text",
      x = 0.3,
      y = 0.1,
      label = auc_text,
      hjust = 0,
      size = 2.8
    ) +
    labs(
      title = title,
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Signature"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      aspect.ratio = 1,
      legend.position = "none"
    ) +
    scale_color_manual(
      values = signature_colors
    )
  return(p)
}

generate_plots <- function(object, time_filter, enrichment_filter, title_suffix) {
  res <- plot_reactive_top(
    object,
    pivot_cols = as.character(signatures_to_compare),
    top_n = top_clonotypes,
    time_filter = time_filter,
    enrichment_filter = enrichment_filter
  )
  a <- res$plot
  print(res$reactive_id_table)
  b <- plot_reactive_tcrs_vs_tested(
    object = object %>%
    subset(
      time == time_filter &
      enrichment %in% enrichment_filter
    ),
    title = paste0(
      time_filter,
      "_",
      paste(enrichment_filter, collapse = "+")
    ),
    step = 25,
    pivot_cols = as.character(signatures_to_compare),
    max_n = length(unique(object$key))
  )
  c <- plot_ROC_curves(
    object = object %>%
    subset(
      time == time_filter &
      enrichment %in% enrichment_filter
    ),
    pivot_cols = as.character(signatures_to_compare),
    title = paste0(
      time_filter,
      "_",
      paste(enrichment_filter, collapse = "+")
    )
  )
  list(
    a = a,
    b = b,
    c = c
  )
}

generate_plots_all <- function(
  object,
  time_filter = NULL,
  enrichment_filter,
  title_suffix
) {
  time_values <- if (is.null(time_filter)) {
    unique(object$time)
  } else {
    time_filter
  }
  object_filtered <- object %>%
    subset(
      time %in% time_values &
      enrichment %in% enrichment_filter
    )
  # Violin plots
  res <- plot_reactive_top(
    object,
    pivot_cols = as.character(signatures_to_compare),
    top_n = top_clonotypes,
    time_filter = time_values,
    enrichment_filter = enrichment_filter
  )
  a <- res$plot +
    ggtitle(
      paste0(
        paste(toupper(enrichment_filter), collapse = "_"),
        "_",
        toupper(title_suffix),
        "_top_",
        top_clonotypes
      )
    )
  print(res$reactive_id_table)
  # Curves
  b <- plot_reactive_tcrs_vs_tested(
    object = object_filtered,
    title = paste0(
      title_suffix,
      "_",
      paste(enrichment_filter, collapse = "+")
    ),
    step = 25,
    pivot_cols = as.character(signatures_to_compare),
    max_n = length(unique(object_filtered$key))
  )
  # AUC
  c <- plot_ROC_curves(
    object = object_filtered,
    pivot_cols = as.character(signatures_to_compare),
    title = paste0(
      title_suffix,
      "_",
      paste(enrichment_filter, collapse = "+")
    )
  )
  list(
    a = a,
    b = b,
    c = c
  )
}

radar_compare_reactive <- function(
  object,
  gene_list,
  reactive_color = "darkred",
  other_color = "gray15",
  scale = c("linear", "log2", "zscore"),
  filter_up = FALSE,
  filter_direction = c("Reactive", "Other")
) {
  scale <- match.arg(scale)
  filter_direction <- match.arg(filter_direction)
  genes <- gene_list
  expr <- FetchData(object, vars = genes)
  expr$Group <- ifelse(object$category == "Reactive", "Reactive", "Other")
  expr$Group[is.na(expr$Group)] <- "Other"
  df <- expr %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "Gene",
      values_to = "Value"
    )
  df <- df %>%
    group_by(Gene) %>%
    mutate(
      Value_scaled = case_when(
        scale == "linear" ~ Value,
        scale == "log2"   ~ log2(Value + 1),
        scale == "zscore" ~ (Value - mean(Value, na.rm = TRUE)) / sd(Value, na.rm = TRUE)
      )
    ) %>%
    ungroup()
  data_summary <- df %>%
    group_by(Group, Gene) %>%
    summarise(
      Mean = mean(Value_scaled, na.rm = TRUE),
      SEM  = sd(Value_scaled, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    filter(Gene %in% genes)
  filtered_genes <- genes
  if (filter_up) {
    data_summary <- data_summary %>%
      group_by(Gene) %>%
      mutate(
        mean_R = Mean[Group == "Reactive"],
        mean_O = Mean[Group == "Other"],
        sem_R  = SEM[Group == "Reactive"],
        sem_O  = SEM[Group == "Other"]
      ) %>%
      filter(
        if (filter_direction == "Reactive") {
      mean_R > mean_O & (mean_R - mean_O) > 1.5 * sem_R
        } else {
      mean_O > mean_R & (mean_O - mean_R) > 1.5 * sem_R
        }
      ) %>%
      ungroup() %>%
      select(-mean_R, -mean_O, -sem_R, -sem_O)
    filtered_genes <- unique(data_summary$Gene)
  }
  data_wide <- data_summary %>%
    pivot_wider(
      names_from = Group,
      values_from = Mean,
      values_fill = 0
    )
  if (filter_up) {
    if (filter_direction == "Other") {
      data_wide <- data_wide %>% arrange(desc(Other))
    } else {
      data_wide <- data_wide %>% arrange(desc(Reactive))
    }
  } else {
    data_wide <- data_wide %>%
      mutate(diff = Reactive - Other) %>%
      arrange(desc(diff))
  }
  genes <- unique(data_wide$Gene)
  n <- length(genes)
  angles <- (pi/2) - 2*pi * ((0:(n-1)) / n)
  angle_map <- tibble(
    Gene = genes,
    angle = angles
  ) %>%
    mutate(
      angle_deg_raw = angle * 180 / pi,
      text_angle = ifelse(angle_deg_raw < -90 | angle_deg_raw > 90, angle_deg_raw + 180, angle_deg_raw),
      hjust = ifelse(angle_deg_raw < -90 | angle_deg_raw > 90, 1, 0)
    )
  data_plot <- data_summary %>%
    left_join(angle_map %>% select(Gene, angle), by = "Gene") %>%
    group_by(Group) %>%
    arrange(angle, .by_group = TRUE) %>%
    mutate(
      x = Mean * cos(angle),
      y = Mean * sin(angle),
      r_top = Mean + SEM,
      r_bottom = pmax(Mean - SEM, 0),
      x_top = r_top * cos(angle),
      y_top = r_top * sin(angle),
      x_bottom = r_bottom * cos(angle),
      y_bottom = r_bottom * sin(angle)
    ) %>%
    bind_rows(slice(., 1)) %>%
    ungroup()
  max_r <- max(data_summary$Mean + data_summary$SEM, na.rm = TRUE)*1.25
  col_map <- c(
    "Reactive" = reactive_color,
    "Other"    = other_color
  )
  p <- ggplot() +
    geom_circle(
      data = tibble(r = seq(0.5, max_r, by = 0.5)),
      aes(x0 = 0, y0 = 0, r = r),
      color = "gray80"
    ) +
    geom_segment(
      data = angle_map,
      aes(
        x = 0, y = 0,
        xend = max_r * cos(angle),
        yend = max_r * sin(angle)
      ),
      color = "gray80",
      linewidth = 0.25
    ) +
    geom_segment(
      data = data_plot,
      aes(
        x = x_bottom, y = y_bottom,
        xend = x_top, yend = y_top,
        color = Group
      ),
      alpha = 0.5,
      linewidth = 0.5
    ) +
    geom_polygon(
      data = data_plot,
      aes(
        x = x, y = y, group = Group,
        fill = Group, color = Group
      ),
      alpha = 0.25,
      linewidth = 0.7
    ) +
    geom_text(
      data = angle_map,
      aes(
        x = 1.05 * max_r * cos(angle),
        y = 1.05 * max_r * sin(angle),
        label = Gene,
        angle = text_angle,
        hjust = hjust
      ),
      size = 2.2
    ) +
    scale_color_manual(values = col_map) +
    scale_fill_manual(values = col_map) +
    coord_fixed(clip = "off") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position = "right",
      plot.margin = margin(40, 60, 20, 60)
    )
  return(list(
    plot = p,
    filtered_genes = filtered_genes
  ))
}

find_markers <- function(object,
  groups=NULL,
  time_point=NULL,
  assay="RNA",
  logfc.threshold=0.1,
  min.pct=0.25,
  sem_multiplier=2,
  slice = 10,
  pval_cutoff=1e-3,
  logfc_cutoff=1,
  radar_scale="log2",
  xlim_range=c(-7.5,7.5)
) {
  subset_obj <- subset(object, subset = category %in% groups & time %in% time_point)
  DefaultAssay(subset_obj) <- assay
  Idents(subset_obj) <- "category"
  subset_obj <- NormalizeData(subset_obj, layer = "counts", normalization.method = "LogNormalize", scale.factor = 10000)
  subset_obj <- FindVariableFeatures(subset_obj, selection.method="vst", nfeatures=5000)
  subset_obj <- ScaleData(subset_obj)
  markers <- FindMarkers(subset_obj, ident.1=groups[1], ident.2=groups[2],
    logfc.threshold=logfc.threshold, min.pct=min.pct) %>% rownames_to_column("gene")
  markers$significant <- markers$p_val<pval_cutoff & abs(markers$avg_log2FC)>logfc_cutoff
  volcano <- ggplot(markers, aes(avg_log2FC, -log10(p_val), fill=significant)) +
    geom_point(size=3, shape=21, colour="black", alpha=0.7) +
    scale_fill_manual(values=c("gray60","red3")) +
    geom_text_repel(aes(label=ifelse(significant, gene, "")), max.overlaps=10) +
    geom_vline(xintercept=c(-logfc_cutoff,logfc_cutoff), linetype="longdash") +
    geom_hline(yintercept=-log10(pval_cutoff), linetype="longdash") +
    theme_minimal(base_size=14) +
    labs(title=paste(groups[1],"vs",groups[2]), x="log2(fc)", y="-log10(p-value)") +
    xlim(xlim_range) + theme(aspect.ratio=0.5)
  markers_to_test_gr1 <- markers %>% filter(significant & avg_log2FC>0) %>% pull(gene)
  markers_to_test_gr2 <- markers %>% filter(significant & avg_log2FC<0) %>% pull(gene)
  markers_to_test <- unique(c(markers_to_test_gr1, markers_to_test_gr2))
  object_sub <- subset(object, time==time_point)
  DefaultAssay(object_sub) <- assay
  data <- FetchData(object_sub, vars=c("category", markers_to_test)) %>%
    mutate(category_group=ifelse(category %in% groups[1], groups[1], "Other"))
  sem <- function(x) sd(x)/sqrt(length(x))
  stats <- data %>%
    group_by(category_group) %>%
    summarise(across(all_of(markers_to_test), list(mean=~mean(.), sem=~sem(.)), .names="{col}_{fn}"), .groups="drop")
  stats_gr1 <- stats %>% filter(category_group==groups[1])
  stats_gr2 <- stats %>% filter(category_group=="Other")
  if(length(markers_to_test_gr1) > 0) {
    stats_markers_gr1 <- tibble(
      gene=markers_to_test_gr1,
      mean_gr1=sapply(markers_to_test_gr1, function(g) stats_gr1[[paste0(g,"_mean")]]),
      sem_gr1=sapply(markers_to_test_gr1, function(g) stats_gr1[[paste0(g,"_sem")]]),
      mean_gr2=sapply(markers_to_test_gr1, function(g) stats_gr2[[paste0(g,"_mean")]]),
      sem_gr2=sapply(markers_to_test_gr1, function(g) stats_gr2[[paste0(g,"_sem")]])
    )
  } else {
    stats_markers_gr1 <- tibble(
      gene = character(0),
      mean_gr1 = numeric(0),
      sem_gr1  = numeric(0),
      mean_gr2 = numeric(0),
      sem_gr2  = numeric(0)
    )
  }
  if(length(markers_to_test_gr2) > 0) {
    stats_markers_gr2 <- tibble(
      gene = markers_to_test_gr2,
      mean_gr1 = as.numeric(sapply(markers_to_test_gr2, function(g) stats_gr1[[paste0(g, "_mean")]])),
      sem_gr1  = as.numeric(sapply(markers_to_test_gr2, function(g) stats_gr1[[paste0(g, "_sem")]])),
      mean_gr2 = as.numeric(sapply(markers_to_test_gr2, function(g) stats_gr2[[paste0(g, "_mean")]])),
      sem_gr2  = as.numeric(sapply(markers_to_test_gr2, function(g) stats_gr2[[paste0(g, "_sem")]]))
    )
  } else {
    stats_markers_gr2 <- tibble(
      gene = character(0),
      mean_gr1 = numeric(0),
      sem_gr1  = numeric(0),
      mean_gr2 = numeric(0),
      sem_gr2  = numeric(0)
    )
  }
  genes_gr1_filtered <- stats_markers_gr1 %>%
    mutate(
      pass_filter =
      (mean_gr1 > mean_gr2) &
      ((mean_gr1 - sem_multiplier * sem_gr1) > (mean_gr2 + sem_multiplier * sem_gr2))
    ) %>%
    filter(pass_filter)
  genes_gr2_filtered <- stats_markers_gr2 %>%
    mutate(
      pass_filter =
      (mean_gr2 > mean_gr1) &
      ((mean_gr2 - sem_multiplier * sem_gr2) > (mean_gr1 + sem_multiplier * sem_gr1))
    ) %>%
    filter(pass_filter) %>%
    arrange(desc(mean_gr2)) %>%
    slice_head(n = slice)
  plot_group_genes <- function(stats_table, group_name) {
    genes <- stats_table$gene
    if(length(genes)==0) return(list(plot=NULL, genes=character(0)))
    if(length(genes)>=3) {
      radar <- tryCatch(
        radar_compare_reactive(object_sub, gene_list=genes, filter_up=FALSE, scale=radar_scale),
        error=function(e) NULL
      )
      if(is.null(radar)) return(list(plot=NULL, genes=character(0)))
      return(list(plot=if(is.list(radar) && "plot" %in% names(radar)) radar$plot else radar, genes=genes))
    } else {
      plot_data <- stats_table %>% filter(pass_filter) %>% select(gene, mean_gr1, sem_gr1, mean_gr2, sem_gr2) %>%
        pivot_longer(cols=-gene, names_to=c("stat","category"), names_sep="_", values_to="value") %>%
        pivot_wider(names_from=stat, values_from=value) %>%
        mutate(category=recode(category,"gr1"="Reactive","gr2"="Other"))
      fill_colors <- c("Reactive"="red3","Other"="gray60")
      plot <- ggplot(plot_data, aes(x=gene, y=mean, fill=category)) +
        geom_bar(stat="identity", position=position_dodge(0.8), width=0.8, color="black", linewidth=0.5, alpha=0.75) +
        geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position=position_dodge(0.8), width=0.2) +
        scale_fill_manual(values=fill_colors) + theme_minimal(base_size=14) + theme(aspect.ratio=1) +
        labs(x="Gene(s)", y="Expression")
      return(list(plot=plot, genes=genes))
    }
  }
  gr1_result <- plot_group_genes(stats_table = genes_gr1_filtered, group_name = groups[1])
  gr2_result <- plot_group_genes(stats_table = genes_gr2_filtered, group_name = "Other")
  gr1_plot <- if(!is.null(gr1_result$plot)) gr1_result$plot + labs(title="Mean ± SEM for group 1") else NULL
  gr2_plot <- if(!is.null(gr2_result$plot)) gr2_result$plot + labs(title="Mean ± SEM for group 2") else NULL
  list(
    plot = wrap_plots(list(volcano, gr1_plot, gr2_plot)[!sapply(list(volcano, gr1_plot, gr2_plot), is.null)], ncol=3),
    markers_group1 = gr1_result$genes,
    markers_group2 = gr2_result$genes
  )
}

refine_gene_set <- function(genes) {
  genes %>%
    .[!grepl("^TRAV|^TRBV", .)] %>%
    .[order(grepl("-$", .), .)]
}

run_dimred_cluster <- function(x, ndim, res, seed = SEED, n_neighbors = 60, min_dist = 0.5, k_param = 30, reduction_name = "umap") {
  set.seed(seed)
  x <- RunPCA(
    x,
    npcs = ndim,
    seed.use = seed,
    verbose = FALSE
  )
  x <- FindNeighbors(
    x,
    dims = 1:ndim,
    k.param = k_param,
    verbose = FALSE
  )
  set.seed(seed)
  x <- FindClusters(
    x,
    resolution = res,
    random.seed = seed,
    verbose = FALSE
  )
  set.seed(seed)
  x <- RunUMAP(
    x,
    dims = 1:ndim,
    n.neighbors = n_neighbors,
    min.dist = min_dist,
    metric = "cosine",
    nn.method = "rann",
    reduction.name = reduction_name,
    reduction.key = paste0(reduction_name, "_"),
    seed.use = seed,
    verbose = FALSE
  )
  return(x)
}

quietTCRgenes <- function(sc) {
  unwanted_genes <- "^TR[ABDG][VDJ][^D]"
  unwanted_genes <- grep(pattern = unwanted_genes, x = sc, value = TRUE)
  sc <- sc[!sc %in% unwanted_genes]
  return(sc)
}
