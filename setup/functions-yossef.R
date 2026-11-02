# --- Row of tcr_reactivity whose TRA+TRB CDR3s the cell carries for its patient, or NA ---
match_reactive_tcr_yossef <- function(patient, chains, cdr3, tcr_reactivity) {
  chain_split <- strsplit(as.character(chains), ";", fixed = TRUE)
  cdr3_split <- strsplit(as.character(cdr3), ";", fixed = TRUE)
  vapply(seq_along(patient), function(i) {
    if (is.na(cdr3[i]) || is.na(chains[i])) return(NA_integer_)
    tra <- cdr3_split[[i]][chain_split[[i]] == "TRA"]
    trb <- cdr3_split[[i]][chain_split[[i]] == "TRB"]
    hit <- which(tcr_reactivity$patient == patient[i] & tcr_reactivity$TCRa %in% tra & tcr_reactivity$TCRb %in% trb)
    if (length(hit) == 0) NA_integer_ else hit[1]
  }, integer(1))
}

# --- Yossef version of plot_reactive_top(): clonotype-level, filtered by patient ---
plot_reactive_top_yossef <- function(
  object,
  pivot_cols = NULL,
  top_n = NULL,
  sample = NULL
) {
  df_long <- object@meta.data %>%
    select(reactive, all_of(pivot_cols), category, clonotype, TCRb, patient) %>%
    filter(patient %in% sample) %>%
    pivot_longer(
      cols = all_of(pivot_cols),
      names_to = "feature",
      values_to = "value"
    )
  max_per_clonotype <- df_long %>%
    group_by(feature, clonotype) %>%
    summarise(
      max_value = max(value, na.rm = TRUE),
      reactive_clonotype = any(category == "Reactive", na.rm = TRUE),
      .groups = "drop"
    )
  top_info <- max_per_clonotype %>%
    group_by(feature) %>%
    arrange(desc(max_value), clonotype) %>%
    slice_head(n = top_n) %>%
    mutate(is_top = TRUE) %>%
    ungroup()
  thresholds <- top_info %>%
    group_by(feature) %>%
    summarise(thresh_value = min(max_value), .groups = "drop") %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  df_long <- df_long %>%
    left_join(top_info %>% select(feature, clonotype, is_top),
      by = c("feature", "clonotype")) %>%
    mutate(is_top = ifelse(is.na(is_top), FALSE, TRUE)) %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  top_clonotypes_with_key <- df_long %>%
    filter(is_top) %>%
    distinct(feature, clonotype) %>%
    group_by(feature) %>%
    summarise(top_clonotypes = list(clonotype), .groups = "drop") %>%
    { setNames(.$top_clonotypes, .$feature) }
  top_clonotype_table <- df_long %>%
    filter(is_top) %>%
    group_by(feature, clonotype, category) %>%
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
    summarise(max_reactive = n_distinct(clonotype), .groups = "drop") %>%
    mutate(feature = factor(feature, levels = pivot_cols))
  id_top_status <- df_long %>%
    group_by(feature, clonotype) %>%
    summarise(
      id_is_reactive = any(category == "Reactive"),
      id_is_top = any(is_top),
      .groups = "drop"
    )
  reactive_in_top <- id_top_status %>%
    filter(id_is_top & id_is_reactive) %>%
    group_by(feature) %>%
    summarise(reactive_top = n_distinct(clonotype), .groups = "drop") %>%
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
      "_top_", top_n))
  reactive_id_table <- df_long %>%
    left_join(thresholds, by = "feature") %>%
    filter(category == "Reactive") %>%
    group_by(feature, clonotype) %>%
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

# --- Yossef version of plot_reactive_tcrs_vs_tested() for one patient ---
plot_reactive_tcrs_vs_tested_yossef <- function(object,
  step = 250,
  max_n = 50000,
  pivot_cols = NULL,
  sample_id = NULL
) {
  meta <- object@meta.data
  meta <- meta %>% subset(meta$patient == sample_id)
  required_cols <- c("clonotype", "reactive", "category", pivot_cols)
  if (!all(required_cols %in% colnames(meta))) {
    stop("meta.data must contain: ", paste(required_cols, collapse = ", "))
  }
  meta <- meta %>% select(clonotype, category, pivot_cols)
  n_clonotypes <- n_distinct(meta$clonotype)
  max_n <- min(max_n, n_clonotypes)
  x_vals <- seq(0, max_n, by = step)
  all_clonotypes <- meta %>%
    distinct(clonotype, category) %>%
    slice_sample(n = n_clonotypes)
  # --- Random baseline ---
  baseline <- lapply(x_vals, function(n) {
    if (n == 0) {
      return(tibble(
        TCR_tested = 0,
        Reactive_identified = 0,
        Method = "Random"
      ))
    }
    top <- all_clonotypes %>% slice_head(n = n)
    reactive_count <- meta %>%
      filter(clonotype %in% top$clonotype, category == "Reactive") %>%
      summarise(n = n_distinct(clonotype)) %>%
      pull(n)
    tibble(
      TCR_tested = n,
      Reactive_identified = reactive_count,
      Method = "Random"
    )
  }) %>% bind_rows()
  # --- Signature-ranked curves ---
  scored_results <- lapply(pivot_cols, function(score_col) {
    df <- meta
    if (!"cell_id" %in% colnames(df)) {
      df <- tibble::rownames_to_column(df, "cell_id")
    }
    df <- df %>%
      select(cell_id, clonotype, category, !!sym(score_col)) %>%
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
      this_key <- df$clonotype[i]
      if (!(this_key %in% seen_keys)) {
        seen_keys <- c(seen_keys, this_key)
      }
      n_unique_keys <- length(seen_keys)
      if (n_unique_keys %in% x_vals && !(n_unique_keys %in% counted_keys)) {
        counted_keys <- c(counted_keys, n_unique_keys)
        reactive_count <- df %>%
          filter(clonotype %in% seen_keys, category == "Reactive") %>%
          summarise(n = n_distinct(clonotype)) %>%
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
    geom_line(linewidth = 0.5) +
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
    scale_color_manual(values = signature_colors) +
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

# --- One row per clonotype: reactive if any cell is Reactive, score = mean across its cells ---
clonotype_scores_yossef <- function(meta, pivot_cols) {
  meta %>%
    group_by(clonotype) %>%
    summarise(
      is_reactive = as.integer(any(category == "Reactive")),
      across(all_of(pivot_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

# --- Per-patient clonotype-level AUC (as in plot_ROC_curves_yossef) with DeLong 95% CI ---
auc_table_yossef <- function(object, pivot_cols) {
  lapply(unique(object$patient), function(i) {
    clonotype_df <- clonotype_scores_yossef(object@meta.data %>% filter(patient == i), pivot_cols)
    n_reactive <- sum(clonotype_df$is_reactive)
    lapply(pivot_cols, function(score_col) {
      ci <- c(NA_real_, NA_real_, NA_real_)
      if (n_reactive > 0 && n_reactive < nrow(clonotype_df)) {
        r <- roc(clonotype_df$is_reactive, clonotype_df[[score_col]], levels = c(0, 1), direction = "<", quiet = TRUE)
        ci <- as.numeric(ci.auc(r, method = "delong"))
      }
      data.frame(patient = i, signature = score_col, n_clonotypes = nrow(clonotype_df), n_reactive = n_reactive,
        AUC = ci[2], CI_low = max(ci[1], 0), CI_high = min(ci[3], 1))
    }) %>% bind_rows()
  }) %>% bind_rows()
}

# --- Across-patient AUC summary and paired Wilcoxon test of each signature vs each reference ---
summarise_auc_yossef <- function(auc_df, references) {
  auc_df <- auc_df %>% filter(!is.na(AUC))
  summary_df <- auc_df %>%
    group_by(signature) %>%
    summarise(n_patients = n(), median_AUC = median(AUC), min_AUC = min(AUC), max_AUC = max(AUC), .groups = "drop")
  wide <- tidyr::pivot_wider(auc_df, id_cols = patient, names_from = signature, values_from = AUC)
  tests <- expand.grid(signature = setdiff(unique(auc_df$signature), references), reference = references, stringsAsFactors = FALSE)
  tests$median_diff <- NA_real_
  tests$n_higher <- NA_integer_
  tests$p_value <- NA_real_
  for (k in seq_len(nrow(tests))) {
    d <- wide[[tests$signature[k]]] - wide[[tests$reference[k]]]
    d <- d[!is.na(d)]
    tests$median_diff[k] <- median(d)
    tests$n_higher[k] <- sum(d > 0)
    tests$p_value[k] <- suppressWarnings(wilcox.test(d, mu = 0, exact = TRUE)$p.value)
  }
  tests$p_adj_BH <- p.adjust(tests$p_value, method = "BH")
  list(summary = summary_df, tests = tests)
}

# --- Yossef version of plot_ROC_curves() for one patient (clonotype level) ---
plot_ROC_curves_yossef <- function(object,
  pivot_cols = NULL,
  sample_id = NULL
) {
  meta <- object@meta.data
  meta <- meta %>% subset(meta$patient == sample_id)
  required_cols <- c("clonotype", "reactive", "category", pivot_cols)
  if (!all(required_cols %in% colnames(meta))) {
    stop("Missing required columns in meta.data")
  }
  clonotype_df <- clonotype_scores_yossef(meta, pivot_cols)
  # --- No ROC without both classes (e.g. a patient with no reactive clonotype after subsetting) ---
  if (n_distinct(clonotype_df$is_reactive) < 2) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No reactive clonotypes") +
        labs(title = sample_id, x = "1 - Specificity", y = "Sensitivity") +
        theme_minimal(base_size = 13) +
        theme(aspect.ratio = 1)
    )
  }
  roc_objs <- lapply(pivot_cols, function(score_col) {
    roc(clonotype_df$is_reactive, clonotype_df[[score_col]], levels = c(0, 1), direction = "<", quiet = TRUE)
  })
  names(roc_objs) <- pivot_cols
  roc_results <- lapply(pivot_cols, function(score_col) {
    coords_df <- coords(roc_objs[[score_col]], "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
    coords_df$score <- score_col
    coords_df$auc <- auc(roc_objs[[score_col]])
    coords_df
  }) %>%
    bind_rows() %>%
    mutate(score = factor(score, levels = pivot_cols))
  # --- AUC per signature ---
  auc_table <- roc_results %>%
    group_by(score) %>%
    summarise(AUC = unique(round(auc, 3)))
  # --- AUC labels coloured like their curve, so each value can be matched to its signature ---
  auc_table <- auc_table %>%
    mutate(
      score = factor(score, levels = pivot_cols),
      label = paste0("AUC = ", format(AUC, nsmall = 3))
    ) %>%
    arrange(score) %>%
    mutate(y = 0.05 + 0.07 * (n() - row_number()))
  p <- ggplot(roc_results, aes(x = 1 - specificity, y = sensitivity, color = score)) +
    geom_line(linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_text(data = auc_table, aes(x = 0.6, y = y, label = label, color = score), hjust = 0, size = 2.8, show.legend = FALSE, inherit.aes = FALSE) +
    labs(
      title = sample_id,
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Signature"
    ) +
    theme_minimal(base_size = 13) +
    theme(aspect.ratio = 1, legend.position = "bottom") +
    scale_color_manual(values = signature_colors)
  return(p)
}
