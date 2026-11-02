# --- Clean FlowJo gate paths into short column names (e.g. CD8_PD1hi_CD39) ---
clean_flowjo_name_PD1 <- function(x, main_subset = c("CD4", "CD8")) {
  main_subset <- match.arg(main_subset)
  x %>%
    str_remove("^Lymphocytes/Single Cells/Alive/") %>%
    str_remove("\\s*\\|.*$") %>%
    str_replace_all("\\+", "") %>%
    str_squish() %>%
    str_remove("^CD3/") %>%
    str_replace_all("/", " ") %>%
    str_squish() %>%
    str_replace("^CD3 CD8", "CD8") %>%
    str_replace(
      paste0("^", main_subset, " PD1hi ", main_subset, " PD1hi"),
      paste0(main_subset, " PD1hi")
    ) %>%
    str_replace(
      paste0("^", main_subset, " PD1int ", main_subset, " PD1int"),
      paste0(main_subset, " PD1int")
    ) %>%
    str_replace(
      paste0("^", main_subset, " PD1neg ", main_subset, " PD1neg"),
      paste0(main_subset, " PD1neg")
    ) %>%
    str_replace(
      paste0("^", main_subset, " ", main_subset, " "),
      paste0(main_subset, " ")
    ) %>%
    str_replace(paste0("^", main_subset, " PD1hi$"), paste0(main_subset, "_PD1hi")) %>%
    str_replace(paste0("^", main_subset, " PD1int$"), paste0(main_subset, "_PD1int")) %>%
    str_replace(paste0("^", main_subset, " PD1neg$"), paste0(main_subset, "_PD1neg")) %>%
    str_replace(paste0("^", main_subset, " PD1hi (.*)$"), paste0(main_subset, "_PD1hi_\\1")) %>%
    str_replace(paste0("^", main_subset, " PD1int (.*)$"), paste0(main_subset, "_PD1int_\\1")) %>%
    str_replace(paste0("^", main_subset, " PD1neg (.*)$"), paste0(main_subset, "_PD1neg_\\1")) %>%
    str_replace_all(" ", "_") %>%
    str_replace_all(paste0(main_subset, "_PD1hi_", main_subset, "_PD1hi_"), paste0(main_subset, "_PD1hi_")) %>%
    str_replace_all(paste0(main_subset, "_PD1int_", main_subset, "_PD1int_"), paste0(main_subset, "_PD1int_")) %>%
    str_replace_all(paste0(main_subset, "_PD1neg_", main_subset, "_PD1neg_"), paste0(main_subset, "_PD1neg_")) %>%
    str_replace("^CD3_CD8$", "CD8")
}

# --- Read a FlowJo export panel, drop summary rows and add the day column ---
read_panel <- function(file, day_vec, main_subset, dir) {
  path <- file.path(dir, file)
  first_line <- readLines(path, n = 1)
  delim <- ifelse(str_count(first_line, ";") > str_count(first_line, ","), ";", ",")
  df <- read_delim(path, delim = delim, show_col_types = FALSE)
  names(df)[1] <- "sample"
  df <- df %>%
    { if (ncol(.) >= 24) select(., -24) else . } %>%
    filter(!sample %in% c("Mean", "SD")) %>%
    select(-sample)
  if (nrow(df) != length(day_vec)) {
    stop(
      file, " has ", nrow(df), " rows after filtering, but day_vec has ",
      length(day_vec), " values."
    )
  }
  df <- df %>%
    mutate(day = day_vec, .before = 1) %>%
    arrange(day)
  names(df) <- clean_flowjo_name_PD1(names(df), main_subset = main_subset)
  df
}

# --- Suffix the subset and PD-1 base columns with the panel number ---
rename_panel_base <- function(df, main_subset, panel) {
  old_names <- c(
    main_subset,
    paste0(main_subset, "_PD1hi"),
    paste0(main_subset, "_PD1int"),
    paste0(main_subset, "_PD1neg")
  )
  new_names <- c(
    paste0(main_subset, "_pan", panel),
    paste0(main_subset, "_PD1hi_pan", panel),
    paste0(main_subset, "_PD1int_pan", panel),
    paste0(main_subset, "_PD1neg_pan", panel)
  )
  names(df) <- recode(names(df), !!!setNames(new_names, old_names))
  df
}

# --- Add per-marker proportions of each PD-1 group (_prop columns) ---
compute_marker_props <- function(df, main_subset, pd1_groups, markers) {
  for (marker in markers) {
    cols <- paste0(main_subset, "_", pd1_groups, "_", marker)
    cols <- cols[cols %in% names(df)]
    if (length(cols) == 0) next
    denom <- rowSums(df[cols], na.rm = TRUE)
    for (col in cols) {
      df[[paste0(col, "_prop")]] <- df[[col]] / denom
    }
  }
  df
}

# --- Reshape the _prop columns to long format (day, PD1_cat, marker, prop) ---
pivot_marker_props <- function(df, main_subset) {
  df %>%
    select(day, ends_with("_prop")) %>%
    pivot_longer(
      -day,
      names_to = c("subset", "PD1_cat", "marker"),
      names_pattern = paste0("^(", main_subset, ")_(PD1hi|PD1int|PD1neg)_(.*)_prop$"),
      values_to = "prop"
    ) %>%
    select(day, PD1_cat, marker, prop)
}

# --- PD-1 group lines + daily pies (Rmd globals: pheno_days, vline_df, pd1_labels_flow, theme_manuscript, plot_font) ---
make_marker_pie_plot <- function(df_long, selected_marker, subset_label = "CD8") {
  plot_data <- df_long %>%
    filter(marker == selected_marker) %>%
    mutate(
      PD1_cat = factor(PD1_cat, levels = c("PD1hi", "PD1int", "PD1neg")),
      day_f = factor(day, levels = pheno_days)
    )
  p_line <- ggplot(plot_data) +
    geom_vline(
      data = vline_df,
      aes(xintercept = day),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.25,
      colour = "grey85"
    ) +
    geom_line(
      aes(x = day, y = 100 * prop, colour = PD1_cat, group = PD1_cat),
      linewidth = 0.45
    ) +
    geom_point(
      aes(x = day, y = 100 * prop, colour = PD1_cat, group = PD1_cat),
      size = 1.4
    ) +
    scale_x_continuous(
      breaks = pheno_days,
      labels = pheno_days,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_colour_discrete(
      breaks = names(pd1_labels_flow),
      labels = unname(pd1_labels_flow)
    ) +
    labs(
      title = selected_marker,
      x = NULL,
      y = paste0("Frequency of ", subset_label, "+ lymphocytes (%)"),
      colour = NULL
    ) +
    theme_manuscript +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    )
  p_pies <- ggplot(plot_data, aes(x = "", y = prop, fill = PD1_cat)) +
    geom_col(width = 1, color = NA) +
    coord_polar(theta = "y") +
    facet_wrap(~ day_f, nrow = 1) +
    scale_fill_discrete(
      breaks = names(pd1_labels_flow),
      labels = unname(pd1_labels_flow)
    ) +
    theme_void(base_family = plot_font) +
    theme(
      legend.position = "none",
      strip.text = element_blank(),
      panel.spacing = unit(0, "pt"),
      panel.spacing.x = unit(0, "pt"),
      panel.spacing.y = unit(0, "pt"),
      plot.margin = margin(-10, 5.5, 0, 5.5)
    )
  p_line / p_pies + plot_layout(heights = c(1, 0.18))
}

# --- Z-score change vs day 0 per marker (Rmd globals: pheno_days, vline_df, pd1_labels_counts, theme_manuscript) ---
make_zscore_variation_plot <- function(df, subset_label = "CD4") {
  zscore_marker_order <- c(
    "CCR7", "CD27", "CD28", "CD38",
    "CD39", "CD45RO", "DR", "HLADR", "FoxP3",
    "ICOS", "KI67", "Ki67", "LAG3", "TIM3"
  )
  marker_levels <- if (is.factor(df$marker)) {
    levels(df$marker)[levels(df$marker) %in% as.character(unique(df$marker))]
  } else {
    c(
      zscore_marker_order[zscore_marker_order %in% unique(df$marker)],
      setdiff(unique(df$marker), zscore_marker_order)
    )
  }
  baseline_df <- df %>%
    filter(day == 0) %>%
    group_by(marker, pd1_level) %>%
    summarise(
      baseline_lymphocyte_pc = first(lymphocyte_pc),
      .groups = "drop"
    )
  sd_df <- df %>%
    group_by(marker, pd1_level) %>%
    summarise(
      lymphocyte_pc_sd = sd(lymphocyte_pc, na.rm = TRUE),
      .groups = "drop"
    )
  plot_data <- df %>%
    filter(day %in% pheno_days) %>%
    left_join(baseline_df, by = c("marker", "pd1_level")) %>%
    left_join(sd_df, by = c("marker", "pd1_level")) %>%
    mutate(
      pd1_level = factor(pd1_level, levels = c("PD1_hi", "PD1_int", "PD1_neg")),
      marker = factor(marker, levels = marker_levels),
      zscore_variation = case_when(
        is.na(baseline_lymphocyte_pc) ~ NA_real_,
        is.na(lymphocyte_pc_sd) | lymphocyte_pc_sd == 0 ~ 0,
        TRUE ~ (lymphocyte_pc - baseline_lymphocyte_pc) / lymphocyte_pc_sd
      )
    )
  ggplot(plot_data) +
    geom_vline(
      data = vline_df,
      aes(xintercept = day),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.25,
      colour = "grey85"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.25,
      colour = "grey35"
    ) +
    geom_line(
      aes(x = day, y = zscore_variation, colour = pd1_level, group = pd1_level),
      linewidth = 0.45
    ) +
    geom_point(
      aes(x = day, y = zscore_variation, colour = pd1_level, group = pd1_level),
      size = 1.4
    ) +
    facet_wrap(~ marker, ncol = 5) +
    scale_x_continuous(
      breaks = pheno_days,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_continuous(
      breaks = c(-2, 0, 2),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    coord_cartesian(ylim = c(-3, 3)) +
    scale_colour_discrete(
      breaks = names(pd1_labels_counts),
      labels = unname(pd1_labels_counts)
    ) +
    labs(
      x = "day",
      y = paste0(subset_label, "+ lymphocytes z-score variation vs day 0"),
      colour = NULL
    ) +
    theme_manuscript
}
