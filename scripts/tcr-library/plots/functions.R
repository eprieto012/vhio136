# --- Libraries ---
library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(patchwork)
library(purrr)
library(readxl)
library(scales)
library(stringr)
library(tibble)
library(tidyr)

# --- Helpers and functions used by script(14).Rmd ---

tcr_clonotype_label <- function(clonotype_by) {
  ifelse(clonotype_by == "aminoAcid", "AA", "nt")
}

tcr_format_timepoint <- function(x) {
  ifelse(grepl("^d", as.character(x)), as.character(x), paste0("d", x))
}

tcr_timepoint_days <- function(tp_levels) {
  as.numeric(gsub("[^0-9]", "", tp_levels))
}

tcr_timepoint_levels <- function(patient_list) {
  tp_levels <- names(patient_list)
  tp_days <- tcr_timepoint_days(tp_levels)
  tp_levels[order(tp_days)]
}

tcr_y_scale_log <- function(floor_y = 1e-3) {
  scale_y_log10(limits = c(floor_y, 10),
                breaks = c(1e-3, 1e-2, 1e-1, 1, 10),
                labels = c("0.001", "0.01", "0.1", "1", "10"))
}

tcr_theme_trajectory <- function(base_size = 8, aspect.ratio = 0.5) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(size = 7, colour = "black"),
          axis.title.y = element_text(size = 8, colour = "black"),
          axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 9.5, face = "plain", colour = "black", hjust = 0.5),
          legend.position = "none",
          aspect.ratio = aspect.ratio)
}

tcr_theme_rank <- function(base_size = 9, aspect.ratio = 1, legend.position = "right") {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 11, face = "plain", colour = "black", hjust = 0.5),
          legend.position = legend.position,
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.size = grid::unit(0.45, "cm"),
          aspect.ratio = aspect.ratio)
}

load_data <- function(samples,
                      datapath.pbl,
                      output.list.name,
                      technology = "adaptive",
                      min_freq = 0.0001) {
  min_freq <<- min_freq
  if (!exists(output.list.name, envir = .GlobalEnv)) {
    assign(output.list.name, list(), envir = .GlobalEnv)
  }
  output_list <- get(output.list.name, envir = .GlobalEnv)
  fc_output_list_name <- paste0(output.list.name, "_fc")
  if (!exists(fc_output_list_name, envir = .GlobalEnv)) {
    assign(fc_output_list_name, list(), envir = .GlobalEnv)
  }
  fc_output_list <- get(fc_output_list_name, envir = .GlobalEnv)
  message("⏳ loading data...")
  for (i in samples) {
    if (technology == "adaptive") {
      datapath_i <- file.path(datapath.pbl, i)
      files <- list.files(path = datapath_i, pattern = "\\.tsv$", full.names = TRUE)
      files <- files[order(as.numeric(sub(".*d(\\d+)\\.tsv$", "\\1", files)))]
      tcr_data_i <- lapply(files, function(x) {
        tcr_data <- read.table(file = x, sep = "\t", header = TRUE)
        tcr_data <- subset(tcr_data, subset = tcr_data$sequenceStatus == "In")
        if ("count" %in% colnames(tcr_data)) {
          tcr_data$freq_in <- tcr_data$count / sum(tcr_data$count) * 100
        } else if ("count..reads." %in% colnames(tcr_data)) {
          tcr_data$freq_in <- tcr_data$count..reads. / sum(tcr_data$count..reads.) * 100
        } else if ("count..templates.reads." %in% colnames(tcr_data)) {
          tcr_data$freq_in <- tcr_data$count..templates.reads. / sum(tcr_data$count..templates.reads.) * 100
        } else {
          stop("invalid format")
        }
        tcr_data <- tcr_data %>% filter(freq_in >= min_freq)
        tcr_data$time <- tools::file_path_sans_ext(basename(x))
        if ("count" %in% colnames(tcr_data)) {
          tcr_data <- subset(tcr_data, select = c("nucleotide", "aminoAcid", "time", "count", "freq_in", "vGeneName", "jGeneName"))
        } else if ("count..reads." %in% colnames(tcr_data)) {
          tcr_data <- subset(tcr_data, select = c("nucleotide", "aminoAcid", "time", "count..reads.", "freq_in", "vGeneName", "jGeneName"))
        } else if ("count..templates.reads." %in% colnames(tcr_data)) {
          tcr_data <- subset(tcr_data, select = c("nucleotide", "aminoAcid", "time", "count..templates.reads.", "freq_in", "vGeneName", "jGeneName"))
        } else {
          stop("invalid format")
        }
        tcr_data$nucleotide <- paste0(tcr_data$nucleotide, "_", tcr_data$vGeneName, "_", tcr_data$jGeneName)
        tcr_data$vGeneName <- NULL
        tcr_data$jGeneName <- NULL
        colnames(tcr_data) <- c("nucleotide", "aminoAcid", "time", "count", "freq_in")
        return(tcr_data)
      })
    } else if (technology == "hsd") {
      datapath_i <- file.path(datapath.pbl, i)
      files <- list.files(path = datapath_i, pattern = "\\.xlsx$", full.names = TRUE)
      files <- files[order(as.numeric(sub(".*d(\\d+)\\.xlsx$", "\\1", files)))]
      tcr_data_i <- lapply(files, function(x) {
        tcr_data <- read_excel(x, col_names = TRUE)
        tcr_data <- as.data.frame(tcr_data)
        tcr_data <- subset(tcr_data, subset = tcr_data$`Frame-info` == "ok")
        tcr_data$freq_in <- tcr_data$"total.nr-reads" / sum(tcr_data$"total.nr-reads") * 100
        # Filter by min_freq_in
        tcr_data <- tcr_data %>% filter(freq_in >= min_freq)
        tcr_data$time <- tools::file_path_sans_ext(basename(x))
        tcr_data$`DNA-seq (AA2)` <- paste0(tcr_data$`DNA-seq (AA2)`, "_", tcr_data$`Vsegment-ID`, "_", tcr_data$`Jsegment-ID`)
        tcr_data <- tcr_data %>% select("DNA-seq (AA2)", "AA1: pept-seq Insert", "time", "total.nr-reads", "freq_in")
        colnames(tcr_data) <- c("nucleotide", "aminoAcid", "time", "count", "freq_in")
        return(tcr_data)
      })
    } else if (technology == "SEQTR") {
      datapath_i <- file.path(datapath.pbl, i)
      files <- list.files(path = datapath_i, pattern = "\\.txt$", full.names = TRUE)
      files <- files[order(as.numeric(sub(".*d(\\d+)\\.txt$", "\\1", files)))]
      tcr_data_i <- lapply(files, function(x) {
        tcr_data <- read.table(x, sep = "\t", fill = TRUE, na.strings = c("", "NA"))
        colnames(tcr_data) <- c("CDR3_sequence", "Count", "TRBV",	"TRBJ",	"Frame", "CDR3_aaseq", "CDR3_length",	"UMI", "CDR3_nt_nbr")
        tcr_data <- as.data.frame(tcr_data)
        tcr_data <- subset(tcr_data, subset = tcr_data$Frame == "IN")
        tcr_data <- tcr_data %>% filter(Count > 1)
        tcr_data$freq_in <- tcr_data$Count / sum(tcr_data$Count) * 100
        # Filter by min_freq_in
        tcr_data <- tcr_data %>% filter(freq_in >= min_freq)
        tcr_data$time <- tools::file_path_sans_ext(basename(x))
        tcr_data$CDR3_sequence <- paste0(tcr_data$CDR3_sequence, "_", tcr_data$TRBV, "_", tcr_data$TRBJ)
        tcr_data$aminoAcid_sequence <- paste0(tcr_data$CDR3_aaseq, "_", tcr_data$TRBV, "_", tcr_data$TRBJ)
        tcr_data <- tcr_data %>% select("CDR3_sequence", "aminoAcid_sequence", "time", "Count", "freq_in")
        colnames(tcr_data) <- c("nucleotide", "aminoAcid", "time", "count", "freq_in")
        tcr_data <- tcr_data[!is.na(rownames(tcr_data)), ]
        return(tcr_data)
      })
    } else if (technology == "MiXCR") {
      datapath_i <- file.path(datapath.pbl, i)
      files <- list.files(path = datapath_i, pattern = "\\.tsv$", full.names = TRUE)
      files <- files[order(as.numeric(sub(".*d(\\d+)\\.tsv$", "\\1", files)))]
      tcr_data_i <- lapply(files, function(x) {
        tcr_data <- read.table(file = x, sep = "\t", header = TRUE)
        colnames(tcr_data) <- c("aaSeqCDR3", "vGene", "jGene", "nSamples", "readCountAggregated", "readFractionAggregated", "uniqueUMICountAggregated", "uniqueUMIFractionAggregated", "tagCounts", "bestVGene", "bestJGene")
        tcr_data <- as.data.frame(tcr_data)
        tcr_data <- tcr_data %>% filter(readCountAggregated > 10)
        tcr_data$freq_in <- tcr_data$readFractionAggregated / sum(tcr_data$readFractionAggregated) * 100
        # Filter by min_freq_in
        tcr_data <- tcr_data %>% filter(freq_in >= min_freq)
        tcr_data$time <- tools::file_path_sans_ext(basename(x))
        tcr_data$CDR3_sequence <- paste0(tcr_data$aaSeqCDR3, "_", tcr_data$bestVGene, "_", tcr_data$bestJGene)
        tcr_data$aminoAcid_sequence <- paste0(tcr_data$aaSeqCDR3, "_", tcr_data$bestVGene, "_", tcr_data$bestJGene)
        tcr_data <- tcr_data %>% select("CDR3_sequence", "aminoAcid_sequence", "time", "readCountAggregated", "freq_in")
        colnames(tcr_data) <- c("nucleotide", "aminoAcid", "time", "count", "freq_in")
        return(tcr_data)
      })
    } else {
      message("🚧 invalid technology; only adaptive, hsd, SEQTR or MiXCR data are supported...")
      return(NULL)
    }
    # check if data has been loaded
    if (length(tcr_data_i) == 0) {
      message("🚧 invalid technology or input format; only only adaptive, hsd, SEQTR or MiXCR data are supported...")
      return(NULL)
    }
    names(tcr_data_i) <- tools::file_path_sans_ext(basename(files))
    output_list[[i]] <- tcr_data_i
    data <- Reduce(rbind, tcr_data_i)
    data$count <- NULL
    data <- aggregate(freq_in ~ nucleotide + time, data, sum)
    data <- as.data.frame(tidyr::pivot_wider(data, names_from = time, values_from = freq_in))
    data[is.na(data)] <- min_freq
    freq_data <- data
    timepoint_values <- sub("d", "", colnames(data)[2:ncol(data)])
    timepoint_numeric <- as.numeric(timepoint_values)
    min_timepoint_index <- which.min(timepoint_numeric)
    fc <- data[, 2:ncol(data)] / data[, min_timepoint_index + 1]
    fc <- log2(fc)
    fc$nucleotide <- data$nucleotide
    timepoints <- colnames(fc)[1:(ncol(fc) - 1)]
    fc_long <- as.data.frame(tidyr::pivot_longer(fc, cols = all_of(timepoints), names_to = "timepoints", values_to = "fc"))
    freq_long <- tidyr::pivot_longer(freq_data, cols = all_of(timepoints), names_to = "timepoints", values_to = "freq_in")
    fc_long <- dplyr::left_join(fc_long, freq_long, by = c("nucleotide", "timepoints"))
    fc_output_list[[i]] <- fc_long
    message(paste0(grep(i, samples), "/", length(samples), " ", i))
  }
  assign(output.list.name, output_list, envir = .GlobalEnv)
  assign(fc_output_list_name, fc_output_list, envir = .GlobalEnv)
  message("✅ done!")
}

load_data_biopsy <- function(samples,
                             datapath.biopsy,
                             output.list.name,
                             min_freq = 0.0001,
                             file_pattern = "TCRB in frame\\.tsv$") {
  min_freq <<- min_freq
  if (!exists(output.list.name, envir = .GlobalEnv)) {
    assign(output.list.name, list(), envir = .GlobalEnv)
  }
  output_list <- get(output.list.name, envir = .GlobalEnv)
  message("⏳ loading biopsy data...")
  for (i in samples) {
    datapath_i <- file.path(datapath.biopsy, i)
    files <- list.files(path = datapath_i, pattern = file_pattern, full.names = TRUE)
    files <- files[order(basename(files))]
    if (length(files) == 0) {
      message(paste0("🚧 no valid biopsy files found for ", i))
      next
    }
    tcr_data_i <- lapply(files, function(x) {
      tcr_data <- read.table(file = x, sep = "\t", header = TRUE, check.names = TRUE)
      if ("sequenceStatus" %in% colnames(tcr_data)) {
        tcr_data <- subset(tcr_data, subset = tcr_data$sequenceStatus == "In")
      }
      if ("count" %in% colnames(tcr_data)) {
        tcr_data$freq_in <- tcr_data$count / sum(tcr_data$count) * 100
      } else if ("count..reads." %in% colnames(tcr_data)) {
        tcr_data$freq_in <- tcr_data$count..reads. / sum(tcr_data$count..reads.) * 100
      } else if ("count..templates.reads." %in% colnames(tcr_data)) {
        tcr_data$freq_in <- tcr_data$count..templates.reads. / sum(tcr_data$count..templates.reads.) * 100
      } else {
        stop(paste0("invalid count format in file: ", basename(x)))
      }
      tcr_data <- tcr_data %>% filter(freq_in >= min_freq)
      tcr_data$time <- tools::file_path_sans_ext(basename(x))
      if ("count" %in% colnames(tcr_data)) {
        tcr_data <- subset(tcr_data, select = c("nucleotide", "aminoAcid", "time", "count", "freq_in", "vGeneName", "jGeneName"))
      } else if ("count..reads." %in% colnames(tcr_data)) {
        tcr_data <- subset(tcr_data, select = c("nucleotide", "aminoAcid", "time", "count..reads.", "freq_in", "vGeneName", "jGeneName"))
      } else if ("count..templates.reads." %in% colnames(tcr_data)) {
        tcr_data <- subset(tcr_data, select = c("nucleotide", "aminoAcid", "time", "count..templates.reads.", "freq_in", "vGeneName", "jGeneName"))
      }
      tcr_data$nucleotide <- paste0(tcr_data$nucleotide, "_", tcr_data$vGeneName, "_", tcr_data$jGeneName)
      tcr_data$vGeneName <- NULL
      tcr_data$jGeneName <- NULL
      colnames(tcr_data) <- c("nucleotide", "aminoAcid", "time", "count", "freq_in")
      return(tcr_data)
    })
    names(tcr_data_i) <- tools::file_path_sans_ext(basename(files))
    output_list[[i]] <- tcr_data_i
    message(paste0(grep(i, samples), "/", length(samples), " ", i))
  }
  assign(output.list.name, output_list, envir = .GlobalEnv)
  message("✅ done!")
}

plot_top_trajectories <- function(patient_list, time = 0, top_n = 30, floor_y = 1e-3, clonotype_by = "aminoAcid") {
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  selected_time <- tcr_format_timepoint(time)
  tp_levels <- tcr_timepoint_levels(patient_list)
  tp_days <- tcr_timepoint_days(tp_levels)
  if (!selected_time %in% tp_levels) {
    stop("time = ", selected_time, " does not exist in patient_list. Available timepoints: ", paste(tp_levels, collapse = ", "))
  }
  all_data <- bind_rows(patient_list, .id = "timepoint") %>%
    mutate(timepoint = factor(timepoint, levels = tp_levels),
           day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
           clone_id = as.character(.data[[clonotype_by]])) %>%
    filter(!is.na(clone_id),
           clone_id != "")
  top_clones_all <- all_data %>%
    filter(as.character(timepoint) == selected_time,
           freq_in > 0) %>%
    group_by(clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(freq_in), desc(count))
  top_clones <- top_clones_all %>%
    slice_head(n = top_n) %>%
    mutate(rank_selected = row_number(),
           clone_label = as.character(rank_selected)) %>%
    select(clone_id, rank_selected, clone_label)
  if (nrow(top_clones) == 0) {
    stop("No clonotypes found for ", selected_time, ".")
  }
  plot_data <- all_data %>%
    semi_join(top_clones, by = "clone_id") %>%
    group_by(timepoint, day, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop") %>%
    right_join(expand_grid(timepoint = factor(tp_levels, levels = tp_levels),
                           clone_id = top_clones$clone_id),
               by = c("timepoint", "clone_id")) %>%
    mutate(day = as.numeric(gsub("[^0-9]", "", as.character(timepoint)))) %>%
    left_join(top_clones, by = "clone_id") %>%
    mutate(count = replace_na(count, 0),
           freq_in = replace_na(freq_in, 0),
           freq_pct = pmax(freq_in, floor_y),
           clone_label = factor(clone_label, levels = as.character(seq_len(nrow(top_clones)))))
  clonotype_label <- tcr_clonotype_label(clonotype_by)
  n_plotted <- length(unique(plot_data$clone_id))
  n_total <- nrow(top_clones_all)
  n_label <- ifelse(n_plotted < n_total,
                    paste0("shown = ", n_plotted, " / ", n_total),
                    paste0("n = ", n_total))
  ggplot(plot_data, aes(x = day, y = freq_pct, group = clone_id)) +
    geom_vline(xintercept = tp_days,
               linetype = "dashed",
               linewidth = 0.25,
               colour = "grey85") +
    geom_line(linewidth = 0.25,
              alpha = 0.75,
              colour = "black") +
    tcr_y_scale_log(floor_y = floor_y) +
    scale_x_continuous(breaks = tp_days) +
    labs(x = NULL,
         y = paste0("Frequency of unique\nCDR3β ", clonotype_label, " seq (%)"),
         title = paste0("Top-", n_plotted, " CDR3β ", clonotype_label, " trajectories at ", selected_time, " | ", n_label)) +
    tcr_theme_trajectory()
}

plot_newly_detected_trajectories <- function(patient_list, target_time = 42, top_n = 30, floor_y = 1e-3, clonotype_by = "aminoAcid", skip_empty = TRUE) {
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  selected_time <- tcr_format_timepoint(target_time)
  tp_levels <- tcr_timepoint_levels(patient_list)
  tp_days <- tcr_timepoint_days(tp_levels)
  if (!selected_time %in% tp_levels) {
    stop("target_time = ", selected_time, " does not exist in patient_list. Available timepoints: ", paste(tp_levels, collapse = ", "))
  }
  selected_day <- tp_days[match(selected_time, tp_levels)]
  previous_timepoints <- tp_levels[tp_days < selected_day]
  all_data <- bind_rows(patient_list, .id = "timepoint") %>%
    mutate(timepoint = factor(timepoint, levels = tp_levels),
           day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
           clone_id = as.character(.data[[clonotype_by]])) %>%
    filter(!is.na(clone_id),
           clone_id != "")
  clone_freqs <- all_data %>%
    group_by(timepoint, day, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop")
  previous_clones <- clone_freqs %>%
    filter(as.character(timepoint) %in% previous_timepoints,
           freq_in > 0) %>%
    distinct(clone_id)
  target_clones_all <- clone_freqs %>%
    filter(as.character(timepoint) == selected_time,
           freq_in > 0) %>%
    anti_join(previous_clones, by = "clone_id") %>%
    arrange(desc(freq_in), desc(count)) %>%
    distinct(clone_id, .keep_all = TRUE)
  n_total_new <- nrow(target_clones_all)
  target_clones <- target_clones_all
  if (!is.null(top_n)) {
    target_clones <- target_clones %>%
      slice_head(n = top_n)
  }
  if (nrow(target_clones) == 0) {
    if (skip_empty) {
      message("Skipping ", selected_time, ": no newly detected clonotypes found.")
      return(NULL)
    } else {
      stop("No newly detected clonotypes found for ", selected_time, ".")
    }
  }
  target_clones <- target_clones %>%
    mutate(rank_target = row_number(),
           clone_label = as.character(rank_target)) %>%
    select(clone_id, rank_target, clone_label)
  plot_data <- clone_freqs %>%
    semi_join(target_clones, by = "clone_id") %>%
    right_join(expand_grid(timepoint = factor(tp_levels, levels = tp_levels),
                           clone_id = target_clones$clone_id),
               by = c("timepoint", "clone_id")) %>%
    mutate(day = as.numeric(gsub("[^0-9]", "", as.character(timepoint)))) %>%
    left_join(target_clones, by = "clone_id") %>%
    mutate(count = replace_na(count, 0),
           freq_in = replace_na(freq_in, 0),
           freq_pct = pmax(freq_in, floor_y),
           clone_label = factor(clone_label, levels = as.character(seq_len(nrow(target_clones)))))
  clonotype_label <- tcr_clonotype_label(clonotype_by)
  n_plotted <- length(unique(plot_data$clone_id))
  n_label <- ifelse(!is.null(top_n) && n_plotted < n_total_new,
                    paste0("shown = ", n_plotted, " / ", n_total_new),
                    paste0("n = ", n_total_new))
  title_prefix <- ifelse(!is.null(top_n),
                         paste0("Top-", n_plotted, " newly detected"),
                         "Newly detected")
  ggplot(plot_data, aes(x = day, y = freq_pct, group = clone_id)) +
    geom_vline(xintercept = tp_days,
               linetype = "dashed",
               linewidth = 0.25,
               colour = "grey85") +
    geom_line(linewidth = 0.25,
              alpha = 0.75,
              colour = "black") +
    tcr_y_scale_log(floor_y = floor_y) +
    scale_x_continuous(breaks = tp_days) +
    labs(x = NULL,
         y = paste0("Frequency of unique\nCDR3β ", clonotype_label, " seq (%)"),
         title = paste0(title_prefix, " CDR3β ", clonotype_label, " trajectories at ", selected_time, " | ", n_label)) +
    tcr_theme_trajectory()
}

get_newly_detected_fc_all <- function(patient_list, floor_y = 1e-3, clonotype_by = "aminoAcid", min_count = 2) {
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  tp_levels <- names(patient_list)
  tp_days <- as.numeric(gsub("[^0-9]", "", tp_levels))
  tp_levels <- tp_levels[order(tp_days)]
  tp_days <- tp_days[order(tp_days)]
  all_data <- bind_rows(patient_list, .id = "timepoint") %>%
    mutate(timepoint = factor(timepoint, levels = tp_levels),
           day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
           clone_id = .data[[clonotype_by]]) %>%
    filter(!is.na(clone_id))
  clone_freqs <- all_data %>%
    group_by(timepoint, day, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop")
  fc_all <- bind_rows(lapply(seq_along(tp_levels), function(i) {
    target_time <- tp_levels[i]
    target_day <- tp_days[i]
    previous_timepoints <- tp_levels[tp_days < target_day]
    if (length(previous_timepoints) == 0) {
      return(NULL)
    }
    previous_clones <- clone_freqs %>%
      filter(as.character(timepoint) %in% previous_timepoints,
             freq_in > 0) %>%
      distinct(clone_id)
    clone_freqs %>%
      filter(as.character(timepoint) == target_time,
             freq_in > 0,
             count >= min_count) %>%
      anti_join(previous_clones, by = "clone_id") %>%
      mutate(prev_freq = 0,
             fc_vs_previous_floor = freq_in / floor_y,
             log2_fc_vs_previous_floor = log2(fc_vs_previous_floor),
             target_time = as.character(timepoint),
             target_day = day)
  }))
  if (nrow(fc_all) == 0) {
    stop("No newly detected clonotypes found.")
  }
  fc_all %>%
    arrange(desc(fc_vs_previous_floor), desc(freq_in), desc(count)) %>%
    mutate(global_rank = row_number(),
           clonotype_by = clonotype_by,
           min_count = min_count,
           floor_y = floor_y)
}

plot_newly_detected_fc_rank_cutoff <- function(cutoff_result,
                                                fc_col = "fc_vs_previous_floor",
                                                point_size = 0.18,
                                                selected_point_size = 2.8,
                                                shape_by_time = TRUE,
                                                include_pie = TRUE,
                                                pie_width = 0.36) {
  required_cutoff_elements <- c(
    "data",
    "cutoff",
    "knee_rank",
    "n_selected"
  )

  missing_cutoff_elements <- setdiff(
    required_cutoff_elements,
    names(cutoff_result)
  )

  if (length(missing_cutoff_elements) > 0) {
    stop(
      "cutoff_result is missing: ",
      paste(
        missing_cutoff_elements,
        collapse = ", "
      )
    )
  }

  required_columns <- c(
    "rank",
    "clonotype_by",
    "min_count",
    "target_time",
    "target_day",
    fc_col
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(cutoff_result$data)
  )

  if (length(missing_columns) > 0) {
    stop(
      "cutoff_result$data is missing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  if (
    length(cutoff_result$cutoff) != 1 ||
    is.na(cutoff_result$cutoff) ||
    !is.finite(cutoff_result$cutoff) ||
    cutoff_result$cutoff <= 0
  ) {
    stop(
      "cutoff_result$cutoff must be one positive finite value."
    )
  }

  plot_data <- cutoff_result$data %>%
    dplyr::mutate(
      significant_fc = .data[[fc_col]] >= cutoff_result$cutoff
    ) %>%
    dplyr::filter(
      !is.na(.data[[fc_col]]),
      is.finite(.data[[fc_col]]),
      .data[[fc_col]] > 0,
      !is.na(rank),
      is.finite(rank)
    )

  if (nrow(plot_data) == 0) {
    stop(
      "No valid positive finite values are available for ",
      fc_col,
      "."
    )
  }

  clonotype_by <- unique(plot_data$clonotype_by)[1]
  clonotype_label <- ifelse(
    clonotype_by == "aminoAcid",
    "AA",
    "nt"
  )
  min_count <- unique(plot_data$min_count)[1]

  selected_data <- plot_data %>%
    dplyr::filter(significant_fc) %>%
    dplyr::arrange(target_day)

  detected_levels <- selected_data %>%
    dplyr::arrange(target_day) %>%
    dplyr::distinct(target_time) %>%
    dplyr::pull(target_time) %>%
    as.character()

  selected_data <- selected_data %>%
    dplyr::mutate(
      target_time = factor(
        target_time,
        levels = detected_levels
      )
    )

  shape_lookup <- c(
    d14 = 16,
    d42 = 17,
    d56 = 8,
    d84 = 15,
    d112 = 18,
    d140 = 25,
    d154 = 3,
    d166 = 4,
    d182 = 5,
    d196 = 6,
    d210 = 0
  )

  color_lookup <- c(
    d14 = "#668CF6",
    d42 = "#EC6A2D",
    d56 = "#33A66B",
    d84 = "#F4B43F",
    d112 = "#705CE6",
    d140 = "#CC3B7B",
    d154 = "#2DB7B5",
    d166 = "#8E9A3A",
    d182 = "#8C6D62",
    d196 = "#4C78A8",
    d210 = "#7A7A7A"
  )

  fallback_shapes <- c(
    16,
    17,
    8,
    15,
    18,
    25,
    3,
    4,
    5,
    6,
    0,
    1,
    2
  )

  fallback_colors <- c(
    "#668CF6",
    "#EC6A2D",
    "#33A66B",
    "#F4B43F",
    "#705CE6",
    "#CC3B7B",
    "#2DB7B5",
    "#8E9A3A",
    "#8C6D62",
    "#4C78A8",
    "#7A7A7A",
    "#A06CD5"
  )

  missing_levels <- setdiff(
    detected_levels,
    names(shape_lookup)
  )

  if (length(missing_levels) > 0) {
    shape_lookup[missing_levels] <- rep(
      fallback_shapes,
      length.out = length(missing_levels)
    )

    color_lookup[missing_levels] <- rep(
      fallback_colors,
      length.out = length(missing_levels)
    )
  }

  time_shapes <- shape_lookup[detected_levels]
  time_colors <- color_lookup[detected_levels]

  p <- ggplot(
    plot_data,
    aes(
      x = rank,
      y = .data[[fc_col]]
    )
  ) +
    geom_line(
      linewidth = 0.15,
      lineend = "round",
      colour = "gray30"
    ) +
    geom_point(
      data = plot_data %>%
        dplyr::filter(!significant_fc),
      size = point_size,
      colour = "gray30"
    )

  if (shape_by_time && nrow(selected_data) > 0) {
    p <- p +
      geom_point(
        data = selected_data,
        aes(
          shape = target_time,
          colour = target_time,
          fill = target_time
        ),
        size = selected_point_size,
        stroke = 0.25
      ) +
      scale_shape_manual(
        values = time_shapes,
        breaks = detected_levels,
        drop = TRUE,
        name = "Detected at"
      ) +
      scale_colour_manual(
        values = time_colors,
        breaks = detected_levels,
        drop = TRUE,
        name = "Detected at"
      ) +
      scale_fill_manual(
        values = time_colors,
        breaks = detected_levels,
        drop = TRUE,
        name = "Detected at"
      )
  } else {
    p <- p +
      geom_point(
        data = selected_data,
        size = selected_point_size,
        colour = "black"
      )
  }

  p_rank <- p +
    geom_vline(
      xintercept = cutoff_result$knee_rank,
      linewidth = 0.35,
      linetype = "dashed",
      colour = "gray25",
      na.rm = TRUE
    ) +
    geom_hline(
      yintercept = cutoff_result$cutoff,
      linewidth = 0.35,
      linetype = "dashed",
      colour = "gray25"
    ) +
    scale_y_log10() +
    scale_x_continuous(
      labels = scales::comma
    ) +
    labs(
      x = "Rank",
      y = "Expansion FC",
      title = paste0(
        "New CDR3β ",
        clonotype_label,
        " rank plot | min count ≥ ",
        min_count,
        " | cutoff = ",
        signif(
          cutoff_result$cutoff,
          3
        ),
        " | selected = ",
        cutoff_result$n_selected
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text = element_text(
        size = 8,
        colour = "black"
      ),
      axis.title = element_text(
        size = 10,
        colour = "black"
      ),
      axis.line = element_line(
        linewidth = 0.3,
        colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.3,
        colour = "black"
      ),
      panel.grid = element_blank(),
      plot.title = element_text(
        size = 11,
        face = "plain",
        colour = "black",
        hjust = 0.5
      ),
      legend.position = ifelse(
        shape_by_time &&
          nrow(selected_data) > 0,
        "right",
        "none"
      ),
      legend.title = element_text(
        size = 8
      ),
      legend.text = element_text(
        size = 7
      ),
      legend.key.size = grid::unit(
        0.45,
        "cm"
      ),
      aspect.ratio = 1
    )

  if (include_pie && nrow(selected_data) > 0) {
    pie_data <- selected_data %>%
      dplyr::mutate(
        target_time = as.character(target_time),
        target_day = as.numeric(
          gsub(
            "[^0-9]",
            "",
            target_time
          )
        )
      ) %>%
      dplyr::count(
        target_time,
        target_day,
        name = "n_expanded"
      ) %>%
      dplyr::arrange(target_day) %>%
      dplyr::mutate(
        target_time = factor(
          target_time,
          levels = unique(target_time)
        ),
        percent = 100 * n_expanded / sum(n_expanded),
        label = paste0(
          target_time,
          "   ",
          n_expanded,
          " (",
          round(
            percent,
            1
          ),
          "%)"
        ),
        label_y = rev(seq_len(dplyr::n()))
      )

    pie_colors <- color_lookup[
      as.character(
        pie_data$target_time
      )
    ]

    missing_pie_colors <- is.na(
      pie_colors
    )

    if (any(missing_pie_colors)) {
      pie_colors[missing_pie_colors] <- rep(
        fallback_colors,
        length.out = sum(
          missing_pie_colors
        )
      )
    }

    names(pie_colors) <- as.character(
      pie_data$target_time
    )

    p_pie <- ggplot(
      pie_data,
      aes(
        x = "",
        y = n_expanded,
        fill = target_time
      )
    ) +
      geom_col(
        width = 1,
        colour = "black",
        linewidth = 0.15
      ) +
      coord_polar(
        theta = "y"
      ) +
      scale_fill_manual(
        values = pie_colors,
        drop = FALSE
      ) +
      labs(
        title = paste0(
          "Expanded\nn = ",
          sum(
            pie_data$n_expanded
          )
        ),
        x = NULL,
        y = NULL
      ) +
      theme_void(
        base_size = 8
      ) +
      theme(
        plot.title = element_text(
          size = 9.5,
          face = "plain",
          hjust = 0.5,
          colour = "black"
        ),
        legend.position = "none",
        plot.margin = margin(
          5,
          5,
          0,
          5
        )
      )

    p_pie_labels <- ggplot(
      pie_data,
      aes(y = label_y)
    ) +
      geom_point(
        aes(
          x = 0,
          colour = target_time
        ),
        shape = 15,
        size = 2.4
      ) +
      geom_text(
        aes(
          x = 0.08,
          label = label
        ),
        hjust = 0,
        size = 2.5,
        colour = "black"
      ) +
      scale_colour_manual(
        values = pie_colors,
        guide = "none"
      ) +
      coord_cartesian(
        xlim = c(
          0,
          1.05
        ),
        ylim = c(
          0.5,
          nrow(pie_data) + 0.5
        ),
        clip = "off"
      ) +
      theme_void(
        base_size = 8
      ) +
      theme(
        plot.margin = margin(
          0,
          5,
          5,
          5
        )
      )

    p_pie_panel <- (
      p_pie /
        p_pie_labels
    ) +
      patchwork::plot_layout(
        heights = c(
          1,
          0.28
        )
      )

    p_pie_panel <- patchwork::wrap_elements(
      full = p_pie_panel
    )

    return(
      p_pie_panel +
        p_rank +
        patchwork::plot_layout(
          widths = c(
            pie_width,
            1
          )
        )
    )
  }

  p_rank
}


plot_expanded_trajectories <- function(patient_list, target_time = 42, expanded_data, top_n = NULL, floor_y = 1e-3, clonotype_by = "aminoAcid", skip_empty = TRUE) {
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  selected_time <- tcr_format_timepoint(target_time)
  tp_levels <- tcr_timepoint_levels(patient_list)
  tp_days <- tcr_timepoint_days(tp_levels)
  if (!selected_time %in% tp_levels) {
    stop("target_time = ", selected_time, " does not exist in patient_list. Available timepoints: ", paste(tp_levels, collapse = ", "))
  }
  all_data <- bind_rows(patient_list, .id = "timepoint") %>%
    mutate(timepoint = factor(timepoint, levels = tp_levels),
           day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
           clone_id = as.character(.data[[clonotype_by]])) %>%
    filter(!is.na(clone_id),
           clone_id != "")
  expanded_clones <- expanded_data
  if ("selected_by_knee" %in% colnames(expanded_clones)) {
    expanded_clones <- expanded_clones %>%
      filter(selected_by_knee)
  }
  if ("target_time" %in% colnames(expanded_clones)) {
    expanded_clones <- expanded_clones %>%
      filter(as.character(.data$target_time) == selected_time)
  }
  if ("clonotype_by" %in% colnames(expanded_clones)) {
    expanded_clones <- expanded_clones %>%
      filter(.data$clonotype_by == !!clonotype_by)
  }
  expanded_clones <- expanded_clones %>%
    arrange(desc(freq_in), desc(count)) %>%
    distinct(clone_id, .keep_all = TRUE)
  if (!is.null(top_n)) {
    expanded_clones <- expanded_clones %>%
      slice_head(n = top_n)
  }
  if (nrow(expanded_clones) == 0) {
    if (skip_empty) {
      message("Skipping ", selected_time, ": no expanded clonotypes found.")
      return(NULL)
    } else {
      stop("No expanded clonotypes found for ", selected_time, ".")
    }
  }
  expanded_clones <- expanded_clones %>%
    mutate(rank_target = row_number(),
           clone_label = as.character(rank_target)) %>%
    select(clone_id, rank_target, clone_label)
  plot_data <- all_data %>%
    semi_join(expanded_clones, by = "clone_id") %>%
    group_by(timepoint, day, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop") %>%
    right_join(expand_grid(timepoint = factor(tp_levels, levels = tp_levels),
                           clone_id = expanded_clones$clone_id),
               by = c("timepoint", "clone_id")) %>%
    mutate(day = as.numeric(gsub("[^0-9]", "", as.character(timepoint)))) %>%
    left_join(expanded_clones, by = "clone_id") %>%
    mutate(count = replace_na(count, 0),
           freq_in = replace_na(freq_in, 0),
           freq_pct = pmax(freq_in, floor_y),
           clone_label = factor(clone_label, levels = as.character(seq_len(nrow(expanded_clones)))))
  clonotype_label <- tcr_clonotype_label(clonotype_by)
  n_total <- length(unique(plot_data$clone_id))
  ggplot(plot_data, aes(x = day, y = freq_pct, group = clone_id)) +
    geom_vline(xintercept = tp_days,
               linetype = "dashed",
               linewidth = 0.25,
               colour = "grey85") +
    geom_line(linewidth = 0.25,
              alpha = 1,
              colour = "black") +
    tcr_y_scale_log(floor_y = floor_y) +
    scale_x_continuous(breaks = tp_days) +
    labs(x = NULL,
         y = paste0("Frequency of unique\nCDR3β ", clonotype_label, " seq (%)"),
         title = paste0("Expanded newly detected CDR3β ", clonotype_label, " trajectories at ", selected_time, " | n = ", n_total)) +
    tcr_theme_trajectory()
}

get_clonotype_set <- function(x, clonotype_by = "aminoAcid", min_count = 1, time_filter = NULL) {
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  if (!is.null(time_filter)) {
    selected_time <- tcr_format_timepoint(time_filter)
  }
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(time_filter)) {
      if (!selected_time %in% names(x)) {
        stop("time_filter = ", selected_time, " does not exist. Available timepoints: ", paste(names(x), collapse = ", "))
      }
      x <- x[selected_time]
    }
    dat <- bind_rows(x, .id = "source_id")
  } else {
    dat <- as_tibble(x) %>%
      mutate(source_id = "sample")
    if (!is.null(time_filter)) {
      if ("timepoint" %in% colnames(dat)) {
        dat <- dat %>% filter(as.character(timepoint) == selected_time)
      } else if ("time" %in% colnames(dat)) {
        dat <- dat %>% filter(as.character(time) == selected_time)
      } else {
        stop("time_filter was provided, but no 'timepoint' or 'time' column was found.")
      }
    }
  }
  if (!clonotype_by %in% colnames(dat)) {
    stop("Column ", clonotype_by, " not found.")
  }
  if (!"count" %in% colnames(dat)) {
    dat$count <- 1
  }
  if (!"freq_in" %in% colnames(dat)) {
    dat$freq_in <- NA_real_
  }
  dat %>%
    mutate(clone_id = as.character(.data[[clonotype_by]])) %>%
    filter(!is.na(clone_id),
           clone_id != "") %>%
    group_by(source_id, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = ifelse(all(is.na(freq_in)), NA_real_, sum(freq_in, na.rm = TRUE)),
              .groups = "drop") %>%
    filter(count >= min_count,
           is.na(freq_in) | freq_in > 0) %>%
    distinct(clone_id) %>%
    pull(clone_id)
}

get_pbl_tumor_overlap <- function(pbl_patient, biopsy_patient, clonotype_by = "aminoAcid", pbl_time = NULL, min_count_pbl = 1, min_count_tumor = 1) {
  pbl_clones <- get_clonotype_set(pbl_patient,
                                  clonotype_by = clonotype_by,
                                  min_count = min_count_pbl,
                                  time_filter = pbl_time)
  tumor_clones <- get_clonotype_set(biopsy_patient,
                                    clonotype_by = clonotype_by,
                                    min_count = min_count_tumor,
                                    time_filter = NULL)
  overlap_clones <- intersect(pbl_clones, tumor_clones)
  list(
    sets = list(PBL = pbl_clones,
                Tumor = tumor_clones),
    overlap = tibble(clone_id = overlap_clones),
    summary = tibble(category = c("PBL only", "Tumor only", "PBL ∩ Tumor"),
                     n = c(length(setdiff(pbl_clones, tumor_clones)),
                           length(setdiff(tumor_clones, pbl_clones)),
                           length(overlap_clones))),
    clonotype_by = clonotype_by,
    pbl_time = pbl_time,
    min_count_pbl = min_count_pbl,
    min_count_tumor = min_count_tumor
  )
}

plot_pbl_tumor_venn <- function(overlap_result) {
  if (!requireNamespace("ggvenn", quietly = TRUE)) {
    stop("Package 'ggvenn' is required. Install it with install.packages('ggvenn').")
  }
  clonotype_label <- tcr_clonotype_label(overlap_result$clonotype_by)
  pbl_label <- ifelse(is.null(overlap_result$pbl_time),
                      "PBL",
                      paste0("PBL ", tcr_format_timepoint(overlap_result$pbl_time)))
  sets <- overlap_result$sets
  names(sets)[names(sets) == "PBL"] <- pbl_label
  set_colors <- c("#668CF6", "#CC3B7B")
  names(set_colors) <- names(sets)
  ggvenn::ggvenn(sets,
                 fill_color = set_colors,
                 stroke_color = "black",
                 stroke_size = 0.25,
                 set_name_size = 3,
                 text_size = 3,
                 show_percentage = TRUE) +
    labs(title = paste0(pbl_label, " vs tumor CDR3β ", clonotype_label, " overlap")) +
    theme(plot.title = element_text(size = 10, face = "plain", hjust = 0.5, colour = "black"),
          legend.position = "none")
}

plot_overlapped_trajectories <- function(patient_list, overlap_result, top_n = 30, floor_y = 1e-3, clonotype_by = NULL, skip_empty = TRUE) {
  if (is.null(clonotype_by)) {
    clonotype_by <- overlap_result$clonotype_by
  }
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  overlap_ids <- overlap_result$overlap %>%
    distinct(clone_id) %>%
    pull(clone_id)
  tp_levels <- tcr_timepoint_levels(patient_list)
  tp_days <- tcr_timepoint_days(tp_levels)
  all_data <- bind_rows(patient_list, .id = "timepoint") %>%
    mutate(timepoint = factor(timepoint, levels = tp_levels),
           day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
           clone_id = as.character(.data[[clonotype_by]])) %>%
    filter(!is.na(clone_id),
           clone_id != "")
  clone_freqs <- all_data %>%
    group_by(timepoint, day, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop")
  overlapped_clones_all <- clone_freqs %>%
    filter(clone_id %in% overlap_ids,
           freq_in > 0) %>%
    group_by(clone_id) %>%
    summarise(max_freq = max(freq_in, na.rm = TRUE),
              total_count = sum(count, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(max_freq), desc(total_count))
  n_total_overlap <- nrow(overlapped_clones_all)
  overlapped_clones <- overlapped_clones_all
  if (!is.null(top_n)) {
    overlapped_clones <- overlapped_clones %>%
      slice_head(n = top_n)
  }
  if (nrow(overlapped_clones) == 0) {
    if (skip_empty) {
      message("No overlapped clonotypes found.")
      return(NULL)
    } else {
      stop("No overlapped clonotypes found.")
    }
  }
  overlapped_clones <- overlapped_clones %>%
    mutate(rank_overlap = row_number(),
           clone_label = as.character(rank_overlap)) %>%
    select(clone_id, rank_overlap, clone_label)
  plot_data <- clone_freqs %>%
    semi_join(overlapped_clones, by = "clone_id") %>%
    right_join(expand_grid(timepoint = factor(tp_levels, levels = tp_levels),
                           clone_id = overlapped_clones$clone_id),
               by = c("timepoint", "clone_id")) %>%
    mutate(day = as.numeric(gsub("[^0-9]", "", as.character(timepoint)))) %>%
    left_join(overlapped_clones, by = "clone_id") %>%
    mutate(count = replace_na(count, 0),
           freq_in = replace_na(freq_in, 0),
           freq_pct = pmax(freq_in, floor_y),
           clone_label = factor(clone_label, levels = as.character(seq_len(nrow(overlapped_clones)))))
  clonotype_label <- tcr_clonotype_label(clonotype_by)
  pbl_label <- ifelse(is.null(overlap_result$pbl_time),
                      "PBL",
                      paste0("PBL ", tcr_format_timepoint(overlap_result$pbl_time)))
  n_plotted <- length(unique(plot_data$clone_id))
  n_label <- ifelse(!is.null(top_n) && n_plotted < n_total_overlap,
                    paste0("shown = ", n_plotted, " / ", n_total_overlap),
                    paste0("n = ", n_total_overlap))
  ggplot(plot_data, aes(x = day, y = freq_pct, group = clone_id)) +
    geom_vline(xintercept = tp_days,
               linetype = "dashed",
               linewidth = 0.25,
               colour = "grey85") +
    geom_line(linewidth = 0.25,
              alpha = 0.75,
              colour = "black") +
    tcr_y_scale_log(floor_y = floor_y) +
    scale_x_continuous(breaks = tp_days) +
    labs(x = NULL,
         y = paste0("Frequency of unique\nCDR3β ", clonotype_label, " seq (%)"),
         title = paste0("Top-", n_plotted, " ", pbl_label, "–tumor overlapping CDR3β ", clonotype_label, " trajectories | ", n_label)) +
    tcr_theme_trajectory()
}

plot_pbl_tumor_overlap_panel <- function(pbl_patient, biopsy_patient, clonotype_by = "aminoAcid", pbl_time = NULL, top_n = 30, floor_y = 1e-3, min_count_pbl = 1, min_count_tumor = 1) {
  overlap_result <- get_pbl_tumor_overlap(pbl_patient,
                                          biopsy_patient,
                                          clonotype_by = clonotype_by,
                                          pbl_time = pbl_time,
                                          min_count_pbl = min_count_pbl,
                                          min_count_tumor = min_count_tumor)
  p_traj <- plot_overlapped_trajectories(pbl_patient,
                                         overlap_result,
                                         top_n = top_n,
                                         floor_y = floor_y,
                                         clonotype_by = clonotype_by)
  p_venn <- plot_pbl_tumor_venn(overlap_result)
  plot <- if (is.null(p_traj)) {
    p_venn
  } else {
    p_traj + p_venn + patchwork::plot_layout(widths = c(1, 1))
  }
  list(overlap = overlap_result,
       summary = overlap_result$summary,
       p_traj = p_traj,
       p_venn = p_venn,
       plot = plot)
}

plot_non_neutral_rankplot_zoom <- function(cutoff_result, fc_col = "fc_vs_previous_floor", point_size = 2.4, zoom_factor = 1.15) {
  plot_data <- cutoff_result$data %>%
    mutate(non_neutral = .data[[fc_col]] >= cutoff_result$cutoff)
  selected_data <- plot_data %>%
    filter(non_neutral) %>%
    arrange(target_day)
  detected_levels <- selected_data %>%
    arrange(target_day) %>%
    distinct(target_time) %>%
    pull(target_time) %>%
    as.character()
  selected_data <- selected_data %>%
    mutate(target_time = factor(target_time, levels = detected_levels))
  shape_lookup <- c(
    d14 = 16,
    d42 = 17,
    d56 = 8,
    d84 = 15,
    d112 = 18,
    d140 = 25,
    d154 = 3,
    d166 = 4,
    d182 = 5,
    d196 = 6,
    d210 = 0
  )
  color_lookup <- c(
    d14 = "#668CF6",
    d42 = "#EC6A2D",
    d56 = "#33A66B",
    d84 = "#F4B43F",
    d112 = "#705CE6",
    d140 = "#CC3B7B",
    d154 = "#2DB7B5",
    d166 = "#8E9A3A",
    d182 = "#8C6D62",
    d196 = "#4C78A8",
    d210 = "#7A7A7A"
  )
  fallback_shapes <- c(16, 17, 8, 15, 18, 25, 3, 4, 5, 6, 0, 1, 2)
  fallback_colors <- c("#668CF6", "#EC6A2D", "#33A66B", "#F4B43F", "#705CE6", "#CC3B7B", "#2DB7B5", "#8E9A3A", "#8C6D62", "#4C78A8", "#7A7A7A", "#A06CD5")
  missing_levels <- setdiff(detected_levels, names(shape_lookup))
  if (length(missing_levels) > 0) {
    shape_lookup[missing_levels] <- rep(fallback_shapes, length.out = length(missing_levels))
    color_lookup[missing_levels] <- rep(fallback_colors, length.out = length(missing_levels))
  }
  time_shapes <- shape_lookup[detected_levels]
  time_colors <- color_lookup[detected_levels]
  ggplot(plot_data, aes(x = rank, y = .data[[fc_col]])) +
    geom_line(linewidth = 0.15, colour = "gray30") +
    geom_point(data = plot_data %>% filter(!non_neutral),
               size = 0.15,
               colour = "gray30") +
    geom_point(data = selected_data,
               aes(shape = target_time, colour = target_time, fill = target_time),
               size = point_size,
               stroke = 0.25) +
    scale_shape_manual(values = time_shapes,
                       breaks = detected_levels,
                       drop = TRUE,
                       name = "Detected at") +
    scale_colour_manual(values = time_colors,
                        breaks = detected_levels,
                        drop = TRUE,
                        name = "Detected at") +
    scale_fill_manual(values = time_colors,
                      breaks = detected_levels,
                      drop = TRUE,
                      name = "Detected at") +
    geom_hline(yintercept = cutoff_result$cutoff,
               linewidth = 0.35,
               linetype = "dashed",
               colour = "gray25") +
    geom_vline(xintercept = cutoff_result$knee_rank,
               linewidth = 0.35,
               linetype = "dashed",
               colour = "gray25") +
    coord_cartesian(xlim = c(0, cutoff_result$knee_rank * zoom_factor)) +
    scale_y_log10() +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Rank",
         y = "Expansion FC",
         title = "Zoom on non-neutral clonotypes") +
    theme_classic(base_size = 9) +
    theme(axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 11, face = "plain", colour = "black", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.size = grid::unit(0.45, "cm"),
          aspect.ratio = 1)
}

summarise_overlap_ci_by_time <- function(patient_list, overlap_result, clonotype_by = NULL, top_n = NULL, include_zeros = TRUE) {
  if (is.null(clonotype_by)) {
    clonotype_by <- overlap_result$clonotype_by
  }
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  overlap_ids <- overlap_result$overlap %>%
    distinct(clone_id) %>%
    pull(clone_id)
  tp_levels <- names(patient_list)
  tp_days <- as.numeric(gsub("[^0-9]", "", tp_levels))
  tp_levels <- tp_levels[order(tp_days)]
  tp_days <- tp_days[order(tp_days)]
  all_data <- bind_rows(patient_list, .id = "timepoint") %>%
    mutate(timepoint = factor(timepoint, levels = tp_levels),
           day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
           clone_id = as.character(.data[[clonotype_by]])) %>%
    filter(!is.na(clone_id),
           clone_id != "")
  clone_freqs <- all_data %>%
    group_by(timepoint, day, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop")
  overlap_rank <- clone_freqs %>%
    filter(clone_id %in% overlap_ids,
           freq_in > 0) %>%
    group_by(clone_id) %>%
    summarise(max_freq = max(freq_in, na.rm = TRUE),
              total_count = sum(count, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(max_freq), desc(total_count))
  if (!is.null(top_n)) {
    overlap_rank <- overlap_rank %>%
      slice_head(n = top_n)
  }
  selected_ids <- overlap_rank$clone_id
  if (include_zeros) {
    ci_data <- clone_freqs %>%
      filter(clone_id %in% selected_ids) %>%
      right_join(expand_grid(timepoint = factor(tp_levels, levels = tp_levels),
                             clone_id = selected_ids),
                 by = c("timepoint", "clone_id")) %>%
      mutate(day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
             count = replace_na(count, 0),
             freq_in = replace_na(freq_in, 0))
  } else {
    ci_data <- clone_freqs %>%
      filter(clone_id %in% selected_ids,
             freq_in > 0)
  }
  ci_data %>%
    group_by(timepoint, day) %>%
    summarise(n_clonotypes = n_distinct(clone_id),
              n_detected = sum(freq_in > 0),
              mean_freq_pct = mean(freq_in, na.rm = TRUE),
              sd_freq_pct = sd(freq_in, na.rm = TRUE),
              se_freq_pct = sd_freq_pct / sqrt(n_clonotypes),
              ci95_low = mean_freq_pct - qt(0.975, df = n_clonotypes - 1) * se_freq_pct,
              ci95_high = mean_freq_pct + qt(0.975, df = n_clonotypes - 1) * se_freq_pct,
              median_freq_pct = median(freq_in, na.rm = TRUE),
              q25_freq_pct = quantile(freq_in, 0.25, na.rm = TRUE),
              q75_freq_pct = quantile(freq_in, 0.75, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(ci95_low = pmax(ci95_low, 0),
           clonotype_by = clonotype_by,
           included_zeros = include_zeros,
           top_n = ifelse(is.null(top_n), NA_integer_, top_n),
           .before = 1) %>%
    arrange(day)
}

summarise_high_freq_overlapping_clones <- function(patient_list, overlap_result, clonotype_by = "aminoAcid", threshold = 0.1) {
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  if (is.data.frame(overlap_result)) {
    overlap_ids <- overlap_result %>%
      distinct(clone_id) %>%
      pull(clone_id)
  } else {
    overlap_ids <- overlap_result$overlap %>%
      distinct(clone_id) %>%
      pull(clone_id)
  }
  tp_levels <- names(patient_list)
  tp_days <- as.numeric(gsub("[^0-9]", "", tp_levels))
  tp_levels <- tp_levels[order(tp_days)]
  tp_days <- tp_days[order(tp_days)]
  clone_time_freqs <- bind_rows(patient_list, .id = "timepoint") %>%
    mutate(timepoint = factor(timepoint, levels = tp_levels),
           day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
           clone_id = as.character(.data[[clonotype_by]])) %>%
    filter(!is.na(clone_id),
           clone_id != "",
           clone_id %in% overlap_ids) %>%
    group_by(timepoint, day, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop")
  high_freq_clones <- clone_time_freqs %>%
    group_by(clone_id) %>%
    summarise(max_freq_pct = max(freq_in, na.rm = TRUE),
              max_timepoint = paste(as.character(timepoint)[freq_in == max(freq_in, na.rm = TRUE)], collapse = ";"),
              n_timepoints_above_threshold = sum(freq_in > threshold, na.rm = TRUE),
              timepoints_above_threshold = paste(as.character(timepoint)[freq_in > threshold], collapse = ";"),
              max_count = max(count, na.rm = TRUE),
              .groups = "drop") %>%
    filter(max_freq_pct > threshold) %>%
    arrange(desc(max_freq_pct))
  list(
    summary = tibble(
      threshold_pct = threshold,
      n_overlapping_clonotypes = length(unique(overlap_ids)),
      n_above_threshold_at_least_once = nrow(high_freq_clones),
      max_freq_pct = max(high_freq_clones$max_freq_pct, na.rm = TRUE),
      max_clone_id = high_freq_clones$clone_id[which.max(high_freq_clones$max_freq_pct)],
      max_timepoint = high_freq_clones$max_timepoint[which.max(high_freq_clones$max_freq_pct)]
    ),
    clones = high_freq_clones,
    all_timepoint_freqs = clone_time_freqs
  )
}

plot_tested_tcr_trajectories <- function(patient_list,
                                         tcr_lib,
                                         reactivity_col,
                                         color_col,
                                         title,
                                         clonotype_by = "aminoAcid",
                                         floor_y = 1e-3,
                                         include_reactivities = c("Yes", "Auto-reactive", "No"),
                                         skip_empty = TRUE,
                                         show_legend = FALSE,
                                         facet_col = NULL) {
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  if (!reactivity_col %in% colnames(tcr_lib)) {
    stop("Column ", reactivity_col, " not found in tcr_lib.")
  }
  if (!color_col %in% colnames(tcr_lib)) {
    stop("Column ", color_col, " not found in tcr_lib.")
  }
  if (!is.null(facet_col) && !facet_col %in% colnames(tcr_lib)) {
    stop("Column ", facet_col, " not found in tcr_lib.")
  }
  facet_cols <- if (is.null(facet_col)) character(0) else facet_col
  tp_levels <- tcr_timepoint_levels(patient_list)
  tp_days <- tcr_timepoint_days(tp_levels)
  selected_clones <- tcr_lib %>%
    mutate(clone_id = sub("_.*$", "", V1),
           reactivity = .data[[reactivity_col]],
           plot_color = .data[[color_col]]) %>%
    filter(!is.na(clone_id),
           clone_id != "",
           !is.na(reactivity),
           reactivity %in% include_reactivities,
           !is.na(plot_color)) %>%
    distinct(across(all_of(c("clone_id", "reactivity", "plot_color", "V2", facet_cols))))
  if (nrow(selected_clones) == 0) {
    if (skip_empty) {
      message("No tested TCRs found for ", title, ".")
      return(NULL)
    } else {
      stop("No tested TCRs found for ", title, ".")
    }
  }
  selected_clones <- selected_clones %>%
    group_by(clone_id) %>%
    summarise(TCR_ID = paste(unique(V2), collapse = "; "),
              reactivity = dplyr::first(reactivity),
              plot_color = dplyr::first(plot_color),
              across(all_of(facet_cols), ~ dplyr::first(.x)),
              .groups = "drop")
  all_data <- bind_rows(patient_list, .id = "timepoint") %>%
    mutate(timepoint = factor(timepoint, levels = tp_levels),
           day = as.numeric(gsub("[^0-9]", "", as.character(timepoint))),
           clone_id = as.character(.data[[clonotype_by]])) %>%
    filter(!is.na(clone_id),
           clone_id != "")
  clone_freqs <- all_data %>%
    group_by(timepoint, day, clone_id) %>%
    summarise(count = sum(count, na.rm = TRUE),
              freq_in = sum(freq_in, na.rm = TRUE),
              .groups = "drop")
  plot_data <- clone_freqs %>%
    semi_join(selected_clones, by = "clone_id") %>%
    right_join(expand_grid(timepoint = factor(tp_levels, levels = tp_levels),
                           clone_id = selected_clones$clone_id),
               by = c("timepoint", "clone_id")) %>%
    mutate(day = as.numeric(gsub("[^0-9]", "", as.character(timepoint)))) %>%
    left_join(selected_clones, by = "clone_id") %>%
    mutate(count = replace_na(count, 0),
           freq_in = replace_na(freq_in, 0),
           freq_pct = pmax(freq_in, floor_y),
           reactivity = factor(reactivity, levels = include_reactivities))
  color_values <- selected_clones %>%
    distinct(reactivity, plot_color) %>%
    filter(!is.na(reactivity),
           !is.na(plot_color))
  color_values <- stats::setNames(color_values$plot_color, color_values$reactivity)
  clonotype_label <- tcr_clonotype_label(clonotype_by)
  n_total <- length(unique(plot_data$clone_id))
  p <- ggplot(plot_data, aes(x = day, y = freq_pct, group = clone_id, colour = reactivity)) +
    geom_vline(xintercept = tp_days,
               linetype = "dashed",
               linewidth = 0.25,
               colour = "grey85") +
    geom_line(linewidth = 0.25,
              alpha = 0.85) +
    scale_colour_manual(values = color_values,
                        breaks = include_reactivities,
                        drop = TRUE,
                        name = NULL) +
    tcr_y_scale_log(floor_y = floor_y) +
    scale_x_continuous(breaks = tp_days) +
    labs(x = NULL,
         y = paste0("Frequency of unique\nCDR3β ", clonotype_label, " seq (%)"),
         title = paste0(title, " | n = ", n_total)) +
    tcr_theme_trajectory()
  if (!is.null(facet_col)) {
    p <- p +
      facet_wrap(stats::as.formula(paste("~", facet_col)))
  }
  if (show_legend) {
    p <- p +
      theme(legend.position = "right",
            legend.title = element_blank(),
            legend.text = element_text(size = 7),
            legend.key.size = grid::unit(0.35, "cm"))
  }
  p
}

apply_newly_detected_fc_cutoff <- function(fc_all, fc_cutoff = 10, fc_col = "fc_vs_previous_floor") {
  plot_data <- fc_all %>%
    filter(!is.na(.data[[fc_col]]),
           .data[[fc_col]] > 0) %>%
    arrange(desc(.data[[fc_col]]), desc(freq_in), desc(count)) %>%
    mutate(rank = row_number(),
           selected_by_cutoff = .data[[fc_col]] >= fc_cutoff,
           selected_by_knee = selected_by_cutoff)
  cutoff_rank <- plot_data %>%
    filter(selected_by_cutoff) %>%
    summarise(cutoff_rank = max(rank, na.rm = TRUE)) %>%
    pull(cutoff_rank)
  if (length(cutoff_rank) == 0 || is.infinite(cutoff_rank)) {
    cutoff_rank <- NA_real_
  }
  list(cutoff = fc_cutoff,
       knee_rank = cutoff_rank,
       n_selected = sum(plot_data$selected_by_cutoff, na.rm = TRUE),
       data = plot_data)
}

make_sc_trb_long <- function(sc_tcr) {
  sc_tcr %>%
    mutate(cell_row = row_number()) %>%
    separate_rows(v_gene, j_gene, cdr3, chain, sep = ";") %>%
    mutate(
      across(c(v_gene, j_gene, cdr3, chain), str_trim),
      CDR3B = cdr3,
      trb_clonotype = paste(CDR3B, v_gene, j_gene, sep = "_")
    ) %>%
    filter(
      chain == "TRB",
      !is.na(CDR3B),
      CDR3B != "",
      CDR3B != "NA"
    ) %>%
    distinct(cell_row, barcode, sample, CDR3B, v_gene, j_gene, trb_clonotype, .keep_all = TRUE)
}

collapse_uniq <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != "" & x != "NA"]
  if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")
}

get_sc_reference_fc_rank <- function(patient_list,
                                     single_cell_clonotypes,
                                     floor_y = 1e-3,
                                     clonotype_by = "aminoAcid",
                                     timepoint_0 = "d0",
                                     timepoint_42 = "d42") {
  if (!clonotype_by %in% c("aminoAcid", "nucleotide")) {
    stop("clonotype_by must be 'aminoAcid' or 'nucleotide'")
  }
  if (!"CDR3B" %in% colnames(single_cell_clonotypes)) {
    stop("single_cell_clonotypes must contain a CDR3B column.")
  }
  if (!timepoint_0 %in% names(patient_list)) {
    stop("timepoint_0 = ", timepoint_0, " not found. Available timepoints: ", paste(names(patient_list), collapse = ", "))
  }
  if (!timepoint_42 %in% names(patient_list)) {
    stop("timepoint_42 = ", timepoint_42, " not found. Available timepoints: ", paste(names(patient_list), collapse = ", "))
  }
  sc_input <- single_cell_clonotypes %>%
    mutate(CDR3B = as.character(CDR3B)) %>%
    filter(!is.na(CDR3B), CDR3B != "") %>%
    distinct(CDR3B, .keep_all = TRUE)
  clone_freqs <- bind_rows(patient_list[c(timepoint_0, timepoint_42)], .id = "timepoint") %>%
    mutate(
      clone_raw = as.character(.data[[clonotype_by]]),
      CDR3B = str_extract(clone_raw, "^[^_]+")
    ) %>%
    filter(!is.na(CDR3B), CDR3B != "", CDR3B %in% sc_input$CDR3B) %>%
    group_by(CDR3B, timepoint) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      freq_in = sum(freq_in, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    complete(
      CDR3B = sc_input$CDR3B,
      timepoint = c(timepoint_0, timepoint_42),
      fill = list(count = 0, freq_in = 0)
    )
  fc_rank <- clone_freqs %>%
    select(CDR3B, timepoint, count, freq_in) %>%
    pivot_wider(
      names_from = timepoint,
      values_from = c(count, freq_in),
      values_fill = list(count = 0, freq_in = 0)
    ) %>%
    mutate(
      reference_count = .data[[paste0("count_", timepoint_0)]],
      reference_freq = .data[[paste0("freq_in_", timepoint_0)]],
      count = .data[[paste0("count_", timepoint_42)]],
      freq_in = .data[[paste0("freq_in_", timepoint_42)]],
      fc_signed = pmax(freq_in, floor_y) / pmax(reference_freq, floor_y),
      log2FC_signed = log2(fc_signed),
      FC_abs = 2^abs(log2FC_signed),
      log2FC_abs = abs(log2FC_signed),
      FC_direction = case_when(
        log2FC_signed > 0 ~ "Expansion",
        log2FC_signed < 0 ~ "Contraction",
        TRUE ~ "No change"
      ),
      timepoint = timepoint_42,
      reference_timepoint = timepoint_0
    ) %>%
    left_join(sc_input, by = "CDR3B") %>%
    arrange(desc(FC_abs), desc(freq_in), desc(count), CDR3B) %>%
    mutate(
      rank = row_number(),
      floor_y = floor_y
    )
  fc_rank
}

plot_sc_clonotype_fc_rank <- function(fc_rank,
                                      point_size = 1.8,
                                      label_size = 2.3,
                                      fc_cutoff = 10) {
  plot_data <- fc_rank %>%
    mutate(
      TCRlib_reactivity_plot = case_when(
        !is.na(TCRlib_reactivity) & TCRlib_reactivity != "" ~ TCRlib_reactivity,
        TCRlib_hit %in% TRUE ~ "TCRlib: unannotated",
        TRUE ~ "Not in TCRlib"
      ),
      highlight_yes = str_detect(TCRlib_reactivity_plot, regex("(^| \\| )Yes($| \\| )")),
      above_fc_cutoff = FC_abs >= fc_cutoff
    )
  
  label_data <- plot_data %>%
    filter(highlight_yes) %>%
    arrange(desc(FC_abs))
  
  cutoff_rank <- if (any(plot_data$above_fc_cutoff, na.rm = TRUE)) {
    max(plot_data$rank[plot_data$above_fc_cutoff], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  p <- ggplot(plot_data, aes(x = rank, y = FC_abs)) +
    geom_line(linewidth = 0.15, colour = "gray35") +
    geom_point(data = plot_data %>% filter(!highlight_yes),
               aes(shape = FC_direction),
               colour = "grey80",
               size = point_size * 0.7,
               alpha = 0.45,
               stroke = 0.15) +
    geom_point(data = plot_data %>% filter(above_fc_cutoff, !highlight_yes),
               aes(shape = FC_direction),
               colour = "black",
               fill = "black",
               size = point_size * 1.1,
               alpha = 0.85,
               stroke = 0.25) +
    geom_point(data = plot_data %>% filter(highlight_yes),
               aes(shape = FC_direction),
               colour = "red3",
               fill = "red3",
               size = point_size * 1.5,
               alpha = 1,
               stroke = 0.35) +
    ggrepel::geom_text_repel(data = label_data,
                             aes(label = CDR3B),
                             colour = "red3",
                             size = label_size,
                             max.overlaps = Inf,
                             min.segment.length = 0,
                             linewidth = 0.15,
                             show.legend = FALSE) +
    scale_y_log10() +
    geom_hline(yintercept = fc_cutoff,
               linetype = "dashed",
               linewidth = 0.35,
               colour = "gray50") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Rank",
         y = "Absolute FC vs max(baseline, d0)",
         shape = "Direction",
         title = paste0("Single-cell clonotypes | FC ≥ ", fc_cutoff,
                        " | n = ", nrow(plot_data),
                        " | selected = ", sum(plot_data$above_fc_cutoff, na.rm = TRUE),
                        " | Yes = ", nrow(label_data))) +
    theme_classic(base_size = 9) +
    theme(axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 11, face = "plain", colour = "black", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.size = grid::unit(0.45, "cm"),
          aspect.ratio = 1)
  if (!is.na(cutoff_rank)) {
    p <- p +
      geom_vline(xintercept = cutoff_rank,
                 linetype = "dashed",
                 linewidth = 0.25,
                 colour = "gray50")
  }
  p
}

criteria_colors <- c(
  "None" = "gray70",
  "Frequent" = "#9ECAE1",
  "Expanding" = "#F16913",
  "Tumor" = "#CB181D",
  "PD1" = "#6A51A3",
  "Tumor | Expanding" = "#A63603",
  "PD1 | Expanding" = "#C51B8A",
  "PD1 | Tumor" = "#54278F",
  "PD1 | Tumor | Expanding" = "#3F007D"
)

plot_sc_clonotype_fc_rank_tcrlib_criteria <- function(fc_rank,
                                                      tcr_lib_annotated,
                                                      fc_cutoff = 10,
                                                      point_size = 1.8,
                                                      label_size = 2.3,
                                                      label_tcrlib = FALSE) {
  collapse_uniq <- function(x) {
    x <- unique(na.omit(as.character(x)))
    x <- x[x != "" & x != "NA"]
    if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")
  }
  
  criteria_levels <- names(criteria_colors)
  
  tcr_lib_annotated <- tcr_lib_annotated %>%
    as_tibble()
  
  if (!"CDR3B" %in% colnames(tcr_lib_annotated)) {
    tcr_lib_annotated <- tcr_lib_annotated %>%
      mutate(CDR3B = stringr::str_extract(V1, "^[^_]+"))
  }
  
  criteria_by_CDR3B <- tcr_lib_annotated %>%
    filter(!is.na(CDR3B), CDR3B != "",
           !is.na(V2), V2 != "") %>%
    group_by(CDR3B) %>%
    summarise(
      TCRlib_ID_plot = collapse_uniq(V2),
      is_PD1 = any(stringr::str_detect(criteria_simplified, stringr::regex("(^| \\| )PD1($| \\| )")), na.rm = TRUE),
      is_tumor = any(stringr::str_detect(criteria_simplified, stringr::regex("(^| \\| )Tumor($| \\| )")), na.rm = TRUE),
      is_expanding = any(stringr::str_detect(criteria_simplified, stringr::regex("(^| \\| )Expanding($| \\| )")), na.rm = TRUE),
      is_frequent = any(criteria_simplified == "Frequent", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      criteria_simplified = purrr::pmap_chr(
        list(is_PD1, is_tumor, is_expanding, is_frequent),
        function(is_PD1, is_tumor, is_expanding, is_frequent) {
          x <- c(
            if (is_PD1) "PD1" else character(0),
            if (is_tumor) "Tumor" else character(0),
            if (is_expanding) "Expanding" else character(0)
          )
          if (length(x) > 0) {
            paste(x, collapse = " | ")
          } else if (is_frequent) {
            "Frequent"
          } else {
            "None"
          }
        }
      )
    ) %>%
    select(CDR3B, TCRlib_ID_plot, criteria_simplified)
  
  plot_data <- fc_rank %>%
    mutate(CDR3B = as.character(CDR3B)) %>%
    left_join(criteria_by_CDR3B, by = "CDR3B") %>%
    mutate(
      in_TCRlib_plot = !is.na(TCRlib_ID_plot),
      criteria_simplified = if_else(
        in_TCRlib_plot,
        criteria_simplified,
        "Not in TCRlib"
      ),
      criteria_simplified = if_else(
        is.na(criteria_simplified) | criteria_simplified == "Not found in single-cell data",
        "Not in TCRlib",
        criteria_simplified
      ),
      criteria_simplified = factor(criteria_simplified, levels = criteria_levels),
      above_fc_cutoff = FC_abs >= fc_cutoff
    )
  
  label_data <- plot_data %>%
    filter(in_TCRlib_plot, label_tcrlib) %>%
    arrange(desc(FC_abs))
  
  cutoff_rank <- plot_data %>%
    filter(above_fc_cutoff) %>%
    summarise(cutoff_rank = max(rank, na.rm = TRUE)) %>%
    pull(cutoff_rank)
  
  if (length(cutoff_rank) == 0 || is.infinite(cutoff_rank)) {
    cutoff_rank <- NA_real_
  }
  
  ggplot(plot_data, aes(x = rank, y = FC_abs)) +
    geom_line(linewidth = 0.15, colour = "gray35") +
    geom_point(data = plot_data %>% filter(!in_TCRlib_plot, !above_fc_cutoff),
               aes(shape = FC_direction),
               colour = "grey88",
               size = point_size * 0.7,
               alpha = 0.45,
               stroke = 0.15) +
    geom_point(data = plot_data %>% filter(!in_TCRlib_plot, above_fc_cutoff),
               aes(shape = FC_direction),
               colour = "black",
               fill = "black",
               size = point_size * 1.0,
               alpha = 0.75,
               stroke = 0.25) +
    geom_point(data = plot_data %>% filter(in_TCRlib_plot),
               aes(shape = FC_direction, colour = criteria_simplified, fill = criteria_simplified),
               size = point_size * 1.45,
               alpha = 0.95,
               stroke = 0.25) +
    ggrepel::geom_text_repel(data = label_data,
                             aes(label = TCRlib_ID_plot, colour = criteria_simplified),
                             size = label_size,
                             max.overlaps = Inf,
                             min.segment.length = 0,
                             linewidth = 0.15,
                             show.legend = FALSE) +
    scale_colour_manual(values = criteria_colors,
                        breaks = criteria_levels,
                        drop = FALSE,
                        name = "Simplified criteria") +
    scale_fill_manual(values = criteria_colors,
                      breaks = criteria_levels,
                      drop = FALSE,
                      name = "Simplified criteria") +
    scale_y_log10() +
    geom_hline(yintercept = fc_cutoff,
               linetype = "dashed",
               linewidth = 0.35,
               colour = "gray50") +
    {if (!is.na(cutoff_rank)) geom_vline(xintercept = cutoff_rank,
                                         linetype = "dashed",
                                         linewidth = 0.25,
                                         colour = "gray50")} +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Ranked single-cell CDR3b clonotypes",
         y = "Fold-change (day 42 vs baseline)",
         shape = "Dynamics",
         title = paste0("Single-cell clonotypes | FC ≥ ", fc_cutoff,
                        " | ranked CDR3B = ", nrow(plot_data),
                        " | FC10 CDR3B = ", sum(plot_data$above_fc_cutoff, na.rm = TRUE),
                        " | TCRlib CDR3B = ", sum(plot_data$in_TCRlib_plot, na.rm = TRUE))) +
    tcr_theme_rank()
}
plot_annotated_bulk_trajectories <- function(patient_list,
                                             annotation,
                                             time,
                                             manual_colors,
                                             compartment,
                                             output_dir,
                                             floor_y = 1e-3,
                                             ncol = 7,
                                             width = 16,
                                             height = 4.5,
                                             dpi = 600) {
  required_annotation_cols <- c(
    "key",
    "TCRb",
    "manual.annotation"
  )
  missing_annotation_cols <- setdiff(
    required_annotation_cols,
    colnames(annotation)
  )
  if (length(missing_annotation_cols) > 0) {
    stop(
      "annotation is missing: ",
      paste(missing_annotation_cols, collapse = ", ")
    )
  }
  required_bulk_cols <- c(
    "aminoAcid",
    "count",
    "freq_in"
  )
  missing_bulk_cols <- setdiff(
    required_bulk_cols,
    unique(unlist(lapply(patient_list, colnames)))
  )
  if (length(missing_bulk_cols) > 0) {
    stop(
      "patient_list is missing: ",
      paste(missing_bulk_cols, collapse = ", ")
    )
  }
  tp_levels <- tcr_timepoint_levels(patient_list)
  tp_days <- tcr_timepoint_days(tp_levels)
  annotation_map_all <- annotation %>%
    select(
      key,
      TCRb,
      manual.annotation
    ) %>%
    filter(
      !is.na(TCRb),
      TCRb != "",
      !is.na(manual.annotation)
    ) %>%
    tidyr::separate_rows(
      TCRb,
      sep = ";"
    ) %>%
    mutate(
      TCRb = stringr::str_trim(TCRb),
      aminoAcid = stringr::str_remove(TCRb, "_.*$")
    ) %>%
    filter(
      !is.na(aminoAcid),
      aminoAcid != ""
    ) %>%
    distinct(
      key,
      TCRb,
      aminoAcid,
      manual.annotation
    )
  ambiguous_annotation <- annotation_map_all %>%
    group_by(aminoAcid) %>%
    filter(
      n_distinct(manual.annotation) > 1
    ) %>%
    ungroup() %>%
    arrange(aminoAcid, manual.annotation)
  annotation_map <- annotation_map_all %>%
    group_by(aminoAcid) %>%
    filter(
      n_distinct(manual.annotation) == 1
    ) %>%
    summarise(
      key = paste(unique(key), collapse = ";"),
      TCRb = paste(unique(TCRb), collapse = ";"),
      manual.annotation = first(manual.annotation),
      .groups = "drop"
    )
  bulk_freq <- bind_rows(
    patient_list,
    .id = "timepoint"
  ) %>%
    transmute(
      timepoint = factor(
        timepoint,
        levels = tp_levels
      ),
      day = as.numeric(
        stringr::str_remove_all(
          as.character(timepoint),
          "[^0-9]"
        )
      ),
      aminoAcid = as.character(aminoAcid),
      count,
      freq_in
    ) %>%
    filter(
      !is.na(aminoAcid),
      aminoAcid != ""
    ) %>%
    group_by(
      timepoint,
      day,
      aminoAcid
    ) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      freq_in = sum(freq_in, na.rm = TRUE),
      .groups = "drop"
    )
  matched_annotation <- annotation_map %>%
    semi_join(
      bulk_freq,
      by = "aminoAcid"
    )
  unmatched_annotation <- annotation_map %>%
    anti_join(
      bulk_freq,
      by = "aminoAcid"
    )
  if (nrow(matched_annotation) == 0) {
    stop(
      "No annotated amino-acid clonotypes were found in patient_list."
    )
  }
  plot_data <- bulk_freq %>%
    semi_join(
      matched_annotation,
      by = "aminoAcid"
    ) %>%
    select(
      timepoint,
      day,
      aminoAcid,
      count,
      freq_in
    ) %>%
    right_join(
      tidyr::expand_grid(
        timepoint = factor(
          tp_levels,
          levels = tp_levels
        ),
        aminoAcid = matched_annotation$aminoAcid
      ),
      by = c(
        "timepoint",
        "aminoAcid"
      )
    ) %>%
    mutate(
      day = tcr_timepoint_days(
        as.character(timepoint)
      ),
      count = tidyr::replace_na(
        count,
        0
      ),
      freq_in = tidyr::replace_na(
        freq_in,
        0
      ),
      freq_plot = pmax(
        freq_in,
        floor_y
      )
    ) %>%
    left_join(
      matched_annotation,
      by = "aminoAcid"
    ) %>%
    arrange(
      manual.annotation,
      aminoAcid,
      day
    )
  missing_colors <- setdiff(
    unique(as.character(plot_data$manual.annotation)),
    names(manual_colors)
  )
  if (length(missing_colors) > 0) {
    stop(
      "manual_colors is missing: ",
      paste(missing_colors, collapse = ", ")
    )
  }
  p <- ggplot(
    plot_data,
    aes(
      x = day,
      y = freq_plot,
      group = aminoAcid,
      color = manual.annotation
    )
  ) +
    geom_vline(
      xintercept = tp_days,
      linetype = "dashed",
      linewidth = 0.25,
      color = "grey85"
    ) +
    geom_line(
      linewidth = 0.35,
      alpha = 0.75
    ) +
    facet_wrap(
      ~ manual.annotation,
      ncol = ncol
    ) +
    scale_x_continuous(
      breaks = tp_days,
      expand = expansion(
        mult = c(0.02, 0.02)
      )
    ) +
    scale_y_log10(
      limits = c(
        floor_y,
        10
      ),
      breaks = c(
        1e-3,
        1e-2,
        1e-1,
        1,
        10
      ),
      labels = c(
        "0.001",
        "0.01",
        "0.1",
        "1",
        "10"
      )
    ) +
    scale_color_manual(
      values = manual_colors,
      drop = FALSE
    ) +
    labs(
      x = "Time",
      y = "Clonotype frequency (%)"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(
        size = 10,
        color = "black"
      ),
      axis.text = element_text(
        size = 8,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.25,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.15,
        color = "black"
      ),
      axis.ticks.length = grid::unit(
        1.5,
        "pt"
      ),
      panel.border = element_rect(
        fill = NA,
        linewidth = 0.4,
        color = "black"
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 9,
        face = "plain",
        color = "black"
      ),
      legend.position = "none",
      panel.spacing.x = grid::unit(
        0.8,
        "lines"
      ),
      aspect.ratio = 0.5,
      plot.margin = margin(
        5,
        5,
        5,
        5
      )
    )
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  output_prefix <- paste0(
    "bulk_",
    time
  )
  write.table(
    plot_data,
    file = file.path(
      output_dir,
      paste0(output_prefix, ".tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  ggsave(
    filename = file.path(
      output_dir,
      paste0(output_prefix, ".png")
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  list(
    data = plot_data,
    annotation_map = annotation_map,
    unmatched_annotation = unmatched_annotation,
    ambiguous_annotation = ambiguous_annotation,
    plot = p
  )
}
