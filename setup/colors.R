# --- Colors for plotting ---
binary_colors <- c(
  "yes" = "#3c4b99",
  "no" = "gray80"
)
criteria_colors <- c(
  "Frequent" = "#9ECAE1",
  "Expanding" = "#F16913",
  "Tumor" = "#CB181D",
  "PD1" = "#6A51A3",
  "Tumor | Frequent" = "#FB6A4A",
  "Tumor | Expanding" = "#A63603",
  "PD1 | Frequent" = "#9E9AC8",
  "PD1 | Expanding" = "#C51B8A",
  "PD1 | Tumor" = "#54278F",
  "PD1 | Tumor | Frequent" = "#756BB1",
  "PD1 | Tumor | Expanding" = "#3F007D",
  "Non tested" = "gray80",
  "Not found" = "gray60",
  "None" = "gray20"
)
reactivity_colors <- c(
  "Yes" = "red3",
  "No" = "black",
  "Auto-reactive" = "#E69F00",
  "EBV" = "#3c4b99",
  "CMV" = "seagreen3",
  "TCR library tested" = "gray60",
  "NA" = "gray80"
)
manual_colors <- c(
  "CD8 TCM" = "#0B6F61",
  "CD8 LEF1+ ZNF683+" = "#84D8C5",
  "CD8 activated TEM GZMK+" = "#E45C95",
  "CD8 TEM KLRG1+" = "#B08BCE",
  "CD8 TEMRA FGFBP2+" = "#123FD1",
  "CD8 TEMRA LILRB1+" = "#78B7F9",
  "CD8 activated/proliferating" = "#4B1D73"
)
clone_colors <- c(
  "Single (0 < X <= 1)" = "gray80",
  "Small (1 < X <= 5)" = "#4DBBD5",
  "Medium (5 < X <= 20)" = "#3C5488",
  "Large (20 < X <= 100)" = "#F39B7F",
  "Hyperexpanded (100 < X <= 500)" = "#E64B35"
)
phenotype_colors <- c(
  "CD8.NaiveLike" = "gray80",
  "CD4.NaiveLike" = "gray80",
  "CD8.MAIT" = "#ADFF2F",
  "CD8.CM" = "#f4c40f",
  "CD8.EM" = "#FF69B4",
  "CD8.TEMRA" = "#8A2BE2",
  "CD8.TEX" = "black",
  "CD8.TPEX" = "#fe9b60",
  "CD4.CTL_EOMES" = "#32CD32",
  "CD4.CTL_Exh" = "#4A4999",
  "CD4.CTL_GNLY" = "#7FFFD4",
  "CD4.Tfh" = "#1E90FF",
  "CD4.Th17" = "#556B2F",
  "CD4.Treg" = "#d8443c"
)
signature_colors <- c(
  "Random" = "gray40",
  "TRTPBL_UCell" = "purple",
  "TRTPBL_NeoTCRPBL_UCell" = "purple4",
  "Yossef_NeoTCR_PBL_UCell" = "darkblue",
  "iTAG_pre_UCell" = "#E41A1C",
  "iTAG_pre_v2_UCell" = "red2",
  "iTAG_post1_UCell" = "#7A0019",
  "iTAG_post1_v2_UCell" = "red4",
  "iTAG_post2_UCell" = "#FF8C00",
  "iTAG_post2_v2_UCell" = "orange3",
  "iTAG_post6_UCell" = "#A6761D",
  "iTAG_all_UCell" = "green2"
)
category_colors_yossef <- c(
  "Reactive" = "darkred",
  "Non-reactive" = "gray15",
  "CD8.TEMRA" = "#8A2BE2",
  "CD8.EM" = "#FF69B4",
  "CD8.TEX" = "black",
  "Non-hits" = "#8C6D5C",
  "NaiveLike" = "#7A8F6B",
  "Other" = "#C6B8A6",
  "Viral_databases" = "#2F6F8F"
)

# --- Configure plotting style and visualization settings ---
theme_ic <- theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 0.5,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "plain"),
    plot.title = element_text(hjust = 0),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  )
