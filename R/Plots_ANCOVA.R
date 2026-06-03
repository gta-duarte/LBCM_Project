#!/usr/bin/env Rscript
# =========================================================

BASE_DIR   <- "~/project_inescid/"
CSV_DIR    <- file.path(BASE_DIR, "output/csv")
DATA_DIR   <- file.path(BASE_DIR, "data")

PLOT_DIR_BIN <- file.path(BASE_DIR, "output/plots_ancova_cnv_clusters_BINARY")
REP_DIR_BIN  <- file.path(BASE_DIR, "output/reports_ancova_cnv_clusters_BINARY")

PLOT_DIR_MUL <- file.path(BASE_DIR, "output/plots_ancova_cnv_clusters_MULTILEVEL")
REP_DIR_MUL  <- file.path(BASE_DIR, "output/reports_ancova_cnv_clusters_MULTILEVEL")

dir.create(PLOT_DIR_BIN, recursive = TRUE, showWarnings = FALSE)
dir.create(REP_DIR_BIN,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR_MUL, recursive = TRUE, showWarnings = FALSE)
dir.create(REP_DIR_MUL,  recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(e1071)
  library(tibble)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# -----------------------------
# Parameters
# -----------------------------
FDR_THRESHOLD <- 0.05
TOP_HITS_BOX  <- 30
LABEL_MAX     <- 50

CNV_ORDER        <- c("Deletion", "Loss", "Neutral", "Gain", "Amplification")
CNV_BINARY_ORDER <- c("Neutral", "Altered")

SIGNED_LINE <- 0.01
DELTA_LINE  <- 0.2
SHAPIRO_MAX_N <- 5000

TISSUE_VAR <- "OncotreeLineage"

# -----------------------------
# Input paths
# -----------------------------
ancova_bin_path <- file.path(CSV_DIR, "ANCOVA_featurecentric_CNV_clusters_BINARY.csv")
contr_bin_path  <- file.path(CSV_DIR, "ANCOVA_featurecentric_CNV_clusters_contrasts_BINARY.csv")
map_bin_path    <- file.path(CSV_DIR, "CNV_clusters_map_BINARY.csv")

ancova_mul_path <- file.path(CSV_DIR, "ANCOVA_featurecentric_CNV_clusters_MULTILEVEL.csv")
contr_mul_path  <- file.path(CSV_DIR, "ANCOVA_featurecentric_CNV_clusters_contrasts_MULTILEVEL.csv")
map_mul_path    <- file.path(CSV_DIR, "CNV_clusters_map_MULTILEVEL.csv")

crispr_long_path <- file.path(DATA_DIR, "CRISPR_long_all.csv")
cnv_mapped_path  <- file.path(DATA_DIR, "CNV_mapped_collapsed.csv")
model_path       <- file.path(DATA_DIR, "Model.csv")

stopifnot(file.exists(ancova_bin_path))
stopifnot(file.exists(contr_bin_path))
stopifnot(file.exists(map_bin_path))
stopifnot(file.exists(ancova_mul_path))
stopifnot(file.exists(contr_mul_path))
stopifnot(file.exists(map_mul_path))
stopifnot(file.exists(crispr_long_path))
stopifnot(file.exists(cnv_mapped_path))
stopifnot(file.exists(model_path))

# -----------------------------
# Load common files
# -----------------------------
crispr_long <- read_csv(crispr_long_path, show_col_types = FALSE)
cnv_mapped  <- read_csv(cnv_mapped_path, show_col_types = FALSE)
model       <- read_csv(model_path, show_col_types = FALSE)

crispr_long <- crispr_long %>%
  mutate(
    ModelID = as.character(ModelID),
    geneY = as.character(geneY),
    essentiality = as.numeric(essentiality)
  ) %>%
  filter(is.finite(essentiality))

tissue_map <- model %>%
  select(ModelID, !!sym(TISSUE_VAR)) %>%
  rename(tissue = !!sym(TISSUE_VAR)) %>%
  filter(!is.na(tissue), tissue != "") %>%
  mutate(
    ModelID = as.character(ModelID),
    tissue  = factor(tissue)
  ) %>%
  distinct()

cnv_mapped <- cnv_mapped %>%
  mutate(
    ModelID = as.character(ModelID),
    geneX = as.character(geneX),
    cn_category = as.character(cn_category)
  ) %>%
  distinct()

# =========================================================
# Helpers
# =========================================================
safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  if (length(x) > SHAPIRO_MAX_N) x <- sample(x, SHAPIRO_MAX_N)
  shapiro.test(x)$p.value
}

add_labels_if_available <- function(p, label_df) {
  if (has_ggrepel && nrow(label_df) > 0) {
    p + ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 3
    )
  } else {
    p
  }
}

prepare_cluster_meta <- function(cluster_map) {
  cluster_id_col <- intersect(c("cluster_id", "geneX"), names(cluster_map))
  if (length(cluster_id_col) == 0) stop("cluster_map missing cluster_id/geneX column")
  if (!("gene_list" %in% names(cluster_map))) stop("cluster_map missing gene_list column")

  cluster_map %>%
    transmute(
      geneX = as.character(.data[[cluster_id_col[1]]]),
      gene_list = as.character(gene_list)
    ) %>%
    distinct()
}

prepare_ancova_tbl <- function(ancova_tbl, cluster_meta) {
  req <- c("geneX", "geneY", "p_value", "signed_eta2", "FDR_by_geneX", "rep_gene", "n_genes")
  miss <- setdiff(req, names(ancova_tbl))
  if (length(miss) > 0) stop("ANCOVA table missing columns: ", paste(miss, collapse = ", "))

  ancova_tbl %>%
    mutate(
      geneX = as.character(geneX),
      geneY = as.character(geneY),
      p_value = as.numeric(p_value),
      signed_eta2 = as.numeric(signed_eta2),
      FDR_by_geneX = as.numeric(FDR_by_geneX),
      rep_gene = as.character(rep_gene),
      n_genes = as.integer(n_genes),
      neglog10_p = -log10(p_value),
      is_hit = is.finite(FDR_by_geneX) & (FDR_by_geneX <= FDR_THRESHOLD)
    ) %>%
    mutate(gene_list = as.character(gene_list))
}

make_rank_tbl <- function(ancova_tbl) {
  ancova_tbl %>%
    group_by(geneX) %>%
    summarise(
      n_tests   = n(),
      n_hits    = sum(is_hit, na.rm = TRUE),
      best_FDR  = suppressWarnings(min(FDR_by_geneX, na.rm = TRUE)),
      best_p    = suppressWarnings(min(p_value, na.rm = TRUE)),
      top_geneY = geneY[which.min(p_value)],
      max_abs_delta = suppressWarnings(max(abs(delta_mean), na.rm = TRUE)),
      rep_gene  = dplyr::first(na.omit(rep_gene)),
      n_genes   = dplyr::first(na.omit(n_genes)),
      gene_list = dplyr::first(na.omit(gene_list)),
      .groups = "drop"
    ) %>%
    mutate(
      rep_gene  = ifelse(is.na(rep_gene)  | rep_gene  == "", "NA", rep_gene),
      gene_list = ifelse(is.na(gene_list) | gene_list == "", "NA", gene_list)
    ) %>%
    arrange(desc(n_hits), best_FDR, desc(max_abs_delta), best_p)
}

# =========================================================
# 1) MULTI-LEVEL BRANCH (MAIN)
# =========================================================
ancova_mul   <- read_csv(ancova_mul_path, show_col_types = FALSE)
contrast_mul <- read_csv(contr_mul_path, show_col_types = FALSE)
cluster_mul  <- read_csv(map_mul_path, show_col_types = FALSE)

cluster_meta_mul <- prepare_cluster_meta(cluster_mul)
ancova_mul <- prepare_ancova_tbl(ancova_mul, cluster_meta_mul)

req_mul_cols <- c("delta_mean", "contrast", "level", "contrast_p_value")
miss_mul_cols <- setdiff(req_mul_cols, names(ancova_mul))
if (length(miss_mul_cols) > 0) {
  stop("Main multilevel ANCOVA table missing columns: ", paste(miss_mul_cols, collapse = ", "))
}

ancova_mul <- ancova_mul %>%
  mutate(
    delta_mean = as.numeric(delta_mean),
    contrast = as.character(contrast),
    level = as.character(level),
    contrast_p_value = as.numeric(contrast_p_value)
  )

cnv_for_plot_mul <- ancova_mul %>%
  select(geneX, rep_gene) %>%
  distinct() %>%
  filter(!is.na(rep_gene), rep_gene != "") %>%
  inner_join(cnv_mapped, by = c("rep_gene" = "geneX")) %>%
  transmute(
    ModelID = as.character(ModelID),
    geneX = as.character(geneX),
    cn_category = factor(cn_category, levels = CNV_ORDER)
  ) %>%
  distinct()

# -----------------------------
# Main rank report (MULTI-LEVEL)
# -----------------------------
rank_mul <- make_rank_tbl(ancova_mul)

writeLines(
  c(
    "TOP RANK CLUSTERS (MULTI-LEVEL CNV branch)",
    "=========================================",
    paste0("Date: ", Sys.time()),
    paste0("Input ANCOVA: ", basename(ancova_mul_path)),
    paste0("Input contrasts: ", basename(contr_mul_path)),
    paste0("Input cluster map: ", basename(map_mul_path)),
    paste0("FDR threshold: ", FDR_THRESHOLD),
    "",
    "This is the MAIN ranking branch.",
    "Gatekeeper: FDR_by_geneX on cn_category.",
    "Main effect interpretation: best adjusted contrast vs Neutral.",
    "",
    "Top 50 clusters:",
    paste(capture.output(print(rank_mul %>% slice_head(n = 50), n = 50, width = Inf)), collapse = "\n"),
    ""
  ),
  con = file.path(REP_DIR_MUL, "top_rank_clusters_multilevel_report.txt")
)

# -----------------------------
# Multi-level QQ
# -----------------------------
pvals_mul <- ancova_mul$p_value
pvals_mul <- pvals_mul[pvals_mul > 0 & is.finite(pvals_mul)]

qq_mul <- tibble(
  obs = -log10(sort(pvals_mul)),
  exp = -log10(ppoints(length(pvals_mul)))
)

p_qq_mul <- ggplot(qq_mul, aes(exp, obs)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.5) +
  labs(
    title = "QQ plot: ANCOVA p-values (cn_category term) [CNV clusters]",
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(PLOT_DIR_MUL, "qqplot_ancova_pvalues_cnv_clusters_multilevel.png"),
  p_qq_mul, width = 6.5, height = 6, dpi = 300
)

# -----------------------------
# Multi-level volcano: signed effect
# -----------------------------
label_mul_signed <- ancova_mul %>%
  filter(is_hit) %>%
  mutate(abs_signed = abs(signed_eta2)) %>%
  arrange(FDR_by_geneX, desc(abs_signed)) %>%
  slice_head(n = LABEL_MAX) %>%
  mutate(
    label = ifelse(
      !is.na(rep_gene) & !is.na(n_genes),
      paste0(geneY, " | ", geneX, " [", rep_gene, "; n=", n_genes, "]"),
      paste0(geneY, " | ", geneX)
    )
  )

p_vol_mul_signed <- ggplot(ancova_mul, aes(signed_eta2, neglog10_p)) +
  geom_point(data = subset(ancova_mul, !is_hit), alpha = 0.25, size = 1) +
  geom_point(data = subset(ancova_mul,  is_hit), alpha = 0.9, size = 1.2) +
  geom_vline(xintercept = c(-SIGNED_LINE, SIGNED_LINE), linewidth = 0.4) +
  labs(
    title = paste0("Volcano (multi-level CNV): signed effect proxy vs -log10(p) (FDR ≤ ", FDR_THRESHOLD, ")"),
    x = "Signed effect proxy",
    y = "-log10(p_value)"
  ) +
  theme_classic(base_size = 12)

p_vol_mul_signed <- add_labels_if_available(p_vol_mul_signed, label_mul_signed)

ggsave(
  file.path(PLOT_DIR_MUL, "volcano_signed_effect_multilevel.png"),
  p_vol_mul_signed, width = 8.5, height = 6.5, dpi = 300
)

# -----------------------------
# Multi-level volcano: best delta_mean (MAIN volcano)
# -----------------------------
vol_mul_delta <- ancova_mul %>% filter(is.finite(delta_mean))

label_mul_delta <- vol_mul_delta %>%
  filter(is_hit) %>%
  mutate(abs_delta = abs(delta_mean)) %>%
  arrange(FDR_by_geneX, desc(abs_delta), contrast_p_value) %>%
  slice_head(n = LABEL_MAX) %>%
  mutate(
    label = ifelse(
      !is.na(rep_gene) & !is.na(n_genes),
      paste0(geneY, " | ", geneX, " [", rep_gene, "; n=", n_genes, "]"),
      paste0(geneY, " | ", geneX)
    )
  )

p_vol_mul_delta <- ggplot(vol_mul_delta, aes(delta_mean, neglog10_p)) +
  geom_point(data = subset(vol_mul_delta, !is_hit), alpha = 0.25, size = 1) +
  geom_point(data = subset(vol_mul_delta,  is_hit), alpha = 0.9, size = 1.2) +
  geom_vline(xintercept = c(-DELTA_LINE, DELTA_LINE), linewidth = 0.4) +
  labs(
    title = paste0("Volcano (multi-level CNV): best adjusted Δmean vs -log10(p) (FDR ≤ ", FDR_THRESHOLD, ")"),
    x = "Best adjusted delta_mean vs Neutral",
    y = "-log10(p_value)"
  ) +
  theme_classic(base_size = 12)

p_vol_mul_delta <- add_labels_if_available(p_vol_mul_delta, label_mul_delta)

ggsave(
  file.path(PLOT_DIR_MUL, "volcano_delta_mean_multilevel.png"),
  p_vol_mul_delta, width = 8.5, height = 6.5, dpi = 300
)

# -----------------------------
# Multi-level boxplots + Shapiro (MAIN)
# -----------------------------
hits_mul <- ancova_mul %>%
  filter(is_hit) %>%
  arrange(FDR_by_geneX, desc(abs(delta_mean)), contrast_p_value) %>%
  slice_head(n = TOP_HITS_BOX)

report_mul_path <- file.path(REP_DIR_MUL, "ancova_residual_normality_shapiro_report_multilevel.txt")

writeLines(
  c(
    "ANCOVA residual normality report (Shapiro-Wilk) [CNV clusters multilevel + PC1..PC10]",
    paste0("Date: ", Sys.time()),
    paste0("Hit definition: FDR_by_geneX ≤ ", FDR_THRESHOLD),
    paste0("Top hits evaluated: ", TOP_HITS_BOX),
    paste0("Diagnostic model used here: essentiality ~ cn_category + ", TISSUE_VAR),
    paste0("Input ANCOVA: ", basename(ancova_mul_path)),
    paste0("CNV for plots: rep_gene used as cluster proxy, raw multilevel CNV"),
    ""
  ),
  con = report_mul_path
)

walk(seq_len(nrow(hits_mul)), function(i) {
  gx <- hits_mul$geneX[i]
  gy <- hits_mul$geneY[i]
  rg <- hits_mul$rep_gene[i]
  ng <- hits_mul$n_genes[i]
  best_contrast <- hits_mul$contrast[i]

  rg <- ifelse(is.na(rg) | rg == "", "NA", rg)
  ng <- ifelse(is.na(ng), NA_integer_, ng)
  best_contrast <- ifelse(is.na(best_contrast) | best_contrast == "", "NA", best_contrast)

  df <- cnv_for_plot_mul %>%
    filter(geneX == gx) %>%
    inner_join(crispr_long %>% filter(geneY == gy), by = "ModelID") %>%
    inner_join(tissue_map, by = "ModelID") %>%
    filter(!is.na(cn_category), !is.na(tissue))

  if (nrow(df) < 10) return(NULL)
  if (n_distinct(df$cn_category) < 2) return(NULL)

  fit <- aov(essentiality ~ cn_category + tissue, data = df)
  r   <- residuals(fit)

  sh_p <- safe_shapiro(r)
  sk   <- e1071::skewness(r, na.rm = TRUE, type = 2)
  ku   <- e1071::kurtosis(r, na.rm = TRUE, type = 2)

  cat(
    paste(
      "Hit", i, "of", nrow(hits_mul),
      "| cluster_id:", gx,
      if (rg != "NA") paste0("| rep_gene:", rg) else "",
      if (!is.na(ng)) paste0("| n_genes:", ng) else "",
      "| geneY:", gy,
      if (best_contrast != "NA") paste0("| best_contrast:", best_contrast) else "",
      "| n:", nrow(df),
      "| Shapiro_p:", format(sh_p, digits = 4),
      "| skew:", format(sk, digits = 3),
      "| kurt:", format(ku, digits = 3)
    ),
    "\n",
    file = report_mul_path,
    append = TRUE
  )

  p_box <- ggplot(df, aes(cn_category, essentiality)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.18, alpha = 0.35, size = 0.9) +
    labs(
      title = paste0(
        gy, " dependency by CNV cluster ", gx,
        if (rg != "NA") paste0(" [rep=", rg, "]") else "",
        if (!is.na(ng)) paste0(" [n_genes=", ng, "]") else "",
        " [best=", best_contrast, "; rep_gene raw multilevel CNV; multilevel ANCOVA used for p-values]"
      ),
      x = "CNV category",
      y = "CRISPR gene effect"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(
      PLOT_DIR_MUL,
      paste0(
        "boxplot_",
        str_replace_all(gx, "[^A-Za-z0-9]+", "_"), "_",
        str_replace_all(gy, "[^A-Za-z0-9]+", "_"),
        ".png"
      )
    ),
    p_box, width = 8, height = 5, dpi = 300
  )
})

write_csv(
  hits_mul,
  file.path(PLOT_DIR_MUL, "hit_list_for_plots_multilevel.csv")
)

