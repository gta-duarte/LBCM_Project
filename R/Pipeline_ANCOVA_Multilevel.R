#!/usr/bin/env Rscript
# =========================================================

BASE_DIR <- "~/project_inescid/"
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "output/csv/")
REP_DIR  <- file.path(BASE_DIR, "output/reports/")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REP_DIR, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(parallel)
library(purrr)
library(emmeans)
library(tibble)
library(broom)
library(e1071)
library(cluster)

N_CORES <- 1

# =========================================================
# 0) Params
# =========================================================
# K_FEATURES <- 100
MIN_SAMPLES_FEATURE <- 200
MIN_GROUP_N         <- 10
MIN_Y_N             <- 300
MIN_Y_SKEW          <- 1.0

CNV_LEVELS_KEEP <- c("Loss", "Neutral", "Gain", "Amplification", "Deletion")
MEDIA_LEVELS    <- c("RPMI", "DMEM", "MEM", "IMDM", "F12", "OTHER", "unknown")

TISSUE_VAR  <- "OncotreeLineage"
SEX_VAR     <- "Sex"
CULTURE_VAR <- "GrowthPattern"

# ---------------------------------------------------------
# Optional gene filter
# ---------------------------------------------------------
REMOVE_SEX_LINKED_GENES <- TRUE

# ---------------------------------------------------------
# CNV clustering params
# ---------------------------------------------------------
CLUSTER_MATCH_THRESHOLD <- 0.95
MIN_CLUSTER_OVERLAP     <- 100
HC_METHOD               <- "complete"
CLUSTER_CUT_HEIGHT      <- 1 - CLUSTER_MATCH_THRESHOLD

# ---------------------------------------------------------
# PCA params
# ---------------------------------------------------------
PCA_PATH     <- file.path(DATA_DIR, "gexp_PCA_scores_top10.csv")
PCA_PC_NAMES <- paste0("PC", 1:5)

# ---------------------------------------------------------
# Gene annotation
# ---------------------------------------------------------
GENE_ANNOT_PATH <- file.path(DATA_DIR, "Gene.csv")

# ---------------------------------------------------------
# Outputs
# ---------------------------------------------------------
OUT_GLOBAL   <- file.path(OUT_DIR, "ANCOVA_featurecentric_CNV_clusters_MULTILEVEL.csv")
OUT_CONTR    <- file.path(OUT_DIR, "ANCOVA_featurecentric_CNV_clusters_contrasts_MULTILEVEL.csv")
OUT_SELECTED <- file.path(OUT_DIR, "CNV_selected_clusters_MULTILEVEL.csv")
OUT_MAP      <- file.path(OUT_DIR, "CNV_clusters_map_MULTILEVEL.csv")
OUT_REPORT   <- file.path(REP_DIR, "CNV_clusters_report_MULTILEVEL.txt")
OUT_TOP50 <- file.path(OUT_DIR, "CNV_top50_ranked_clusters_MULTILEVEL.csv")

# =========================================================
# 1) Load
# =========================================================
crispr <- read_csv(file.path(DATA_DIR, "CRISPRGeneEffect.csv"), show_col_types = FALSE)
cnv    <- read_csv(file.path(DATA_DIR, "WES_pureCN_CNV_genes_20250207.csv"), show_col_types = FALSE)
model  <- read_csv(file.path(DATA_DIR, "Model.csv"), show_col_types = FALSE)
gene_annot_raw <- read_csv(GENE_ANNOT_PATH, show_col_types = FALSE)
pca_raw <- read_csv(PCA_PATH, show_col_types = FALSE)

# =========================================================
# 2) Gene annotation
# =========================================================
symbol_col <- intersect(c("symbol", "gene_symbol", "GeneSymbol", "hgnc_symbol"), names(gene_annot_raw))
if (length(symbol_col) == 0) stop("Could not find gene symbol column in Gene.csv")
symbol_col <- symbol_col[1]

locus_col <- intersect(c("locus_group", "LocusGroup"), names(gene_annot_raw))
if (length(locus_col) == 0) stop("Could not find locus_group column in Gene.csv")
locus_col <- locus_col[1]

essential_col <- intersect(c("essentiality", "Essentiality"), names(gene_annot_raw))
if (length(essential_col) == 0) stop("Could not find essentiality column in Gene.csv")
essential_col <- essential_col[1]

location_col <- intersect(c("location", "Location"), names(gene_annot_raw))
if (length(location_col) == 0) stop("Could not find location column in Gene.csv")
location_col <- location_col[1]

gene_annot <- gene_annot_raw %>%
  transmute(
    gene = as.character(.data[[symbol_col]]),
    locus_group = as.character(.data[[locus_col]]),
    essentiality_class = as.character(.data[[essential_col]]),
    location = as.character(.data[[location_col]])
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct() %>%
  mutate(
    is_sex_linked = str_detect(location, "^(X|Y)[pqPQ]?")
  )

protein_coding_genes <- gene_annot %>%
  filter(locus_group == "protein-coding gene") %>%
  pull(gene) %>%
  unique()

non_commonessential_genes <- gene_annot %>%
  filter(is.na(essentiality_class) | essentiality_class != "common essential") %>%
  pull(gene) %>%
  unique()

sex_linked_genes <- gene_annot %>%
  filter(is_sex_linked) %>%
  pull(gene) %>%
  unique()

if (REMOVE_SEX_LINKED_GENES) {
  protein_coding_genes <- setdiff(protein_coding_genes, sex_linked_genes)
  non_commonessential_genes <- setdiff(non_commonessential_genes, sex_linked_genes)
}

# =========================================================
# 3) SIDM -> ACH (1:1)
# =========================================================
id_map <- model %>%
  select(ModelID, SangerModelID) %>%
  filter(!is.na(ModelID), !is.na(SangerModelID), SangerModelID != "") %>%
  distinct() %>%
  group_by(SangerModelID) %>%
  slice(1) %>%
  ungroup()

# =========================================================
# 3.1) Covariate maps
# =========================================================
tissue_map <- model %>%
  select(ModelID, !!sym(TISSUE_VAR)) %>%
  rename(tissue = !!sym(TISSUE_VAR)) %>%
  filter(!is.na(tissue), tissue != "") %>%
  mutate(tissue = factor(tissue)) %>%
  distinct()

sex_map <- model %>%
  select(ModelID, !!sym(SEX_VAR)) %>%
  rename(sex_raw = !!sym(SEX_VAR)) %>%
  mutate(
    sex = case_when(
      is.na(sex_raw) | sex_raw == "" ~ "Unknown",
      sex_raw %in% c("Male", "Female") ~ sex_raw,
      TRUE ~ "Unknown"
    ),
    sex = factor(sex, levels = c("Male", "Female", "Unknown"))
  ) %>%
  select(ModelID, sex) %>%
  distinct()

culture_map <- model %>%
  select(ModelID, !!sym(CULTURE_VAR)) %>%
  rename(culture_raw = !!sym(CULTURE_VAR)) %>%
  mutate(
    culture = case_when(
      is.na(culture_raw) | culture_raw == "" ~ "unknown",
      str_to_lower(culture_raw) %in% c("adherent") ~ "adherent",
      str_to_lower(culture_raw) %in% c("suspension") ~ "suspension",
      TRUE ~ "unknown"
    ),
    culture = factor(culture, levels = c("adherent", "suspension", "unknown"))
  ) %>%
  select(ModelID, culture) %>%
  distinct()

classify_base <- function(x) {
  if (is.na(x) || str_trim(x) == "") return("unknown")
  x_up <- str_to_upper(x)

  hits <- c(
    RPMI = str_detect(x_up, "\\bRPMI\\b"),
    DMEM = str_detect(x_up, "\\bDMEM\\b"),
    MEM  = str_detect(x_up, "\\bMEM\\b") | str_detect(x_up, "\\bEMEM\\b"),
    IMDM = str_detect(x_up, "\\bIMDM\\b"),
    F12  = str_detect(x_up, "\\bF12\\b") | str_detect(x_up, "F-12K")
  )

  if (sum(hits) == 1) return(names(hits)[which(hits)])
  "OTHER"
}

media_map <- model %>%
  transmute(
    ModelID,
    media_base = map_chr(as.character(FormulationID), classify_base)
  ) %>%
  mutate(media_base = factor(media_base, levels = MEDIA_LEVELS)) %>%
  distinct()

# =========================================================
# 3.2) PCA covariates
# =========================================================
required_pca_cols <- c("SangerModelID", "ModelID", PCA_PC_NAMES)
missing_pca_cols <- setdiff(required_pca_cols, names(pca_raw))
if (length(missing_pca_cols) > 0) {
  stop("Missing PCA columns: ", paste(missing_pca_cols, collapse = ", "))
}

pca_map <- pca_raw %>%
  transmute(
    SangerModelID = as.character(SangerModelID),
    ModelID       = as.character(ModelID),
    across(all_of(PCA_PC_NAMES), as.numeric)
  ) %>%
  filter(!is.na(ModelID), ModelID != "") %>%
  distinct()

# =========================================================
# 4) CRISPR long
# =========================================================
crispr_long <- crispr %>%
  rename(ModelID = `...1`) %>%
  pivot_longer(cols = -ModelID, names_to = "gene_raw", values_to = "essentiality") %>%
  mutate(
    geneY = str_extract(gene_raw, "^[^ ]+"),
    essentiality = suppressWarnings(as.numeric(essentiality))
  ) %>%
  select(ModelID, geneY, essentiality) %>%
  filter(geneY %in% protein_coding_genes) %>%
  filter(geneY %in% non_commonessential_genes) %>%
  group_by(geneY) %>%
  filter(
    sum(is.finite(essentiality)) >= MIN_Y_N,
    abs(e1071::skewness(essentiality, na.rm = TRUE, type = 2)) >= MIN_Y_SKEW
  ) %>%
  ungroup()

# =========================================================
# 5) CNV mapped + collapsed
# =========================================================
dist_tbl <- tibble(
  cn_category = CNV_LEVELS_KEEP,
  dist_neutral = case_when(
    cn_category == "Neutral" ~ 0L,
    cn_category %in% c("Loss", "Gain") ~ 1L,
    cn_category %in% c("Deletion", "Amplification") ~ 2L,
    TRUE ~ NA_integer_
  )
)

cnv_mapped <- cnv %>%
  transmute(
    SangerModelID = as.character(model_id),
    geneX         = as.character(symbol),
    cn_category   = as.character(cn_category),
    source        = as.character(source)
  ) %>%
  filter(
    !is.na(SangerModelID), SangerModelID != "",
    !is.na(geneX), geneX != "",
    !is.na(cn_category), cn_category != "",
    cn_category %in% CNV_LEVELS_KEEP
  ) %>%
  filter(geneX %in% protein_coding_genes) %>%
  distinct() %>%
  left_join(dist_tbl, by = "cn_category") %>%
  mutate(is_broad = if_else(source == "Broad", 1L, 0L)) %>%
  group_by(SangerModelID, geneX) %>%
  arrange(dist_neutral, desc(is_broad)) %>%
  slice(1) %>%
  ungroup() %>%
  select(SangerModelID, geneX, cn_category) %>%
  inner_join(id_map, by = "SangerModelID") %>%
  transmute(
    ModelID = ModelID,
    geneX = geneX,
    cn_category = factor(cn_category, levels = CNV_LEVELS_KEEP)
  ) %>%
  distinct() %>%
  mutate(
    cn_category = droplevels(cn_category),
    cn_category = relevel(cn_category, ref = "Neutral")
  )

# =========================================================
# 6) Feature selection at gene level before clustering
# =========================================================
feature_stats <- cnv_mapped %>%
  group_by(geneX) %>%
  summarise(
    n_models    = n_distinct(ModelID),
    n_levels    = n_distinct(cn_category),
    non_neutral = sum(cn_category != "Neutral"),
    min_level_n = min(table(droplevels(cn_category))),
    .groups = "drop"
  )

feature_list <- feature_stats %>%
  filter(
    n_models >= MIN_SAMPLES_FEATURE,
    n_levels >= 2,
    min_level_n >= MIN_GROUP_N
  ) %>%
  arrange(desc(non_neutral), desc(n_models)) %>%
  pull(geneX)

cnv_features <- cnv_mapped %>%
  filter(geneX %in% feature_list) %>%
  mutate(cn_category = droplevels(cn_category))

# =========================================================
# 7) Approximate CNV clusters via hierarchical clustering
# =========================================================
cnv_wide <- cnv_features %>%
  select(geneX, ModelID, cn_category) %>%
  mutate(
    ModelID = as.character(ModelID),
    cn_category = as.character(cn_category)
  ) %>%
  pivot_wider(names_from = ModelID, values_from = cn_category)

sig_cols <- setdiff(names(cnv_wide), "geneX")
sig_cols <- sort(sig_cols)

cnv_mat_df <- cnv_wide %>%
  select(all_of(sig_cols)) %>%
  mutate(across(everything(), as.factor))

gene_names <- cnv_wide$geneX

obs_mat <- !is.na(as.matrix(cnv_wide %>% select(all_of(sig_cols))))
overlap_mat <- obs_mat %*% t(obs_mat)

gower_dist <- cluster::daisy(cnv_mat_df, metric = "gower")
dmat <- as.matrix(gower_dist)
rownames(dmat) <- gene_names
colnames(dmat) <- gene_names

dmat[overlap_mat < MIN_CLUSTER_OVERLAP] <- 1
diag(dmat) <- 0

hc <- hclust(as.dist(dmat), method = HC_METHOD)
cluster_membership <- cutree(hc, h = CLUSTER_CUT_HEIGHT)

cluster_tbl <- tibble(
  geneX = gene_names,
  cluster_group = cluster_membership
) %>%
  arrange(cluster_group, geneX) %>%
  group_by(cluster_group) %>%
  mutate(cluster_id = sprintf("CLU%06d", cur_group_id())) %>%
  ungroup() %>%
  select(geneX, cluster_id)

cluster_map <- cluster_tbl %>%
  group_by(cluster_id) %>%
  summarise(
    n_genes   = n(),
    rep_gene  = sort(geneX)[1],
    gene_list = paste(sort(geneX), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genes), cluster_id)

rep_genes <- cluster_map %>% select(cluster_id, rep_gene)

cnv_clusters <- cnv_features %>%
  inner_join(rep_genes, by = c("geneX" = "rep_gene")) %>%
  transmute(
    ModelID = ModelID,
    geneX = cluster_id,
    cn_category = factor(as.character(cn_category), levels = CNV_LEVELS_KEEP)
  ) %>%
  distinct() %>%
  mutate(cn_category = droplevels(cn_category))

cluster_stats <- cnv_clusters %>%
  group_by(geneX) %>%
  summarise(
    n_models    = n_distinct(ModelID),
    n_levels    = n_distinct(cn_category),
    non_neutral = sum(cn_category != "Neutral"),
    min_level_n = min(table(droplevels(cn_category))),
    .groups = "drop"
  ) %>%
  left_join(
    cluster_map %>% transmute(geneX = cluster_id, n_genes, rep_gene, gene_list),
    by = "geneX"
  )

cluster_list <- cluster_stats %>%
  filter(
    n_models >= MIN_SAMPLES_FEATURE,
    n_levels >= 2,
    min_level_n >= MIN_GROUP_N
  ) %>%
  arrange(desc(non_neutral), desc(n_models)) %>%
  pull(geneX)

cnv_clusters_features <- cnv_clusters %>%
  filter(geneX %in% cluster_list) %>%
  mutate(cn_category = droplevels(cn_category))

# =========================================================
# 8) Feature-centric ANCOVA + contrasts (MULTI-LEVEL MAIN)
# =========================================================
ancova_one_feature <- function(gX, cnv_df, crispr_df, tissue_df, sex_df, culture_df, media_df, pca_df, clu_meta) {

  x <- cnv_df %>%
    filter(geneX == gX) %>%
    select(ModelID, cn_category)

  if (nrow(x) == 0) return(list(global = tibble(), contrasts = tibble()))

  df <- x %>%
    inner_join(crispr_df,  by = "ModelID") %>%
    inner_join(tissue_df,  by = "ModelID") %>%
    inner_join(sex_df,     by = "ModelID") %>%
    inner_join(culture_df, by = "ModelID") %>%
    inner_join(media_df,   by = "ModelID") %>%
    inner_join(pca_df,     by = "ModelID") %>%
    filter(
      is.finite(essentiality),
      !is.na(cn_category),
      !is.na(tissue)
    ) %>%
    mutate(
      cn_category = droplevels(cn_category),
      cn_category = relevel(cn_category, ref = "Neutral"),
      tissue      = droplevels(tissue),
      sex         = droplevels(sex),
      culture     = droplevels(culture),
      media_base  = droplevels(media_base)
    )

  if (nrow(df) == 0) return(list(global = tibble(), contrasts = tibble()))

  counts_wide <- df %>%
    group_by(geneY, cn_category) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = cn_category, values_from = n, values_fill = 0)

  level_cols <- setdiff(names(counts_wide), "geneY")
  counts_mat <- as.matrix(counts_wide %>% select(all_of(level_cols)))

  n_total  <- rowSums(counts_mat)
  n_groups <- rowSums(counts_mat > 0)

  min_group <- apply(counts_mat, 1, function(v) {
    v_pos <- v[v > 0]
    if (length(v_pos) == 0) 0 else min(v_pos)
  })

  neutral_n <- if ("Neutral" %in% level_cols) counts_wide$Neutral else rep(0, nrow(counts_wide))
  altered_n <- n_total - neutral_n

  ok_genes <- counts_wide %>%
    mutate(
      n_total   = n_total,
      n_groups  = n_groups,
      min_group = min_group,
      neutral_n = neutral_n,
      altered_n = altered_n
    ) %>%
    filter(
      n_total >= MIN_Y_N,
      n_groups >= 2,
      min_group >= MIN_GROUP_N,
      neutral_n >= MIN_GROUP_N,
      altered_n >= MIN_GROUP_N
    ) %>%
    pull(geneY)

  if (length(ok_genes) == 0) return(list(global = tibble(), contrasts = tibble()))

  df_ok <- df %>% filter(geneY %in% ok_genes)

  ancova_formula <- as.formula(
    paste(
      "essentiality ~ tissue + sex + culture + media_base +",
      paste(PCA_PC_NAMES, collapse = " + "),
      "+ cn_category"
    )
  )

  fit_tbl <- df_ok %>%
    group_by(geneY) %>%
    summarise(
      n_samples = n(),
      fit = list(tryCatch(
        aov(ancova_formula, data = pick(everything())),
        error = function(e) NULL
      )),
      tidy = list(tryCatch(
        broom::tidy(fit[[1]]),
        error = function(e) NULL
      )),
      .groups = "drop"
    )

  global_tbl <- fit_tbl %>%
    unnest(tidy) %>%
    group_by(geneY, n_samples, fit) %>%
    summarise(
      p_value    = suppressWarnings(min(p.value[term == "cn_category"], na.rm = TRUE)),
      ss_feat    = suppressWarnings(sum(sumsq[term == "cn_category"], na.rm = TRUE)),
      ss_resid   = suppressWarnings(sum(sumsq[term == "Residuals"], na.rm = TRUE)),
      ss_total   = suppressWarnings(sum(sumsq, na.rm = TRUE)),
      df_between = suppressWarnings(sum(df[term == "cn_category"], na.rm = TRUE)),
      df_within  = suppressWarnings(sum(df[term == "Residuals"], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      effect_size = ss_feat / (ss_feat + ss_resid),
      geneX = gX,
      feature = paste(
        "tissue + sex + culture + media_base +",
        paste(PCA_PC_NAMES, collapse = " + "),
        "+ cn_category"
      )
    )

  emm_tbl <- fit_tbl %>%
    mutate(
      contrast = map(fit, function(f) {
        if (is.null(f)) return(tibble())

        emm <- emmeans::emmeans(f, ~ cn_category)

        emmeans::contrast(
          emm,
          method = "trt.vs.ctrl",
          ref = "Neutral"
        ) %>%
          as.data.frame() %>%
          transmute(
            contrast = as.character(contrast),
            level = str_trim(str_remove(as.character(contrast), " - Neutral$")),
            estimate_adj = as.numeric(estimate),
            contrast_p_value = as.numeric(p.value)
          )
      })
    ) %>%
    select(geneY, contrast) %>%
    unnest(contrast)

  neu_stats <- df_ok %>%
    filter(cn_category == "Neutral") %>%
    group_by(geneY) %>%
    summarise(
      n_baseline  = n(),
      mu_baseline = mean(essentiality, na.rm = TRUE),
      .groups = "drop"
    )

  lvl_stats <- df_ok %>%
    filter(cn_category != "Neutral") %>%
    group_by(geneY, cn_category) %>%
    summarise(
      n_level  = n(),
      mu_level = mean(essentiality, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      geneY = geneY,
      level = as.character(cn_category),
      n_level = n_level,
      mu_level = mu_level,
      contrast = paste0(level, " - Neutral")
    )

  contrasts_tbl <- lvl_stats %>%
    inner_join(neu_stats, by = "geneY") %>%
    left_join(emm_tbl, by = c("geneY", "contrast", "level")) %>%
    transmute(
      geneX = gX,
      geneY = geneY,
      contrast = contrast,
      baseline = "Neutral",
      level = level,
      n_baseline = as.integer(n_baseline),
      n_level = as.integer(n_level),
      mu_baseline = as.numeric(mu_baseline),
      mu_level = as.numeric(mu_level),
      delta_mean = as.numeric(estimate_adj),
      contrast_p_value = as.numeric(contrast_p_value)
    ) %>%
    filter(is.finite(delta_mean))

  best_contrast_tbl <- contrasts_tbl %>%
    group_by(geneY) %>%
    arrange(desc(abs(delta_mean)), contrast_p_value, contrast, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(geneY, delta_mean, contrast, level, contrast_p_value)

  global_tbl <- global_tbl %>%
    left_join(best_contrast_tbl, by = "geneY") %>%
    mutate(
      signed_eta2 = sign(delta_mean) * effect_size
    ) %>%
    filter(
      is.finite(p_value),
      is.finite(effect_size),
      is.finite(delta_mean),
      is.finite(contrast_p_value),
      is.finite(signed_eta2),
      ss_total > 0,
      ss_resid > 0,
      df_between > 0,
      df_within > 0
    ) %>%
    select(
      geneX, geneY, feature,
      n_samples,
      p_value,
      effect_size,
      df_between,
      df_within,
      ss_total,
      ss_feat,
      ss_resid,
      signed_eta2,
      delta_mean,
      contrast,
      level,
      contrast_p_value
    )

  meta_one <- clu_meta %>% filter(cluster_id == gX) %>% slice(1)
  if (nrow(meta_one) == 1) {
    global_tbl <- global_tbl %>%
      mutate(
        rep_gene = meta_one$rep_gene,
        n_genes  = meta_one$n_genes
      )

    contrasts_tbl <- contrasts_tbl %>%
      mutate(
        rep_gene = meta_one$rep_gene,
        n_genes  = meta_one$n_genes
      )
  }

  list(global = global_tbl, contrasts = contrasts_tbl)
}

# =========================================================
# 9) Iterate clusters
# =========================================================
res_list <- mclapply(
  cluster_list,
  ancova_one_feature,
  cnv_df = cnv_clusters_features,
  crispr_df = crispr_long,
  tissue_df = tissue_map,
  sex_df = sex_map,
  culture_df = culture_map,
  media_df = media_map,
  pca_df = pca_map,
  clu_meta = cluster_map,
  mc.cores = N_CORES
)

ancova_results   <- map_dfr(res_list, "global")
contrast_results <- map_dfr(res_list, "contrasts")

ancova_results <- ancova_results %>%
  group_by(geneX) %>%
  mutate(FDR_by_geneX = p.adjust(p_value, "fdr")) %>%
  ungroup() %>%
  mutate(FDR_global = p.adjust(p_value, "fdr"))

ancova_results <- ancova_results %>%
  mutate(
    neglog10_p = -log10(p_value),
    is_hit = is.finite(FDR_by_geneX) & FDR_by_geneX <= 0.05
  ) %>%
  left_join(
    cluster_map %>% transmute(geneX = cluster_id, gene_list),
    by = "geneX"
  )

# =========================================================
# 10) Report
# =========================================================
n_geneX_input     <- n_distinct(cnv_mapped$geneX)
n_geneX_feature   <- length(feature_list)
n_clusters_total  <- nrow(cluster_map)
n_clusters_select <- length(cluster_list)

cluster_size_bins <- cluster_map %>%
  mutate(
    size_bin = case_when(
      n_genes == 1 ~ "1",
      n_genes >= 2 & n_genes <= 5 ~ "2-5",
      n_genes >= 6 & n_genes <= 10 ~ "6-10",
      n_genes > 10 ~ ">10",
      TRUE ~ "other"
    )
  ) %>%
  count(size_bin, sort = TRUE)

top_clusters <- ancova_results %>%
  group_by(geneX) %>%
  arrange(desc(abs(delta_mean)), contrast_p_value, p_value, .by_group = TRUE) %>%
  summarise(
    n_tests        = n(),
    n_hits         = sum(is_hit, na.rm = TRUE),
    best_FDR       = suppressWarnings(min(FDR_by_geneX, na.rm = TRUE)),
    best_p         = suppressWarnings(min(p_value, na.rm = TRUE)),
    top_geneY      = first(geneY),
    top_contrast   = first(contrast),
    top_level      = first(level),
    top_delta_mean = first(delta_mean),
    top_contrast_p = first(contrast_p_value),
    max_abs_delta  = suppressWarnings(max(abs(delta_mean), na.rm = TRUE)),
    rep_gene       = dplyr::first(na.omit(rep_gene)),
    n_genes        = dplyr::first(na.omit(n_genes)),
    gene_list      = dplyr::first(na.omit(gene_list)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_hits), best_FDR, desc(abs(top_delta_mean)), top_contrast_p, best_p) %>%
  slice_head(n = 50)

report_lines <- c(
  "CNV CLUSTERS MULTI-LEVEL REPORT (MAIN BRANCH)",
  "============================================",
  "",
  paste0("Date: ", Sys.time()),
  "",
  "Model:",
  paste0(" - essentiality ~ tissue + sex + culture + media_base + ",
         paste(PCA_PC_NAMES, collapse = " + "), " + cn_category"),
  "",
  "Role of this branch:",
  " - main ranking / top-hits branch",
  " - gatekeeper: FDR_by_geneX on cn_category",
  " - main effect interpretation: best adjusted contrast vs Neutral",
  "",
  paste0("GeneY skewness threshold: ", MIN_Y_SKEW),
  paste0("Sex-linked gene filter enabled: ", REMOVE_SEX_LINKED_GENES),
  "",
  "Best-contrast rule:",
  " - choose largest |delta_mean| vs Neutral",
  " - use smallest contrast p-value as tie-break",
  "",
  paste0("Cluster method: hierarchical clustering (", HC_METHOD, " linkage)"),
  paste0("Distance: 1 - concordance (Gower on nominal CNV categories)"),
  paste0("Cluster match threshold: ", CLUSTER_MATCH_THRESHOLD),
  paste0("Cut height: ", CLUSTER_CUT_HEIGHT),
  paste0("Minimum overlap for valid gene-gene comparison: ", MIN_CLUSTER_OVERLAP),
  "",
  "Gene filters:",
  " - geneX: protein-coding only",
  " - geneY: protein-coding only",
  " - geneY common essentials removed only when annotation == 'common essential'",
  if (REMOVE_SEX_LINKED_GENES) " - sex-linked genes removed using Gene.csv location" else " - sex-linked genes retained",
  " - NA annotation not used as exclusion criterion",
  "",
  "Counts:",
  paste0(" - total protein-coding geneX in cnv_mapped: ", n_geneX_input),
  paste0(" - geneX passing feature prefilter: ", n_geneX_feature),
  paste0(" - total clusters built: ", n_clusters_total),
  paste0(" - clusters passing multi-level selection: ", n_clusters_select),
  paste0(" - ANCOVA rows out: ", nrow(ancova_results)),
  paste0(" - contrast rows out: ", nrow(contrast_results)),
  "",
  "Cluster size distribution:"
)

if (nrow(cluster_size_bins) > 0) {
  report_lines <- c(
    report_lines,
    paste0(" - ", cluster_size_bins$size_bin, ": ", cluster_size_bins$n)
  )
}

report_lines <- c(
  report_lines,
  "",
  "Top 50 ranked clusters (main multi-level ranking):",
  paste(capture.output(print(top_clusters, n = 50, width = Inf)), collapse = "\n"),
  "",
  "Files written:",
  paste0(" - ", OUT_GLOBAL),
  paste0(" - ", OUT_CONTR),
  paste0(" - ", OUT_SELECTED),
  paste0(" - ", OUT_TOP50),
  paste0(" - ", OUT_MAP),
  paste0(" - ", OUT_REPORT),
  ""
)

writeLines(report_lines, OUT_REPORT)

# =========================================================
# 11) Output
# =========================================================
write_csv(ancova_results,   OUT_GLOBAL)
write_csv(contrast_results, OUT_CONTR)
write_csv(cluster_stats %>% filter(geneX %in% cluster_list), OUT_SELECTED)
write_csv(cluster_map, OUT_MAP)
write_csv(top_clusters, OUT_TOP50)

message("OK: wrote multi-level ANCOVA global + all contrasts vs Neutral.")
message("OK: multi-level branch is the main ranking / top-hits branch.")