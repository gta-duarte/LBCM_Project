#!/usr/bin/env Rscript

# =========================================================
# CNV_cluster_driver_prioritization_FULL_ALIGNED.R
# =========================================================

BASE_DIR <- "~/project_inescid/"
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "output/csv/")
REP_DIR  <- file.path(BASE_DIR, "output/reports/")
TOP50_PATH <- file.path(OUT_DIR, "CNV_top50_ranked_clusters_MULTILEVEL.csv")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(broom)

PCA_PC_NAMES <- paste0("PC", 1:5)

# ANCOVA thresholds
MIN_Y_N     <- 300
MIN_GROUP_N <- 9

# =========================================================
# 1) LOAD
# =========================================================
top50_tbl   <- read_csv(file.path(OUT_DIR, "CNV_top50_ranked_clusters_MULTILEVEL.csv"))
cluster_map <- read_csv(file.path(OUT_DIR, "CNV_clusters_map_MULTILEVEL.csv"))

crispr <- read_csv(file.path(DATA_DIR, "CRISPRGeneEffect.csv"))
gexp   <- read_csv(file.path(DATA_DIR, "gexp_clines_voom.csv"))
cnv    <- read_csv(file.path(DATA_DIR, "WES_pureCN_CNV_genes_20250207.csv"))
model  <- read_csv(file.path(DATA_DIR, "Model.csv"))
pca    <- read_csv(file.path(DATA_DIR, "gexp_PCA_scores_top10.csv"))


# =========================================================
# 2) TOP geneY por cluster
# =========================================================
top_pairs <- top50_tbl %>%
  transmute(
    cluster_id = as.character(geneX),
    geneY      = as.character(top_geneY)
  ) %>%
  distinct()

cluster_long <- cluster_map %>%
  separate_rows(gene_list, sep = ";") %>%
  transmute(
    cluster_id = as.character(cluster_id),
    geneX      = as.character(gene_list)
  )

pairs <- cluster_long %>%
  inner_join(top_pairs, by = "cluster_id")

# =========================================================
# 3) CRISPR
# =========================================================
crispr_long <- crispr %>%
  rename(ModelID = `...1`) %>%
  pivot_longer(-ModelID, names_to = "gene_raw", values_to = "essentiality") %>%
  mutate(
    geneY = str_extract(gene_raw, "^[^ ]+"),
    essentiality = suppressWarnings(as.numeric(essentiality))
  ) %>%
  select(ModelID, geneY, essentiality)

# =========================================================
# 4) GEXP
# =========================================================
id_map <- model %>%
  select(ModelID, SangerModelID) %>%
  filter(!is.na(ModelID), !is.na(SangerModelID), SangerModelID != "") %>%
  distinct() %>%
  group_by(SangerModelID) %>%
  slice(1) %>%
  ungroup()

gexp_long <- gexp %>%
  rename(geneX = `...1`) %>%
  pivot_longer(
    cols = -geneX,
    names_to = "SangerModelID",
    values_to = "expression"
  ) %>%
  mutate(
    geneX = as.character(geneX),
    SangerModelID = as.character(SangerModelID),
    expression = as.numeric(expression)
  ) %>%
  filter(!is.na(geneX), geneX != "", !is.na(SangerModelID), SangerModelID != "") %>%
  inner_join(id_map, by = "SangerModelID") %>%
  transmute(
    ModelID = ModelID,
    geneX = geneX,
    expression = expression
  ) %>%
  distinct()

# =========================================================
# 5) CNV 
# =========================================================
dist_tbl <- tibble(
  cn_category = c("Loss","Neutral","Gain","Amplification","Deletion"),
  dist_neutral = c(1L,0L,1L,2L,2L)
)

cnv_long <- cnv %>%
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
    cn_category %in% c("Loss","Neutral","Gain","Amplification","Deletion")
  ) %>%
  distinct() %>%
  left_join(dist_tbl, by = "cn_category") %>%
  mutate(is_broad = if_else(source == "Broad", 1L, 0L)) %>%
  group_by(SangerModelID, geneX) %>%
  arrange(dist_neutral, desc(is_broad)) %>%
  slice(1) %>%
  ungroup() %>%
  inner_join(id_map, by = "SangerModelID") %>%
  transmute(
    ModelID = ModelID,
    geneX = geneX,
    cn_category = factor(cn_category,
      levels = c("Loss","Neutral","Gain","Amplification","Deletion"))
  ) %>%
  distinct() %>%
  mutate(cn_category = droplevels(cn_category))

# =========================================================
# 6) COVARIATES
# =========================================================
tissue_map <- model %>%
  select(ModelID, OncotreeLineage) %>%
  rename(tissue = OncotreeLineage) %>%
  filter(!is.na(tissue), tissue != "") %>%
  mutate(tissue = factor(tissue)) %>%
  distinct()

sex_map <- model %>%
  select(ModelID, Sex) %>%
  rename(sex_raw = Sex) %>%
  mutate(
    sex = case_when(
      is.na(sex_raw) | sex_raw == "" ~ "Unknown",
      sex_raw %in% c("Male","Female") ~ sex_raw,
      TRUE ~ "Unknown"
    ),
    sex = factor(sex, levels = c("Male","Female","Unknown"))
  ) %>%
  select(ModelID, sex) %>%
  distinct()

culture_map <- model %>%
  select(ModelID, GrowthPattern) %>%
  rename(culture_raw = GrowthPattern) %>%
  mutate(
    culture = case_when(
      is.na(culture_raw) | culture_raw == "" ~ "unknown",
      str_to_lower(culture_raw) %in% "adherent" ~ "adherent",
      str_to_lower(culture_raw) %in% "suspension" ~ "suspension",
      TRUE ~ "unknown"
    ),
    culture = factor(culture, levels = c("adherent","suspension","unknown"))
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
  mutate(media_base = factor(media_base,
    levels = c("RPMI","DMEM","MEM","IMDM","F12","OTHER","unknown"))) %>%
  distinct()

pca_map <- pca %>%
  select(ModelID, all_of(PCA_PC_NAMES)) %>%
  distinct()

cov_df <- tissue_map %>%
  inner_join(sex_map, by = "ModelID") %>%
  inner_join(culture_map, by = "ModelID") %>%
  inner_join(media_map, by = "ModelID") %>%
  inner_join(pca_map, by = "ModelID")

# =========================================================
# 7) CORE FUNCTION 
# =========================================================
score_one_pair <- function(gX, gY, cluster_id) {

  df <- cnv_long %>%
  filter(geneX == gX) %>%
  inner_join(gexp_long %>% filter(geneX == gX), by = c("ModelID","geneX")) %>%
  inner_join(crispr_long %>% filter(geneY == gY), by = "ModelID") %>%
  inner_join(cov_df, by = "ModelID") %>%
  filter(
    is.finite(expression),
    is.finite(essentiality),
    !is.na(cn_category),
    !is.na(tissue)
  ) %>%
  mutate(cn_category = relevel(cn_category, ref = "Neutral"))

  if (nrow(df) == 0) return(NULL)

  # =========================
  # ANCOVA GUARDRAILS
  # =========================
  counts <- df %>%
    count(cn_category)

  n_total  <- sum(counts$n)
  n_groups <- sum(counts$n > 0)
  min_group <- min(counts$n)
  neutral_n <- sum(counts$n[counts$cn_category == "Neutral"])
  altered_n <- n_total - neutral_n

  if (
    n_total < MIN_Y_N ||
    n_groups < 2 ||
    min_group < MIN_GROUP_N ||
    neutral_n < MIN_GROUP_N ||
    altered_n < MIN_GROUP_N
  ) return(NULL)

  # =========================
  # MODEL 1
  # =========================
  fit1 <- lm(
    as.formula(paste(
      "expression ~ cn_category + tissue + sex + culture + media_base +",
      paste(PCA_PC_NAMES, collapse = " + ")
    )),
    data = df
  )

  t1 <- broom::tidy(fit1) %>%
    filter(str_detect(term, "cn_category"))

  t1_best <- t1 %>%
  arrange(desc(abs(estimate)), p.value) %>%
  slice(1)

  beta1 <- abs(t1_best$estimate)
  p1    <- t1_best$p.value

  # =========================
  # MODEL 2
  # =========================
  fit2 <- lm(
    as.formula(paste(
      "essentiality ~ expression + tissue + sex + culture + media_base +",
      paste(PCA_PC_NAMES, collapse = " + ")
    )),
    data = df
  )

  t2 <- broom::tidy(fit2) %>%
    filter(term == "expression")

  beta2 <- abs(t2$estimate)
  p2    <- t2$p.value

  tibble(
    cluster_id = cluster_id,
    geneX = gX,
    geneY = gY,
    beta1 = beta1,
    beta2 = beta2,
    p1 = p1,
    p2 = p2,
    p_combined = max(p1, p2),
    score = beta1 * beta2
  )
}

# =========================================================
# 8) RUN
# =========================================================
results <- pairs %>%
  mutate(res = pmap(list(geneX, geneY, cluster_id), score_one_pair)) %>%
  select(res) %>%
  unnest(res)

# =========================================================
# 9) FDR + ranking
# =========================================================
results <- results %>%
  group_by(cluster_id) %>%
  mutate(FDR = p.adjust(p_combined, "fdr")) %>%
  arrange(desc(score), FDR, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# =========================================================
# 10) SAVE
# =========================================================
write_csv(results, file.path(OUT_DIR, "CNV_cluster_driver_scores_FINAL.csv"))

# =========================================================
# 11) REPORT
# =========================================================
OUT_REPORT <- file.path(REP_DIR, "CNV_cluster_driver_prioritization_report.txt")

report_tbl <- results %>%
  group_by(cluster_id) %>%
  group_modify(~ {
    sig <- .x %>% filter(FDR < 0.05)

    if (nrow(sig) > 0) {
      sig %>% arrange(FDR, desc(score), p_combined)
    } else {
      .x %>% arrange(FDR, desc(score), p_combined) %>% slice(1)
    }
  }) %>%
  ungroup() %>%
  arrange(cluster_id, FDR, desc(score), p_combined)

report_lines <- c(
  "CNV CLUSTER DRIVER PRIORITIZATION REPORT",
  "========================================",
  paste0("Date: ", Sys.time()),
  "",
  paste0("Input top-50 ranking file: ", basename(TOP50_PATH)),
  paste0("Clusters analysed: ", n_distinct(results$cluster_id)),
  paste0("Total geneX tested: ", n_distinct(results$geneX)),
  paste0("Total tested pairs: ", nrow(results)),
  "",
  "Model 1:",
  paste0(
    " - expression ~ cn_category + tissue + sex + culture + media_base + ",
    paste(PCA_PC_NAMES, collapse = " + ")
  ),
  " - beta1 = absolute strongest cn_category coefficient vs Neutral",
  " - p1 = p-value of that strongest coefficient",
  "",
  "Model 2:",
  paste0(
    " - essentiality ~ expression + tissue + sex + culture + media_base + ",
    paste(PCA_PC_NAMES, collapse = " + ")
  ),
  " - beta2 = absolute coefficient for expression",
  " - p2 = p-value for expression",
  "",
  "Score definition:",
  " - score = beta1 * beta2",
  " - p_combined = max(p1, p2)",
  " - FDR computed by cluster on p_combined",
  "",
  "Reported results:",
  " - all FDR < 0.05 results within each cluster",
  " - if a cluster has no FDR < 0.05 result, report the best one only",
  "",
  "Results:",
  paste(capture.output(print(report_tbl, n = Inf, width = Inf)), collapse = "\n")
)

writeLines(report_lines, OUT_REPORT)

message("DONE: fully aligned + ANCOVA-consistent.")