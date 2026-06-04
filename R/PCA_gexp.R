#!/usr/bin/env Rscript

BASE_DIR <- "~/project_inescid/"
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "output/csv")
PLOT_DIR <- file.path(BASE_DIR, "output/plots_pca_gexp")
REP_DIR  <- file.path(BASE_DIR, "output/reports_pca_gexp")

dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REP_DIR,  recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(stringr)

# =========================================================
# 0) Parâmetros
# =========================================================
EXPR_FILE <- file.path(DATA_DIR, "gexp_clines_voom.csv")
MODEL_FILE <- file.path(DATA_DIR, "Model.csv")

N_PCS <- 10
GENE_ID_COL <- 1       
CENTER_DATA <- TRUE
SCALE_DATA  <- FALSE    

# =========================================================
# 1) Ler dados
# =========================================================
expr_raw <- read_csv(EXPR_FILE, show_col_types = FALSE)
model    <- read_csv(MODEL_FILE, show_col_types = FALSE)

if (ncol(expr_raw) < 3) {
  stop("A matriz de expressão parece ter poucas colunas. Esperava genes nas linhas e SIDM nas colunas.")
}

# =========================================================
# 2) Separar gene IDs e matriz de expressão
# =========================================================
gene_ids <- expr_raw[[GENE_ID_COL]]

expr_mat_df <- expr_raw[, -GENE_ID_COL, drop = FALSE]

sidm_cols <- names(expr_mat_df)
sidm_cols <- sidm_cols[str_detect(sidm_cols, "^SIDM")]

if (length(sidm_cols) == 0) {
  stop("Não encontrei colunas SIDM na matriz de expressão.")
}

expr_mat_df <- expr_mat_df %>%
  select(all_of(sidm_cols))

expr_mat <- as.matrix(expr_mat_df)
mode(expr_mat) <- "numeric"

rownames(expr_mat) <- make.unique(as.character(gene_ids))

# =========================================================
# 4) Transpor: PCA em amostras × genes
# =========================================================
X <- t(expr_mat)

keep_samples <- rowSums(is.finite(X)) > 0
X <- X[keep_samples, , drop = FALSE]

if (nrow(X) < 3) {
  stop("Após filtragem, sobraram poucas amostras para PCA.")
}

col_means <- colMeans(X, na.rm = TRUE)
na_idx <- which(!is.finite(X), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  X[na_idx] <- col_means[na_idx[, 2]]
}

col_var <- apply(X, 2, var)
keep_cols <- is.finite(col_var) & col_var > 0
X <- X[, keep_cols, drop = FALSE]

if (ncol(X) < 2) {
  stop("Após filtragem, sobraram poucas features (genes) para PCA.")
}

# =========================================================
# 5) PCA
# =========================================================
pca_fit <- prcomp(
  X,
  center = CENTER_DATA,
  scale. = SCALE_DATA
)

# =========================================================
# 6) Variância explicada
# =========================================================
var_explained <- (pca_fit$sdev ^ 2) / sum(pca_fit$sdev ^ 2)

var_tbl <- tibble(
  PC = paste0("PC", seq_along(var_explained)),
  variance_explained = var_explained,
  cumulative_variance = cumsum(var_explained)
)

# =========================================================
# 7) Scores dos PCs por amostra
# =========================================================
scores_tbl <- as_tibble(pca_fit$x[, seq_len(min(N_PCS, ncol(pca_fit$x))), drop = FALSE]) %>%
  mutate(SangerModelID = rownames(pca_fit$x)) %>%
  relocate(SangerModelID)

id_map <- model %>%
  select(ModelID, SangerModelID) %>%
  filter(!is.na(ModelID), !is.na(SangerModelID), SangerModelID != "") %>%
  distinct() %>%
  group_by(SangerModelID) %>%
  slice(1) %>%
  ungroup()

scores_tbl <- scores_tbl %>%
  left_join(id_map, by = "SangerModelID") %>%
  relocate(ModelID, .after = SangerModelID)

# =========================================================
# 8) Loadings dos top PCs
# =========================================================
loadings_tbl <- as_tibble(
  pca_fit$rotation[, seq_len(min(N_PCS, ncol(pca_fit$rotation))), drop = FALSE],
  rownames = "gene"
)

# =========================================================
# 9) Scree plot
# =========================================================
plot_n <- min(20, nrow(var_tbl))

p_scree <- var_tbl %>%
  slice_head(n = plot_n) %>%
  mutate(PC = factor(PC, levels = PC)) %>%
  ggplot(aes(x = PC, y = variance_explained)) +
  geom_col() +
  geom_point(aes(y = cumulative_variance), size = 1.5) +
  geom_line(aes(y = cumulative_variance, group = 1), linewidth = 0.5) +
  labs(
    title = "Gene expression PCA scree plot",
    x = "Principal Component",
    y = "Variance explained"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(PLOT_DIR, "scree_plot_gexp_pca.png"),
  p_scree,
  width = 8,
  height = 5.5,
  dpi = 300
)

# =========================================================
# 10) Scatter plots dos PCs
# =========================================================
make_pc_scatter <- function(df, pcx, pcy, var_tbl, out_file) {
  vx <- var_tbl$variance_explained[var_tbl$PC == pcx]
  vy <- var_tbl$variance_explained[var_tbl$PC == pcy]

  p <- ggplot(df, aes_string(x = pcx, y = pcy)) +
    geom_point(alpha = 0.7, size = 1.4) +
    labs(
      title = paste0("PCA of gene expression: ", pcx, " vs ", pcy),
      x = paste0(pcx, " (", round(100 * vx, 2), "%)"),
      y = paste0(pcy, " (", round(100 * vy, 2), "%)")
    ) +
    theme_classic(base_size = 12)

  ggsave(
    out_file,
    p,
    width = 6.5,
    height = 5.5,
    dpi = 300
  )
}

pcs_available <- names(scores_tbl)[str_detect(names(scores_tbl), "^PC[0-9]+$")]

plot_pairs <- list(
  c("PC1", "PC2"),
  c("PC1", "PC3"),
  c("PC2", "PC3"),
  c("PC1", "PC4"),
  c("PC1", "PC5")
)

for (pair in plot_pairs) {
  if (all(pair %in% pcs_available)) {
    make_pc_scatter(
      df = scores_tbl,
      pcx = pair[1],
      pcy = pair[2],
      var_tbl = var_tbl,
      out_file = file.path(PLOT_DIR, paste0("scatter_", pair[1], "_vs_", pair[2], ".png"))
    )
  }
}

# =========================================================
# 11) Report TXT
# =========================================================
n_samples <- nrow(X)
n_genes   <- ncol(X)

top10_var <- var_tbl %>%
  slice_head(n = min(10, nrow(var_tbl)))

lines <- c(
  "PCA GENE EXPRESSION REPORT",
  "==========================",
  "",
  paste0("Input file: ", EXPR_FILE),
  paste0("Model file: ", MODEL_FILE),
  "",
  "Settings:",
  paste0(" - N_PCS: ", N_PCS),
  paste0(" - center: ", CENTER_DATA),
  paste0(" - scale: ", SCALE_DATA),
  "",
  "Dimensions used in PCA:",
  paste0(" - n_samples: ", n_samples),
  paste0(" - n_genes: ", n_genes),
  "",
  "Top variance explained:",
  paste(capture.output(print(top10_var, n = 10)), collapse = "\n"),
  "",
  "Files written:",
  paste0(" - ", file.path(OUT_DIR, "gexp_PCA_scores_top10.csv")),
  paste0(" - ", file.path(OUT_DIR, "gexp_PCA_variance_explained.csv")),
  paste0(" - ", file.path(OUT_DIR, "gexp_PCA_loadings_top10.csv")),
  paste0(" - ", file.path(PLOT_DIR, "scree_plot_gexp_pca.png")),
  paste0(" - scatter plots in: ", PLOT_DIR),
  ""
)

writeLines(lines, file.path(REP_DIR, "PCA_gexp_report.txt"))

# =========================================================
# 12) Escrever outputs
# =========================================================
write_csv(
  scores_tbl,
  file.path(OUT_DIR, "gexp_PCA_scores_top10.csv")
)

write_csv(
  var_tbl,
  file.path(OUT_DIR, "gexp_PCA_variance_explained.csv")
)

write_csv(
  loadings_tbl,
  file.path(OUT_DIR, "gexp_PCA_loadings_top10.csv")
)
