#!/usr/bin/env Rscript
# 07_scatter_plot.R
# Scatter plot comparando a direção do efeito (beta MVP × beta exoma)
# para as coding variants genotipadas no exoma.
#
# Input : resultados/tabelas/lookup_completo.csv  (output do 02_lookup_exoma.R)
# Output: resultados/figuras/scatter_direcao.png

# ─────────────────────────── CONFIG ───────────────────────────
LOOKUP_CSV <- "resultados/tabelas/lookup_completo.csv"
OUT_FIGURAS <- "resultados/figuras"
OUT_PDF <- "resultados/figuras/scatter_direcao.pdf"   # opcional (também salva PDF)

COL_BETA_MVP <- "Beta Population"
COL_PIP      <- "Overall PIP"
# ───────────────────────────────────────────────────────────────

# Caminhos relativos à raiz do repo (o pai do dir deste script), para o script
# funcionar rodado de qualquer CWD (ex.: via PBS em scripts/).
args0 <- commandArgs(trailingOnly = FALSE)
fidx <- grep("^--file=", args0)
SCRIPT_DIR <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
LOOKUP_CSV  <- file.path(REPO_ROOT, LOOKUP_CSV)
OUT_FIGURAS <- file.path(REPO_ROOT, OUT_FIGURAS)
OUT_PDF     <- file.path(REPO_ROOT, OUT_PDF)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
if (!dir.exists(OUT_FIGURAS)) dir.create(OUT_FIGURAS, recursive = TRUE)

lookup <- read.csv(LOOKUP_CSV, stringsAsFactors = FALSE, check.names = FALSE)
cat("== 07_scatter_plot.R ==\n")

plot_data <- lookup %>%
  filter(classificacao != "Nao_genotipada",
         !is.na(.data[[COL_BETA_MVP]]), !is.na(BETA_exoma))

cat(sprintf("Pontos plotados: %d\n", nrow(plot_data)))

if (nrow(plot_data) == 0) {
  stop("Nenhum ponto genotipado com beta disponível. Rode 02_lookup_exoma.R antes.")
}

p <- ggplot(plot_data, aes(x = .data[[COL_BETA_MVP]], y = BETA_exoma,
                            color = classificacao,
                            size = suppressWarnings(as.numeric(.data[[COL_PIP]])))) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +
  scale_color_manual(values = c("Direcao_consistente" = "#1b9e77",
                                 "Discordante" = "#d95f02",
                                 "Indeterminada" = "gray50")) +
  labs(title = "Comparação de direção do efeito: MVP × Exoma",
       x = "Beta (MVP)", y = "Beta (Exoma)",
       color = "Classificação", size = "PIP") +
  theme_minimal() +
  guides(size = guide_legend(override.aes = list(alpha = 1)))

ggsave(file.path(OUT_FIGURAS, "scatter_direcao.png"), p, width = 8, height = 6, dpi = 300)
if (requireNamespace("grDevices", quietly = TRUE)) {
  pdf(OUT_PDF, width = 8, height = 6)
  print(p)
  dev.off()
}

cat("OK. Scatter salvo em ", OUT_FIGURAS, "/scatter_direcao.png\n", sep = "")
