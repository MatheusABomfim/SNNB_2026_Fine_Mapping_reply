#!/usr/bin/env Rscript
# 03_analise_ancestralidade.R
# Análise por população (AFR / AMR / EUR / EAS) da consistência da direção
# do efeito entre MVP e exoma.
#
# Input : resultados/tabelas/lookup_completo.csv  (output do 02_lookup_exoma.R)
# Output: resultados/tabelas/sumario_por_populacao.csv,
#         resultados/tabelas/variantes_multi_pop.csv,
#         resultados/tabelas/matriz_z_por_populacao.csv

# ─────────────────────────── CONFIG ───────────────────────────
LOOKUP_CSV <- "resultados/tabelas/lookup_completo.csv"
OUT_TABELAS <- "resultados/tabelas"

# Colunas (conforme lookup_completo.csv gerado pelo 02)
COL_POP        <- "Population"
COL_RSID       <- "RSID"
COL_CLASS      <- "classificacao"
COL_Z_EXOMA    <- "Z_exoma"
# ───────────────────────────────────────────────────────────────

# Caminhos relativos à raiz do repo (o pai do dir deste script), para o script
# funcionar rodado de qualquer CWD (ex.: via PBS em scripts/).
args0 <- commandArgs(trailingOnly = FALSE)
fidx <- grep("^--file=", args0)
SCRIPT_DIR <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
LOOKUP_CSV   <- file.path(REPO_ROOT, LOOKUP_CSV)
OUT_TABELAS  <- file.path(REPO_ROOT, OUT_TABELAS)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
if (!dir.exists(OUT_TABELAS)) dir.create(OUT_TABELAS, recursive = TRUE)

lookup <- read.csv(LOOKUP_CSV, stringsAsFactors = FALSE, check.names = FALSE)
cat("== 03_analise_ancestralidade.R ==\n")

# 1. Sumário por população (apenas genotipadas)
genotipadas <- lookup %>% filter(classificacao != "Nao_genotipada")

pop_summary <- genotipadas %>%
  group_by(.data[[COL_POP]]) %>%
  summarise(
    n = n(),
    n_consistentes = sum(classificacao %in%
                           c("Replicada_p0.05", "Direcao_consistente")),
    n_replicadas = sum(classificacao == "Replicada_p0.05"),
    n_discordantes = sum(classificacao == "Discordante"),
    taxa_consistencia = n_consistentes / n,
    .groups = "drop"
  )

cat("\n--- Sumário por população ---\n")
print(pop_summary)
write.csv(pop_summary, file.path(OUT_TABELAS, "sumario_por_populacao.csv"),
          row.names = FALSE)

# 2. Variantes fine-mapped em múltiplas populações
multi_pop <- genotipadas %>%
  group_by(.data[[COL_RSID]]) %>%
  filter(n() > 1) %>%
  arrange(.data[[COL_RSID]], .data[[COL_POP]])

cat(sprintf("\nVariantes em múltiplas populações: %d\n",
            length(unique(multi_pop[[COL_RSID]]))))
write.csv(multi_pop, file.path(OUT_TABELAS, "variantes_multi_pop.csv"),
          row.names = FALSE)

# 3. Matriz de Z-score por população (direção do efeito)
if (COL_Z_EXOMA %in% colnames(genotipadas)) {
  beta_matrix <- genotipadas %>%
    select(all_of(c(COL_RSID, COL_POP, COL_Z_EXOMA))) %>%
    distinct() %>%
    pivot_wider(names_from = all_of(COL_POP), values_from = all_of(COL_Z_EXOMA))
  write.csv(beta_matrix, file.path(OUT_TABELAS, "matriz_z_por_populacao.csv"),
            row.names = FALSE)
}

cat("\nOK. Análise por ancestralidade salva em ", OUT_TABELAS, "\n", sep = "")
