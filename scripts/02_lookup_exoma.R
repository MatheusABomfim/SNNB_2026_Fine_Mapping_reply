#!/usr/bin/env Rscript
# 02_lookup_exoma.R
# Look-up das coding variants do MVP nas summary statistics do exoma phs000822.
# Classifica cada variante quanto à consistência da direção do efeito:
#   Replicada_p0.05 / Direcao_consistente / Discordante / Nao_genotipada
#
# Input :
#   - dados/mvp/bc_signal.csv                (output do 01_filtrar_mvp.R)
#   - dados/exoma/gwas_dbgap_+vep.tsv         (sumstats exoma com VEP)
# Output: resultados/tabelas/lookup_completo.csv, resultados/tabelas/top_replicadas.csv

# ─────────────────────────── CONFIG ───────────────────────────
MVP_SIGNAL <- "dados/mvp/bc_signal.csv"
EXOMA_SUMSTATS <- "dados/exoma/gwas_dbgap_sem_filtragem+joint_call+vep.tsv"
OUT_TABELAS <- "resultados/tabelas"

# Mapeamento de colunas do MVP (do bc_signal.csv)
MVP_RSID  <- "RSID"
MVP_POP   <- "Population"
MVP_PIP   <- "Overall PIP"
MVP_BETA  <- "Beta Population"
MVP_SE    <- "SE Population"

# Mapeamento de colunas do exoma (preencher conforme headers reais)
EX_RSID  <- "rsid"
EX_CHR   <- "chr"
EX_BP    <- "bp"
EX_A1    <- "a1"
EX_A2    <- "a2"
EX_MAF   <- "maf"
EX_BETA  <- "beta"
EX_SE    <- "se"
EX_PVAL  <- "pval"
EX_N     <- "n"

# Threshold de significância para "replicada"
PVAL_THRESHOLD <- 0.05
# ───────────────────────────────────────────────────────────────

# Caminhos relativos à raiz do repo (o pai do dir deste script), para o script
# funcionar rodado de qualquer CWD (ex.: via PBS em scripts/).
args0 <- commandArgs(trailingOnly = FALSE)
fidx <- grep("^--file=", args0)
SCRIPT_DIR <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
MVP_SIGNAL    <- file.path(REPO_ROOT, MVP_SIGNAL)
EXOMA_SUMSTATS <- file.path(REPO_ROOT, EXOMA_SUMSTATS)
OUT_TABELAS   <- file.path(REPO_ROOT, OUT_TABELAS)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
if (!dir.exists(OUT_TABELAS)) dir.create(OUT_TABELAS, recursive = TRUE)

# 1. Carregar dados
cat("== 02_lookup_exoma.R ==\n")
mvp_signal <- read.csv(MVP_SIGNAL, stringsAsFactors = FALSE, check.names = FALSE)
exoma <- read.delim(EXOMA_SUMSTATS, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "#")

cat(sprintf("Coding variants MVP: %d\n", nrow(mvp_signal)))
cat(sprintf("SNPs no exoma: %d\n", nrow(exoma)))
cat("Colunas do exoma:\n")
print(colnames(exoma))

# 2. Padronizar nomes do exoma (colunas reais → nomes canônicos)
exoma <- exoma %>%
  rename_with(~ ifelse(.x %in% EX_RSID, "RSID", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_CHR, "CHR", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_BP, "BP", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_A1, "A1", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_A2, "A2", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_MAF, "MAF_exoma", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_BETA, "BETA_exoma", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_SE, "SE_exoma", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_PVAL, "PVAL_exoma", .x)) %>%
  rename_with(~ ifelse(.x %in% EX_N, "N_exoma", .x))

# 3. Match por RSID (left_join: mantém todas as coding do MVP)
merged <- mvp_signal %>%
  left_join(select(exoma, any_of(c("RSID", "CHR", "BP", "A1", "A2",
                                   "MAF_exoma", "BETA_exoma", "SE_exoma",
                                   "PVAL_exoma", "N_exoma"))),
            by = "RSID")

# 4. Classificar
merged <- merged %>%
  mutate(
    Z_exoma = ifelse(is.na(SE_exoma) | SE_exoma == 0, NA,
                     BETA_exoma / SE_exoma),
    Z_MVP = ifelse(is.na(.data[[MVP_SE]]) | .data[[MVP_SE]] == 0, NA,
                   .data[[MVP_BETA]] / .data[[MVP_SE]]),
    classificacao = case_when(
      is.na(BETA_exoma) ~ "Nao_genotipada",
      is.na(.data[[MVP_BETA]]) ~ "Indeterminada",
      sign(BETA_exoma) == sign(.data[[MVP_BETA]]) &
        PVAL_exoma < PVAL_THRESHOLD ~ "Replicada_p0.05",
      sign(BETA_exoma) == sign(.data[[MVP_BETA]]) &
        PVAL_exoma >= PVAL_THRESHOLD ~ "Direcao_consistente",
      sign(BETA_exoma) != sign(.data[[MVP_BETA]]) ~ "Discordante",
      TRUE ~ "Indeterminada"
    )
  )

# 5. Tabela de contingência
cat("\n--- Classificação ---\n")
print(table(merged$classificacao, useNA = "ifany"))

# 6. Salvar lookup completo
write.csv(merged, file.path(OUT_TABELAS, "lookup_completo.csv"), row.names = FALSE)

# 7. Totais por população
cat("\n--- Por população ---\n")
if (MVP_POP %in% colnames(merged)) {
  pop_summary <- merged %>%
    count(.data[[MVP_POP]], classificacao) %>%
    pivot_wider(names_from = classificacao, values_from = n, values_fill = 0)
  print(pop_summary)
  write.csv(pop_summary, file.path(OUT_TABELAS, "lookup_por_populacao.csv"),
            row.names = FALSE)
}

# 8. Top replicadas (PIP > 0.8)
if (MVP_PIP %in% colnames(merged)) {
  top <- merged %>%
    filter(classificacao %in% c("Replicada_p0.05", "Direcao_consistente"),
           suppressWarnings(as.numeric(.data[[MVP_PIP]])) > 0.8) %>%
    arrange(desc(suppressWarnings(as.numeric(.data[[MVP_PIP]]))))
  cat(sprintf("\nTop replicadas com PIP > 0.8: %d\n", nrow(top)))
  write.csv(top, file.path(OUT_TABELAS, "top_replicadas.csv"), row.names = FALSE)
}

cat("\nOK. Look-up salvo em ", OUT_TABELAS, "/lookup_completo.csv\n", sep = "")
