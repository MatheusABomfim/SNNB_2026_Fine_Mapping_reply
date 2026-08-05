#!/usr/bin/env Rscript
# 04_tabelas_finais.R
# Monta as tabelas principais do artigo a partir do lookup e do MAGMA.
#
# Input :
#   - resultados/tabelas/lookup_completo.csv   (output do 02_lookup_exoma.R)
#   - resultados/magma/exoma.genes.out         (MAGMA, colaborador — opcional)
#   - dados/exoma/gwas_dbgap_+vep.tsv          (sumstats exoma)
# Output: resultados/tabelas/tabela{1,2,3,4}_*.csv

# ─────────────────────────── CONFIG ───────────────────────────
LOOKUP_CSV <- "resultados/tabelas/lookup_completo.csv"
MAGMA_OUT  <- "resultados/magma/exoma.genes.out"   # opcional
EXOMA_SUMSTATS <- "dados/exoma/gwas_dbgap_sem_filtragem+joint_call+vep.tsv"
OUT_TABELAS <- "resultados/tabelas"

# Colunas do lookup
COL_PIP  <- "Overall PIP"
COL_BETA <- "Beta Population"
COL_VEP  <- "VEP Annotation"
COL_POP  <- "Population"
COL_NOVEL <- "Previously Unidentified if High PIP"

# Colunas do exoma
EX_RSID <- "rsid"; EX_CHR <- "chr"; EX_BP <- "bp"; EX_MAF <- "maf"
EX_BETA <- "beta"; EX_PVAL <- "pval"
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(dplyr))
if (!dir.exists(OUT_TABELAS)) dir.create(OUT_TABELAS, recursive = TRUE)

lookup <- read.csv(LOOKUP_CSV, stringsAsFactors = FALSE, check.names = FALSE)
cat("== 04_tabelas_finais.R ==\n")

# ─── Tabela 1: Coding variants replicadas / direção consistente ───
tabela1 <- lookup %>%
  filter(classificacao %in% c("Replicada_p0.05", "Direcao_consistente")) %>%
  arrange(desc(suppressWarnings(as.numeric(.data[[COL_PIP]])))) %>%
  select(RSID, any_of(c(COL_VEP, COL_PIP, COL_BETA,
                        "BETA_exoma", "PVAL_exoma", "SE_exoma", COL_POP, COL_NOVEL)))
write.csv(tabela1, file.path(OUT_TABELAS, "tabela1_replicadas.csv"), row.names = FALSE)

# ─── Tabela 2: Discordantes ───
tabela2 <- lookup %>%
  filter(classificacao == "Discordante") %>%
  arrange(PVAL_exoma) %>%
  select(RSID, any_of(c(COL_VEP, COL_PIP, COL_BETA,
                        "BETA_exoma", "PVAL_exoma", COL_POP)))
write.csv(tabela2, file.path(OUT_TABELAS, "tabela2_discordantes.csv"), row.names = FALSE)

# ─── Tabela 3: SNPs do exoma não no MVP (novos candidatos) ───
mvp_rsids <- unique(lookup$RSID)
exoma <- read.delim(EXOMA_SUMSTATS, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "#")
# Padronizar nomes de colunas do exoma
rename_map <- c(EX_RSID = "RSID", EX_CHR = "CHR", EX_BP = "BP",
                EX_MAF = "MAF", EX_BETA = "BETA", EX_PVAL = "PVAL")
exoma <- exoma %>% rename(any_of(rename_map))

tabela3 <- exoma %>%
  filter(PVAL < 0.01, !RSID %in% mvp_rsids) %>%
  arrange(PVAL) %>%
  select(any_of(c("RSID", "CHR", "BP", "MAF", "BETA", "PVAL")))
write.csv(tabela3, file.path(OUT_TABELAS, "tabela3_novos_candidatos.csv"), row.names = FALSE)

# ─── Tabela 4: Genes do MAGMA (se disponível) ───
n_magma <- 0
if (file.exists(MAGMA_OUT)) {
  magma <- read.table(MAGMA_OUT, header = TRUE, stringsAsFactors = FALSE)
  tabela4 <- magma %>%
    filter(P < 0.01) %>%
    arrange(P) %>%
    select(any_of(c("GENE", "CHR", "ZSTAT", "P")))
  write.csv(tabela4, file.path(OUT_TABELAS, "tabela4_magma_genes.csv"), row.names = FALSE)
  n_magma <- nrow(tabela4)
} else {
  message(sprintf("MAGMA output não encontrado em %s — tabela 4 pulada.", MAGMA_OUT))
}

# ─── Resumo estatístico ───
n_genotipadas <- sum(lookup$classificacao != "Nao_genotipada")
n_consistentes <- sum(lookup$classificacao %in% c("Replicada_p0.05", "Direcao_consistente"))
n_replicadas <- sum(lookup$classificacao == "Replicada_p0.05")
n_discordantes <- sum(lookup$classificacao == "Discordante")

cat("\n=== RESUMO FINAL ===\n")
cat(sprintf("Total coding variants no MVP (breast cancer): %d\n", nrow(lookup)))
cat(sprintf("Genotipadas no exoma: %d (%.1f%%)\n",
            n_genotipadas, 100 * n_genotipadas / max(1, nrow(lookup))))
cat(sprintf("Direcao consistente: %d (%.1f%% das genotipadas)\n",
            n_consistentes, 100 * n_consistentes / max(1, n_genotipadas)))
cat(sprintf("Replicadas com p<0.05: %d\n", n_replicadas))
cat(sprintf("Discordantes: %d (%.1f%% das genotipadas)\n",
            n_discordantes, 100 * n_discordantes / max(1, n_genotipadas)))
cat(sprintf("Genes com p<0.01 no MAGMA exoma: %d\n", n_magma))
cat(sprintf("SNPs do exoma nao no MVP (p<0.01): %d\n", nrow(tabela3)))

cat("\nOK. Tabelas salvas em ", OUT_TABELAS, "\n", sep = "")
