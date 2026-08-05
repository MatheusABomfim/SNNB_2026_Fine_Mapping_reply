#!/usr/bin/env Rscript
# 01_filtrar_mvp.R
# Filtrar o fine-mapping do MVP (Verma et al. 2024, Science) para:
#   1) câncer de mama (coluna Description = "Biological Mother: Cancer, Breast")
#
# Input : fine_mapping_GCST90479802.xlsx  (ou GCST90479802.tsv.gz)
# Output: dados/mvp/bc_signal.csv, dados/mvp/rsids.txt

# ─────────────────────────── CONFIG ───────────────────────────
# Ajuste os caminhos e o mapeamento de colunas conforme os headers reais.

MVP_INPUT <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/fine_mapping_GCST90479802.xlsx"  # ou ".tsv.gz"
# Se usar o .tsv.gz (delimitado por tab), mude para TRUE:
MVP_IS_TSV_GZ <- FALSE

# Mapeamento de colunas do MVP (preencher com os nomes reais).
# Use "" se a coluna não existir no arquivo.
COL_TRAIT     <- "Trait"                # nome do trait/doença (buscar "breast")
COL_DESC      <- "Description"          # descrição do fenótipo (ex.: "Biological Mother:  Cancer, Breast")
COL_RSID      <- "RSID"                 # identificador rs
COL_VEP       <- "VEP Annotation"       # tipo de variante (VEP consequence)
COL_PIP       <- "Overall PIP"          # posterior inclusion probability
COL_POP       <- "Population"           # AFR / AMR / EUR / EAS
COL_BETA      <- "Beta Population"      # efeito na população
COL_SE        <- "SE Population"        # erro padrão
COL_EAF       <- "EAF Population"       # frequência do alelo de efeito
COL_LOCUS     <- "Locus"                # região genômica
COL_NOVEL     <- "Previously Unidentified if High PIP"  # novo ou conhecido

OUT_DIR <- "dados/mvp"
# ───────────────────────────────────────────────────────────────

# Caminhos relativos à raiz do repo (a raiz é o pai do dir deste script),
# para o script funcionar rodado de qualquer CWD (ex.: via PBS em scripts/).
# Caminhos absolutos são mantidos como estão.
abs_path <- function(p) {
  if (grepl("^/", p)) p else file.path(REPO_ROOT, p)
}
args0 <- commandArgs(trailingOnly = FALSE)
fidx <- grep("^--file=", args0)
SCRIPT_DIR <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
MVP_INPUT <- abs_path(MVP_INPUT)
OUT_DIR   <- abs_path(OUT_DIR)

# Detectar formato de entrada: se readxl não estiver disponível e o .xlsx não
# existir, tenta o .tsv.gz (mesmo arquivo em formato tab separado, sem readxl).
have_readxl <- requireNamespace("readxl", quietly = TRUE)
MVP_IS_TSV_GZ <- MVP_IS_TSV_GZ || !file.exists(MVP_INPUT)
if (MVP_IS_TSV_GZ) {
  candidates <- c(sub("\\.xlsx$", ".tsv.gz", MVP_INPUT),
                  sub("fine_mapping_", "", sub("\\.xlsx$", ".tsv.gz", MVP_INPUT)))
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("MVP não encontrado. Nem .xlsx (", MVP_INPUT, ") nem .tsv.gz (",
         paste(candidates, collapse = ", "), ") existem — confira dados/mvp/. Se for usar .xlsx, instale readxl (micromamba install -n quali -c conda-forge r-readxl).")
  }
  MVP_INPUT <- found[1]
} else if (!have_readxl) {
  stop("Pacote 'readxl' necessário para ler .xlsx e nenhum .tsv.gz encontrado. Instale com micromamba install -n quali -c conda-forge r-readxl.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})
if (!MVP_IS_TSV_GZ) library(readxl)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 1. Carregar MVP
cat("== 01_filtrar_mvp.R ==\n")
cat(sprintf("Lendo MVP de: %s\n", MVP_INPUT))
if (MVP_IS_TSV_GZ) {
  mvp <- read.delim(MVP_INPUT, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "#")
} else {
  mvp <- read_excel(MVP_INPUT)
  # O arquivo tem uma linha de título antes do header real ("Data S1. ..."),
  # então se a coluna de RSID não aparecer, o header está na linha 2.
  if (COL_RSID != "" && !COL_RSID %in% colnames(mvp)) {
    mvp <- read_excel(MVP_INPUT, skip = 1)
  }
  mvp <- as.data.frame(mvp)
}

cat(sprintf("Total de sinais no MVP: %d\n", nrow(mvp)))
cat("Colunas disponíveis:\n")
print(colnames(mvp))
if (COL_RSID != "" && !COL_RSID %in% colnames(mvp)) {
  warning(sprintf("Coluna '%s' (COL_RSID) não encontrada. Verifique o mapeamento.", COL_RSID))
}

# 2. Filtrar câncer de mama (coluna Description = fenótipo da mãe)
bc <- mvp
if (COL_DESC != "" && COL_DESC %in% colnames(mvp)) {
  bc <- mvp %>% filter(str_detect(tolower(.data[[COL_DESC]]),
                                  "biological mother\\s*:\\s*cancer\\s*,\\s*breast"))
} else if (COL_TRAIT != "" && COL_TRAIT %in% colnames(mvp)) {
  bc <- mvp %>% filter(str_detect(tolower(.data[[COL_TRAIT]]), "breast"))
} else {
  warning("Nenhuma coluna de fenótipo (Description/Trait) encontrada. Mantendo todos os sinais.")
}
cat(sprintf("Sinais de câncer de mama (Biological Mother): %d\n", nrow(bc)))

# 3. Sem filtro por variantes codantes — mantém todos os sinais de mãe-mama
#    (a curadoria por p-value/priorização é feita depois).
cat(sprintf("Sinais mantidos (sem filtro coding): %d\n", nrow(bc)))

# 4. Estatísticas descritivas
cat("\n--- Por tipo de VEP ---\n")
if (COL_VEP %in% colnames(bc)) print(table(bc[[COL_VEP]]))
cat("\n--- Por população ---\n")
if (COL_POP %in% colnames(bc)) print(table(bc[[COL_POP]]))
cat("\n--- Por PIP > 0.8 ---\n")
if (COL_PIP %in% colnames(bc)) {
  pip <- suppressWarnings(as.numeric(bc[[COL_PIP]]))
  print(table(pip > 0.8, useNA = "ifany"))
}

# 5. Salvar
write.csv(bc, file.path(OUT_DIR, "bc_signal.csv"), row.names = FALSE)

# 6. Lista de RSIDs para match
if (COL_RSID %in% colnames(bc)) {
  writeLines(unique(as.character(bc[[COL_RSID]])),
             file.path(OUT_DIR, "rsids.txt"))
}

cat("\nOK. Arquivos gerados em ", OUT_DIR, "\n", sep = "")
