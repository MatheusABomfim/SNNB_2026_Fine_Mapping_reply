#!/usr/bin/env Rscript
# 01_filtrar_mvp.R
# Filtrar o fine-mapping do MVP (Verma et al. 2024, Science) para:
#   1) câncer de mama
#   2) variantes codantes (VEP annotation)
#
# Input : fine_mapping_GCST90479802.xlsx  (ou GCST90479802.tsv.gz)
# Output: dados/mvp/bc_coding.csv, dados/mvp/bc_todos_sinais.csv,
#         dados/mvp/rsids_coding.txt

# ─────────────────────────── CONFIG ───────────────────────────
# Ajuste os caminhos e o mapeamento de colunas conforme os headers reais.

MVP_INPUT <- "dados/mvp/fine_mapping_GCST90479802.xlsx"  # ou ".tsv.gz"
# Se usar o .tsv.gz (delimitado por tab), mude para TRUE:
MVP_IS_TSV_GZ <- FALSE

# Mapeamento de colunas do MVP (preencher com os nomes reais).
# Use "" se a coluna não existir no arquivo.
COL_TRAIT     <- "Trait"                # nome do trait/doença (buscar "breast")
COL_RSID      <- "RSID"                 # identificador rs
COL_VEP       <- "VEP Annotation"       # tipo de variante (VEP consequence)
COL_PIP       <- "Overall PIP"          # posterior inclusion probability
COL_POP       <- "Population"           # AFR / AMR / EUR / EAS
COL_BETA      <- "Beta Population"      # efeito na população
COL_SE        <- "SE Population"        # erro padrão
COL_EAF       <- "EAF Population"       # frequência do alelo de efeito
COL_LOCUS     <- "Locus"                # região genômica
COL_NOVEL     <- "Previously Unidentified if High PIP"  # novo ou conhecido

# Tipos de variante considerados codantes
CODING_TYPES <- c(
  "missense_variant",
  "synonymous_variant",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "splice_donor_variant",
  "splice_acceptor_variant",
  "inframe_deletion",
  "inframe_insertion",
  "protein_altering_variant"
)

OUT_DIR <- "dados/mvp"
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})
if (!MVP_IS_TSV_GZ) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Pacote 'readxl' necessário para ler .xlsx. Instale com install.packages('readxl').")
  }
  library(readxl)
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 1. Carregar MVP
cat("== 01_filtrar_mvp.R ==\n")
cat(sprintf("Lendo MVP de: %s\n", MVP_INPUT))
if (MVP_IS_TSV_GZ) {
  mvp <- read.delim(MVP_INPUT, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "#")
} else {
  mvp <- read_excel(MVP_INPUT)
  mvp <- as.data.frame(mvp)
}

cat(sprintf("Total de sinais no MVP: %d\n", nrow(mvp)))
cat("Colunas disponíveis:\n")
print(colnames(mvp))
if (COL_RSID != "" && !COL_RSID %in% colnames(mvp)) {
  warning(sprintf("Coluna '%s' (COL_RSID) não encontrada. Verifique o mapeamento.", COL_RSID))
}

# 2. Filtrar câncer de mama
bc <- mvp
if (COL_TRAIT != "" && COL_TRAIT %in% colnames(mvp)) {
  bc <- mvp %>% filter(str_detect(tolower(.data[[COL_TRAIT]]), "breast"))
}
cat(sprintf("Sinais de câncer de mama: %d\n", nrow(bc)))

# 3. Filtrar variantes codantes
if (COL_VEP != "" && COL_VEP %in% colnames(bc)) {
  vep_col <- bc[[COL_VEP]]
  # VEP annotation pode vir como string única ou múltiplas (CSQ); detectar coding
  bc_coding <- bc[vapply(vep_col, function(x) {
    any(CODING_TYPES %in% strsplit(as.character(x), "[&,]")[[1]])
  }, logical(1)), , drop = FALSE]
} else {
  warning(sprintf("Coluna '%s' (COL_VEP) não encontrada. Mantendo todos os sinais.", COL_VEP))
  bc_coding <- bc
}
cat(sprintf("Coding variants no câncer de mama: %d\n", nrow(bc_coding)))

# 4. Estatísticas descritivas
cat("\n--- Por tipo de VEP ---\n")
if (COL_VEP %in% colnames(bc_coding)) print(table(bc_coding[[COL_VEP]]))
cat("\n--- Por população ---\n")
if (COL_POP %in% colnames(bc_coding)) print(table(bc_coding[[COL_POP]]))
cat("\n--- Por PIP > 0.8 ---\n")
if (COL_PIP %in% colnames(bc_coding)) {
  pip <- suppressWarnings(as.numeric(bc_coding[[COL_PIP]]))
  print(table(pip > 0.8, useNA = "ifany"))
}

# 5. Salvar
write.csv(bc_coding, file.path(OUT_DIR, "bc_coding.csv"), row.names = FALSE)
write.csv(bc, file.path(OUT_DIR, "bc_todos_sinais.csv"), row.names = FALSE)

# 6. Lista de RSIDs para match
if (COL_RSID %in% colnames(bc_coding)) {
  writeLines(unique(as.character(bc_coding[[COL_RSID]])),
             file.path(OUT_DIR, "rsids_coding.txt"))
}

cat("\nOK. Arquivos gerados em ", OUT_DIR, "\n", sep = "")
