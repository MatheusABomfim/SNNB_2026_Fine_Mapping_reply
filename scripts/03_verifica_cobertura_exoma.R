#!/usr/bin/env Rscript
# 03_verifica_cobertura_exoma.R
# Verifica se as regiões dos sinais fine-mapped do MVP (Biological Mother:
# Cancer, Breast) são cobertas pelo exoma phs000822, medindo os SNPs do exoma
# dentro de janelas crescentes ao redor de cada sinal (posição hg38 = BP38).
#
# Input :
#   - dados/mvp/bc_signal.csv                  (output do 01_filtrar_mvp.R)
#   - gwas/gwas_prod/exoma_sumstats.txt        (sumstats exoma, output do 05.5)
# Output: resultados/tabelas/cobertura_exoma.csv + resumo no log

# ─────────────────────────── CONFIG ───────────────────────────
MVP_SIGNAL <- "dados/mvp/bc_signal.csv"
EXOMA_SUMSTATS <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/gwas_prod/exoma_sumstats.txt"
OUT_TABELAS <- "resultados/tabelas"

# Janelas (bp) ao redor do sinal index para contar SNPs do exoma
WINDOWS <- c(10000, 100000, 250000, 500000, 1000000)
WIN_NAMES <- paste0("n_snps_", formatC(WINDOWS, format = "d"))

# Colunas
MVP_CHR  <- "CHR"
MVP_BP   <- "BP38"      # posição hg38 do sinal MVP
MVP_RSID <- "RSID"
MVP_LOCUS <- "Locus"
# ───────────────────────────────────────────────────────────────

# Caminhos relativos à raiz do repo (o pai do dir deste script), para o script
# funcionar rodado de qualquer CWD (ex.: via PBS em scripts/).
# Caminhos absolutos são mantidos como estão.
abs_path <- function(p) {
  if (grepl("^/", p)) p else file.path(REPO_ROOT, p)
}
args0 <- commandArgs(trailingOnly = FALSE)
fidx <- grep("^--file=", args0)
SCRIPT_DIR <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
MVP_SIGNAL   <- abs_path(MVP_SIGNAL)
EXOMA_SUMSTATS <- abs_path(EXOMA_SUMSTATS)
OUT_TABELAS  <- abs_path(OUT_TABELAS)

suppressPackageStartupMessages(library(data.table))
if (!dir.exists(OUT_TABELAS)) dir.create(OUT_TABELAS, recursive = TRUE)

cat("== 03_verifica_cobertura_exoma.R ==\n")
mvp <- fread(MVP_SIGNAL)
exoma <- fread(EXOMA_SUMSTATS)

# 1. Normalizar colunas do exoma (05.5 grava MAIÚSCULAS) e chr sem prefixo
setnames(exoma, tolower(names(exoma)))
if (!all(c("chr", "bp") %in% names(exoma))) {
  stop("exoma_sumstats.txt sem colunas chr/bp: ", paste(names(exoma), collapse = ", "))
}
exoma[, chr := gsub("^chr", "", as.character(chr), ignore.case = TRUE)]
exoma[, bp := as.numeric(bp)]

# 2. Normalizar colunas do MVP
if (!all(c(MVP_CHR, MVP_BP) %in% names(mvp))) {
  stop("bc_signal.csv sem colunas ", MVP_CHR, "/", MVP_BP,
       ". Colunas: ", paste(names(mvp), collapse = ", "))
}
mvp[, chr_mvp := gsub("^chr", "", as.character(get(MVP_CHR)), ignore.case = TRUE)]
mvp[, bp_mvp := as.numeric(get(MVP_BP))]

# 3. Índice do exoma por chr (para contagem por janela)
setkey(exoma, chr)

mvp2 <- data.table(
  RSID  = if (MVP_RSID %in% names(mvp)) as.character(mvp[[MVP_RSID]]) else rep(NA_character_, nrow(mvp)),
  CHR   = mvp$chr_mvp,
  BP    = mvp$bp_mvp,
  Locus = if (MVP_LOCUS %in% names(mvp)) as.character(mvp[[MVP_LOCUS]]) else rep(NA_character_, nrow(mvp))
)
mvp2[, row_id := .I]

cobertura <- mvp2[, {
  ex_bp <- exoma[.(CHR)]$bp
  dist_min <- if (length(ex_bp) == 0) NA_real_ else min(abs(ex_bp - BP))
  pos_exata <- any(ex_bp == BP)
  counts <- vapply(WINDOWS, function(w) sum(abs(ex_bp - BP) <= w), integer(1))
  c(list(
    RSID  = RSID,
    CHR   = CHR,
    BP38  = BP,
    Locus = Locus,
    dist_snp_min_bp = dist_min,
    posicao_exata_no_exoma = pos_exata
  ), setNames(as.list(counts), WIN_NAMES))
}, by = row_id][, row_id := NULL]

cat(sprintf("Sinais MVP: %d | SNPs exoma: %d\n", nrow(mvp), nrow(exoma)))
cat(sprintf("Posição exata presente no exoma: %d/%d\n",
            sum(cobertura$posicao_exata_no_exoma), nrow(cobertura)))

# 4. Resumo por janela
cat("\n--- Cobertura por janela ---\n")
for (i in seq_along(WINDOWS)) {
  col <- WIN_NAMES[i]
  n_cob <- sum(cobertura[[col]] > 0)
  cat(sprintf("  janela ±%s kb: %d/%d sinais com >=1 SNP do exoma\n",
              formatC(WINDOWS[i] / 1000, format = "d"), n_cob, nrow(cobertura)))
}
cat(sprintf("\nDistância mediana ao SNP exômico mais próximo: %s bp\n",
            ifelse(all(is.na(cobertura$dist_snp_min_bp)), "NA",
                   formatC(median(cobertura$dist_snp_min_bp, na.rm = TRUE), format = "d"))))

# 5. Sinais sem nenhum SNP do exoma mesmo em ±1 Mb
col_max <- WIN_NAMES[which.max(WINDOWS)]
uncovered <- cobertura[is.na(dist_snp_min_bp) | get(col_max) == 0]
if (nrow(uncovered) > 0) {
  cat("\n--- Sinais SEM cobertura do exoma (nem ±1 Mb) ---\n")
  print(uncovered[, .(RSID, CHR, BP38, Locus)])
} else {
  cat("\nTodos os sinais têm pelo menos 1 SNP do exoma em ±1 Mb.\n")
}

# 6. Salvar
fwrite(cobertura, file.path(OUT_TABELAS, "cobertura_exoma.csv"))
cat("\nSalvo: ", file.path(OUT_TABELAS, "cobertura_exoma.csv"), "\n", sep = "")
cat("OK.\n")
