#!/usr/bin/env Rscript
# 04_prepara_glimpse.R
# Prepara os inputs do passo 04 (imputação GLIMPSE2 dos sinais MVP no exoma):
#   1) janelas por cromossomo (± WINDOW_MB ao redor de cada sinal, fundidas)
#   2) lista de BAMs recalibrated (Sarek) dos casos, com ID de amostra
#
# Input :
#   - dados/mvp/bc_signal.csv            (output do 01_filtrar_mvp.R)
#   - BAMs recalibrated do Sarek (path no CONFIG)
# Output: resultados/glimpse/regioes.txt (CHR\tstart\tend, janelas fundidas)
#         resultados/glimpse/sites.txt   (CHR:POS exatos dos sinais)
#         resultados/glimpse/bams.txt    (1 path por linha — GLIMPSE2_phase)
#         resultados/glimpse/bam_samples.tsv (path<TAB>sample, para conferência)

# ─────────────────────────── CONFIG ───────────────────────────
MVP_SIGNAL <- "dados/mvp/bc_signal.csv"
OUT_GLIMPSE <- "resultados/glimpse"

WINDOW_MB <- 1        # meia-janela (bp) ao redor de cada sinal → janela total = 2x
BAM_DIR <- "/storage4/matheusbomfim/programas/sarek/results_CA_Mama_DBGAP_New/variant_calling/preprocessing/recalibrated"
BAM_PATTERN <- ".*\\.recal\\.bam$"   # regex; ajuste conforme os nomes reais
BAM_RECURSIVE <- TRUE

# Colunas do bc_signal.csv
MVP_CHR  <- "CHR"
MVP_BP   <- "BP38"     # posição hg38
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
MVP_SIGNAL <- abs_path(MVP_SIGNAL)
OUT_GLIMPSE <- abs_path(OUT_GLIMPSE)

suppressPackageStartupMessages(library(data.table))
if (!dir.exists(OUT_GLIMPSE)) dir.create(OUT_GLIMPSE, recursive = TRUE)

cat("== 04_prepara_glimpse.R ==\n")
mvp <- fread(MVP_SIGNAL, na.strings = c("NA", "", "N/A"))
stopifnot(all(c(MVP_CHR, MVP_BP) %in% names(mvp)))

mvp[, CHR := gsub("^chr", "", as.character(get(MVP_CHR)), ignore.case = TRUE)]
mvp[, BP  := as.numeric(get(MVP_BP))]
mvp <- mvp[!is.na(BP) & !is.na(CHR)]
if (!any(grepl("^[0-9]+$", mvp$CHR))) stop("CHR não-numérico no bc_signal.csv")

# 1. Janelas por cromossomo (fundir sobreposições)
win <- mvp[, .(start = min(BP) - WINDOW_MB, end = max(BP) + WINDOW_MB), by = CHR]
win <- win[order(as.integer(CHR), start)]
merged <- win[, {
  s <- start[1]; e <- end[1]
  out <- list()
  for (i in seq_along(start)) {
    if (start[i] > e) { out[[length(out) + 1]] <- c(s, e); s <- start[i]; e <- end[i] }
    else e <- max(e, end[i])
  }
  out[[length(out) + 1]] <- c(s, e)
  as.data.table(do.call(rbind, out))
}, by = CHR]
setnames(merged, c("CHR", "start", "end"))
fwrite(merged, file.path(OUT_GLIMPSE, "regioes.txt"), sep = "\t", col.names = FALSE)
cat(sprintf("Janelas: %d (cobrindo %d sinais)\n", nrow(merged), nrow(mvp)))

# 2. Sítios exatos dos sinais (chr:pos)
sites <- mvp[, paste0("chr", CHR, ":", BP)]
writeLines(sites, file.path(OUT_GLIMPSE, "sites.txt"))
cat(sprintf("Sites: %d\n", length(sites)))

# 3. BAMs recalibrated (GLIMPSE2_phase exige 1 path por linha)
if (!dir.exists(BAM_DIR)) {
  warning("BAM_DIR não encontrado: ", BAM_DIR, " — bams.txt não gerado.")
  quit(status = 1)
}
bams <- list.files(BAM_DIR, pattern = BAM_PATTERN, recursive = BAM_RECURSIVE, full.names = TRUE)
bams <- grep("\\.bai$", bams, value = TRUE, invert = TRUE)
if (length(bams) == 0) {
  warning("Nenhum BAM com o padrão '", BAM_PATTERN, "' em ", BAM_DIR)
  quit(status = 1)
}
bam_tbl <- data.table(path = normalizePath(bams),
                      sample = sub("\\.[^.]*$", "", basename(bams)))
# Lista de paths (coluna única) para o GLIMPSE2_phase --bam-list
writeLines(bam_tbl$path, file.path(OUT_GLIMPSE, "bams.txt"))
# Mapa path<TAB>sample (referência; IDs reais vêm do @RG do BAM)
fwrite(bam_tbl, file.path(OUT_GLIMPSE, "bam_samples.tsv"), sep = "\t", col.names = FALSE)
cat(sprintf("BAMs: %d (padrão %s)\n", nrow(bam_tbl), BAM_PATTERN))
cat("OK. Inputs prontos em", OUT_GLIMPSE, "\n")
