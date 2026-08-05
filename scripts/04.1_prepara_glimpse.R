#!/usr/bin/env Rscript
# 04.1_prepara_glimpse.R
# Prepara os inputs do passo 04 (imputação GLIMPSE2 dos sinais MVP no exoma):
#   1) janelas por cromossomo (± WINDOW_MB ao redor de cada sinal, fundidas)
#   2) lista de CRAMs (Sarek, work/ do Nextflow) dos casos, 1 por amostra
#
# Input :
#   - dados/mvp/bc_signal.csv            (output do 01_filtrar_mvp.R)
#   - CRAMs *.md.cram do Sarek (path no CONFIG)
# Output: resultados/glimpse/regioes.txt (CHR\tstart\tend, janelas fundidas)
#         resultados/glimpse/sites.txt   (CHR:POS exatos dos sinais)
#         resultados/glimpse/bams.txt    (path<SPACE>sample — GLIMPSE2_phase)
#         resultados/glimpse/bam_samples.tsv (path<TAB>sample, para conferência)

# ─────────────────────────── CONFIG ───────────────────────────
MVP_SIGNAL <- "dados/mvp/bc_signal.csv"
OUT_GLIMPSE <- "resultados/glimpse"

WINDOW_MB <- 1        # meia-janela (bp) ao redor de cada sinal → janela total = 2x
# BAMs/CRAMs dos casos: saída do Sarek está no work/ do Nextflow (não há
# diretório 'preprocessing/recalibrated' publicado). Usamos os *.md.cram
# (pós-MarkDuplicates, 1 por amostra). O mesmo CRAM pode aparecer em mais de um
# hash do work/ — deduplicamos por nome de amostra abaixo.
BAM_DIR <- "/storage4/matheusbomfim/programas/sarek/work/nextflow_work"
BAM_PATTERN <- ".*\\.md\\.cram$"   # regex; ajuste conforme os nomes reais
BAM_RECURSIVE <- TRUE
# Suaviza o nome da amostra: GLIMPSE2 usa a 2ª coluna do bams.txt como sample
# ID no BCF de saída — precisa casar com os IIDs do GWAS (ex.: SRR2943808).
# Remove o sufixo de MarkDuplicates ('.md') do nome do arquivo.
SAMPLE_STRIP <- "\\.md$"   # regex do sufixo a remover do basename

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

cat("== 04.1_prepara_glimpse.R ==\n")
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

# 3. CRAMs dos casos (GLIMPSE2_phase exige 1 arquivo por linha; aceita CRAM)
if (!dir.exists(BAM_DIR)) {
  warning("BAM_DIR não encontrado: ", BAM_DIR, " — bams.txt não gerado.")
  quit(status = 1)
}
bams <- list.files(BAM_DIR, pattern = BAM_PATTERN, recursive = BAM_RECURSIVE, full.names = TRUE)
# descarta índices (.crai/.bai) e duplicatas do work/ do Nextflow (mesmo
# arquivo pode existir em mais de um hash); escolhe 1 path por amostra
bams <- bams[!grepl("\\.crai$|\\.bai$", bams)]
if (length(bams) == 0) {
  warning("Nenhum CRAM com o padrão '", BAM_PATTERN, "' em ", BAM_DIR)
  quit(status = 1)
}
bam_tbl <- data.table(path = bams, sample = basename(bams))
# SRR2943808.md.cram -> SRR2943808 (remove extensão + sufixo de MarkDuplicates)
bam_tbl[, sample := sub("\\.[^.]*$", "", sample)]
bam_tbl[, sample := sub(SAMPLE_STRIP, "", sample)]
setorder(bam_tbl, sample, path)
n_dups <- sum(duplicated(bam_tbl$sample))
if (n_dups > 0) {
  n_uniq <- length(unique(bam_tbl$sample))
  cat(sprintf("  aviso: %d duplicata(s) no work/ do Nextflow (%d amostra(s) únicas) — mantendo o 1º path por amostra\n",
              n_dups, n_uniq))
}
bam_tbl <- bam_tbl[!duplicated(sample)]
bam_tbl[, path := normalizePath(path)]
# GLIMPSE2_phase --bam-list: 1 arquivo por linha, 2ª coluna (espaço) = sample
fwrite(bam_tbl, file.path(OUT_GLIMPSE, "bams.txt"), sep = " ", col.names = FALSE)
# Mapa path<TAB>sample (referência)
fwrite(bam_tbl, file.path(OUT_GLIMPSE, "bam_samples.tsv"), sep = "\t", col.names = FALSE)
cat(sprintf("CRAMs: %d (padrão %s)\n", nrow(bam_tbl), BAM_PATTERN))
cat("OK. Inputs prontos em", OUT_GLIMPSE, "\n")
