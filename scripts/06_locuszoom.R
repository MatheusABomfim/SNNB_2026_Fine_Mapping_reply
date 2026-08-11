#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# 06_locuszoom.R
# Gera locus zoom plots (estilo LocusZoom) dos hits significativos do GWAS
# usando o pacote locuszoomr, SEM LD (pontos coloridos por significância).
#
# Inputs:
#   - GWAS assoc completo : /storage4/.../gwas/gwas_prod/SRR_gwas.assoc
#   - Hits significativos : resultados/tabelas/gwas_hits_significativos.csv
#     (se vazio/ausente, tenta gwas_hits_sugestivos.csv)
# Output:
#   - resultados/figuras/locuszoom_chr{chr}_{pos}.png  (um por locus)
#
# Dependências R (env quali):
#   install.packages("locuszoomr")
#   BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v112"))
#   # v112 = GRCh38 (hg38). Ajuste ENS_DB conforme o banco instalado.
# ──────────────────────────────────────────────────────────────────────────────

# ─────────────────────────── CONFIG ───────────────────────────
GWAS_ASSOC  <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/gwas_prod/SRR_gwas.assoc"
HITS_SIG    <- "resultados/tabelas/gwas_hits_significativos.csv"
HITS_SUG    <- "resultados/tabelas/gwas_hits_sugestivos.csv"
OUT_FIGURAS <- "resultados/figuras"

ENS_DB      <- "EnsDb.Hsapiens.v112"   # GRCh38 (hg38); precisa estar instalado
WINDOW_BP   <- 5e5                     # janela total (fix_window) centrada no SNP
PCUTOFF     <- 5e-8                    # limiar de significância p/ colorir pontos
SCHEME      <- c("grey", "dodgerblue", "red")  # normal, significativo, index
# ───────────────────────────────────────────────────────────────

# Caminhos relativos à raiz do repo (o pai do dir deste script), para o script
# funcionar rodado de qualquer CWD (ex.: via PBS em scripts/).
# Caminhos absolutos são mantidos como estão.
abs_path <- function(p) if (grepl("^/", p)) p else file.path(REPO_ROOT, p)
args0 <- commandArgs(trailingOnly = FALSE)
fidx <- grep("^--file=", args0)
SCRIPT_DIR <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))
HITS_SIG    <- abs_path(HITS_SIG)
HITS_SUG    <- abs_path(HITS_SUG)
OUT_FIGURAS <- abs_path(OUT_FIGURAS)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(locuszoomr))
suppressPackageStartupMessages(library(ENS_DB, character.only = TRUE))

if (!dir.exists(OUT_FIGURAS)) dir.create(OUT_FIGURAS, recursive = TRUE)
cat("== 06_locuszoom.R ==\n")
cat("GWAS assoc :", GWAS_ASSOC, "\n")
cat("EnsDB      :", ENS_DB, "\n")
cat("Janela     :", WINDOW_BP, "bp (fix_window)\n")
cat("Output     :", OUT_FIGURAS, "\n\n")

# ── 1. Carregar sumstats GWAS ──────────────────────────────────────────────────
if (!file.exists(GWAS_ASSOC)) stop("GWAS assoc não encontrado: ", GWAS_ASSOC)
gwas <- fread(GWAS_ASSOC)
setnames(gwas, tolower(names(gwas)))

required <- c("chr", "pos", "id", "p")
missing  <- setdiff(required, names(gwas))
if (length(missing) > 0) stop("Colunas obrigatórias ausentes no assoc: ",
                              paste(missing, collapse = ", "))
gwas <- gwas[, .(chr, pos, id, p)]
gwas[, p := suppressWarnings(as.numeric(p))]
gwas <- gwas[!is.na(p) & is.finite(p)]
cat(sprintf("SNPs no assoc: %d\n", nrow(gwas)))

# ── 2. Carregar hits (significativos, com fallback p/ sugestivos) ─────────────
read_hits <- function(f) {
  if (!file.exists(f)) return(NULL)
  h <- fread(f)
  if (nrow(h) == 0) return(NULL)
  setnames(h, tolower(names(h)))
  if (!all(c("chr", "pos", "id") %in% names(h))) return(NULL)
  unique(h[, .(chr, pos, id)])
}

hits <- read_hits(HITS_SIG)
if (is.null(hits)) {
  cat("AVISO:", HITS_SIG, "vazio/ausente — usando", HITS_SUG, "\n")
  hits <- read_hits(HITS_SUG)
}
if (is.null(hits) || nrow(hits) == 0) stop("Nenhum hit para plotar.")

cat(sprintf("Loci a plotar: %d\n", nrow(hits)))

# ── 3. Gerar um locus plot por hit ─────────────────────────────────────────────
n_ok <- 0
for (i in seq_len(nrow(hits))) {
  hit_id <- hits$id[i]
  label  <- sprintf("%s:%s", hits$chr[i], hits$pos[i])

  if (!hit_id %in% gwas$id) {
    warning(sprintf("[%d/%d] SNP index '%s' não encontrado no assoc — pulando.",
                    i, nrow(hits), hit_id))
    next
  }

  loc <- tryCatch(
    locus(data = gwas, index_snp = hit_id, fix_window = WINDOW_BP,
          ens_db = ENS_DB, chrom = "chr", pos = "pos", labs = "id", p = "p"),
    error = function(e) { warning(sprintf("locus() falhou para %s: %s", label,
                                          conditionMessage(e))); NULL }
  )
  if (is.null(loc)) next

  out_png <- file.path(OUT_FIGURAS,
                       sprintf("locuszoom_%s_%s.png",
                               sub("^chr", "", hits$chr[i]), hits$pos[i]))
  png(out_png, width = 2400, height = 1400, res = 300)
  locus_plot(loc, labels = "index", pcutoff = PCUTOFF, scheme = SCHEME)
  dev.off()

  n_snp <- if (is.null(loc$data)) 0 else nrow(loc$data)
  cat(sprintf("  [%d/%d] %s  (%d SNPs)  → %s\n",
              i, nrow(hits), label, n_snp, basename(out_png)))
  n_ok <- n_ok + 1
}

cat(sprintf("\nDone. %d locus plots em %s\n", n_ok, OUT_FIGURAS))
