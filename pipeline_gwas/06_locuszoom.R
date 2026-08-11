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
#   - resultados/figuras/locuszoom_chr{chr}_{pos}.png   (um por locus)
#   - resultados/tabelas/gwas_localzoom.tsv.gz + .tbi   (para LocalZoom/browser)
#
# Dependências R (env locuszoom, via micromamba):
#   micromamba create -n locuszoom -c conda-forge -c bioconda \
#       r-locuszoomr bioconductor-ensdb.hsapiens.v86 tabix
#   # v86 = GRCh38 (hg38), único EnsDb.hg38 no bioconda (v112 não existe)
#   # tabix traz bgzip + tabix (htslib) para indexar o arquivo do LocalZoom
# ──────────────────────────────────────────────────────────────────────────────

# ─────────────────────────── CONFIG ───────────────────────────
GWAS_ASSOC  <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/gwas_prod/SRR_gwas.assoc"
HITS_SIG    <- "resultados/tabelas/gwas_hits_significativos.csv"
HITS_SUG    <- "resultados/tabelas/gwas_hits_sugestivos.csv"
OUT_FIGURAS <- "resultados/figuras"
OUT_TABELAS <- "resultados/tabelas"
LZ_BASE     <- "gwas_localzoom"          # base name do arquivo (tsv.gz + .tbi)

ENS_DB      <- "EnsDb.Hsapiens.v86"   # GRCh38 (hg38); via bioconda
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
OUT_TABELAS <- abs_path(OUT_TABELAS)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(locuszoomr))
suppressPackageStartupMessages(library(ENS_DB, character.only = TRUE))

if (!dir.exists(OUT_FIGURAS)) dir.create(OUT_FIGURAS, recursive = TRUE)
if (!dir.exists(OUT_TABELAS)) dir.create(OUT_TABELAS, recursive = TRUE)
cat("== 06_locuszoom.R ==\n")
cat("GWAS assoc :", GWAS_ASSOC, "\n")
cat("EnsDB      :", ENS_DB, "\n")
cat("Janela     :", WINDOW_BP, "bp (fix_window)\n")
cat("Output     :", OUT_FIGURAS, "\n\n")

# ── 1. Carregar sumstats GWAS ──────────────────────────────────────────────────
if (!file.exists(GWAS_ASSOC)) stop("GWAS assoc não encontrado: ", GWAS_ASSOC)
gwas_raw <- fread(GWAS_ASSOC)
setnames(gwas_raw, tolower(names(gwas_raw)))

required <- c("chr", "pos", "id", "p")
missing  <- setdiff(required, names(gwas_raw))
if (length(missing) > 0) stop("Colunas obrigatórias ausentes no assoc: ",
                              paste(missing, collapse = ", "))
gwas_raw[, p := suppressWarnings(as.numeric(p))]
gwas_raw <- gwas_raw[!is.na(p) & is.finite(p)]
cat(sprintf("SNPs no assoc: %d\n", nrow(gwas_raw)))

# Subset mínimo para o locuszoomr (chrom, pos, id, p)
gwas <- gwas_raw[, .(chr, pos, id, p)]

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

# ── 4. Exportar GWAS completo para LocalZoom (tabix) ───────────────────────────
# https://statgen.github.io/localzoom/  → arrastar .tsv.gz + .tbi na página.
cat("\n[4] Exportando GWAS para LocalZoom (bgzip + tabix)...\n")

# Colunas essenciais: chrom pos ref alt id a1_freq beta se p logp
lz_cols <- c("chr", "pos", "ref", "alt", "id", "a1_freq", "beta", "se", "p", "logp")
lz_cols <- intersect(lz_cols, names(gwas_raw))
if (length(lz_cols) < 4 || !all(c("chr", "pos", "p") %in% lz_cols)) {
  cat("  AVISO: colunas essenciais ausentes no assoc; LocalZoom exportado com as disponíveis.\n")
  lz_cols <- intersect(c("chr", "pos", "p", names(gwas_raw)), names(gwas_raw))
}
lz <- gwas_raw[, lz_cols, with = FALSE]

# Normalizar chr para números (tabix requer chr consistente, ex: "1" não "chr1")
lz[, chr_num := suppressWarnings(as.numeric(sub("^chr", "", chr)))]
lz[is.na(chr_num), chr_num := 23]          # X → 23 (aproximação para ordenação)
lz[, chr := as.character(chr_num)]         # gravar chr numérico no arquivo
setkey(lz, chr_num, pos)
lz[, chr_num := NULL]

# Header deve começar com "#" para o tabix ignorá-lo automaticamente
setnames(lz, c("chrom", "pos", "ref", "alt", "id", "a1_freq", "beta", "se", "p", "logp")[seq_along(lz_cols)])

lz_tsv <- file.path(OUT_TABELAS, paste0(LZ_BASE, ".tsv"))
writeLines(paste0("#", paste(names(lz), collapse = "\t")), lz_tsv)
fwrite(lz, lz_tsv, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
cat("  TSV:", lz_tsv, "\n")

# bgzip + tabix (o tabix está no PATH do env locuszoom)
tools_ok <- Sys.which(c("bgzip", "tabix"))
for (t in names(tools_ok)) {
  if (!nzchar(tools_ok[[t]])) {
    cat(sprintf("  AVISO: '%s' não encontrado no PATH do env locuszoom (Sys.which vazio).\n", t))
  } else {
    cat(sprintf("  %-6s -> %s\n", t, tools_ok[[t]]))
  }
}
lz_gz  <- paste0(lz_tsv, ".gz")
err_tmp <- tempfile(fileext = ".err")
status_bgzip <- system2("bgzip", args = c("-f", "-c", lz_tsv),
                        stdout = lz_gz, stderr = err_tmp)
if (!identical(as.integer(status_bgzip), 0L)) {
  msg <- if (file.exists(err_tmp)) paste(readLines(err_tmp, warn = FALSE), collapse = "\n") else ""
  stop("bgzip falhou (exit ", status_bgzip, "):\n", msg)
}
status_tabix <- system2("tabix",
                        args = c("-f", "-s", "1", "-b", "2", "-e", "2", lz_gz),
                        stderr = err_tmp)
if (!identical(as.integer(status_tabix), 0L)) {
  msg <- if (file.exists(err_tmp)) paste(readLines(err_tmp, warn = FALSE), collapse = "\n") else ""
  stop("tabix falhou (exit ", status_tabix, "):\n", msg)
}
file.remove(err_tmp)

file.remove(lz_tsv)   # mantém apenas o .tsv.gz + .tbi
cat("  LocalZoom: ", lz_gz, " + .tbi\n", sep = "")
cat("  Uso: abrir https://statgen.github.io/localzoom/ e arrastar os dois arquivos.\n")
cat("  Mapear colunas: chrom=1, pos=2, ref=3, alt=4, id=5, a1_freq=6, beta=7, se=8, p=9, logp=10\n")

cat("\nDone.\n")
