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
PCUTOFF     <- 5e-8                    # limiar de significância (genome-wide)
SUG_CUTOFF  <- 1e-5                    # limiar de sugestivo
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

# Header padrão (sem "#"), nomes limpos para a auto-detecção do LocalZoom.
# O tabix pula a 1ª linha via --skip-lines 1 (recomendado pela doc do LocalZoom).
# lz já herda os nomes corretos de lz_cols; basta renomear chr -> chrom.
setnames(lz, "chr", "chrom")

lz_tsv <- file.path(OUT_TABELAS, paste0(LZ_BASE, ".tsv"))
writeLines(paste(names(lz), collapse = "\t"), lz_tsv)
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
                        args = c("-f", "-s", "1", "-b", "2", "-e", "2",
                                 "-S", "1", lz_gz),
                        stderr = err_tmp)
if (!identical(as.integer(status_tabix), 0L)) {
  msg <- if (file.exists(err_tmp)) paste(readLines(err_tmp, warn = FALSE), collapse = "\n") else ""
  stop("tabix falhou (exit ", status_tabix, "):\n", msg)
}
file.remove(err_tmp)

file.remove(lz_tsv)   # mantém apenas o .tsv.gz + .tbi
cat("\n=== Export para LocalZoom (https://statgen.github.io/localzoom/) ===\n")
cat("Arquivos: ", lz_gz, " + .tbi\n", sep = "")
cat("Uso: abrir a página, arrastar os DOIS arquivos e mapear as colunas abaixo.\n")
lz_fields <- c(
  chrom   = "cromossomo",
  pos     = "posição (bp)",
  ref     = "alelo de referência",
  alt     = "alelo alternativo (effect allele)",
  id      = "ID do SNP",
  beta    = "efeito (beta)",
  a1_freq = "frequência do effect allele",
  se      = "erro padrão (se)",
  p       = "p-value (bruto)",
  logp    = "p-value em -log10"
)
cat("  Mapeamento das colunas do arquivo:\n")
for (nm in names(lz)) {
  desc <- if (nm %in% names(lz_fields)) lz_fields[[nm]] else "(coluna extra)"
  cat(sprintf("    coluna %-2d  %-8s = %s\n", which(names(lz) == nm), nm, desc))
}
cat("  No 'p-value' do mapeamento: use a coluna 'p' e responda 'Não' a '-log10',\n")
cat("  ou use a coluna 'logp' e responda 'Sim'.\n")

# ── 5. Manhattan dos hits com genes (significativos e sugestivos) ──────────────
# Dois plots, cada um apenas com as variantes do respectivo threshold, rotuladas
# com os genes anotados pelo VEP (04.4). ggplot2/ggrepel já são dependências do
# locuszoomr, então estão disponíveis no env locuszoom.
cat("\n[5] Manhattan plots (hits + genes acima do threshold)...\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
  if (!have_ggrepel) cat("  AVISO: ggrepel indisponível — rótulos sem repulsão.\n")

  # ---- posição genômica acumulada (chr 1-22, com gap) ----
  manh <- gwas_raw[, .(chr, pos, id, p)]
  manh[, chr_num := suppressWarnings(as.numeric(sub("^chr", "", chr)))]
  manh[, pos := as.numeric(pos)]
  manh[, logp := -log10(p)]
  manh <- manh[!is.na(chr_num) & !is.na(pos) & is.finite(logp)]

  chr_max <- manh[, .(maxp = max(pos)), by = chr_num][order(chr_num)]
  gap <- round(0.02 * max(chr_max$maxp))
  chr_max[, off := c(0, cumsum(maxp + gap)[- .N])]
  manh <- merge(manh, chr_max[, .(chr_num, off)], by = "chr_num")
  manh[, x := off + pos]
  chr_center <- chr_max[, .(x = off + maxp / 2)]

  # ---- genes por (chr,pos) a partir dos CSVs de hits (04.4/VEP) ----
  # gene canônico: o CSQ do VEP/AlphaMissense do Sarek não grava CANONICAL,
  # então o transcrito canônico (Ensembl) é obtido de EnsDb v86 por overlap.
  build_canon_ranges <- function(ens_db_pkg) {
    if (!requireNamespace("ensembldb", quietly = TRUE) ||
        !requireNamespace(ens_db_pkg, quietly = TRUE)) {
      cat("  AVISO: ensembldb/", ens_db_pkg, " indisponível — sem gene canônico.\n", sep = "")
      return(NULL)
    }
    edb <- get(ens_db_pkg, envir = asNamespace(ens_db_pkg))
    canon <- tryCatch(as.data.frame(ensembldb::canonicalTranscript(edb)),
                      error = function(e) NULL)
    if (is.null(canon) || nrow(canon) == 0 ||
        !all(c("tx_id", "gene_id") %in% names(canon))) {
      cat("  AVISO: canonicalTranscript() falhou — sem gene canônico.\n")
      return(NULL)
    }
    tx <- tryCatch(ensembldb::transcripts(
            edb,
            filter = AnnotationFilter::GeneIdFilter(unique(as.character(canon$gene_id))),
            columns = c("tx_id", "gene_id", "seq_name", "tx_start", "tx_end")),
          error = function(e) NULL)
    if (is.null(tx) || length(tx) == 0) {
      cat("  AVISO: ranges dos transcritos indisponíveis — sem gene canônico.\n")
      return(NULL)
    }
    tx_dt <- data.table(
      tx_id   = as.character(tx$tx_id),
      gene_id = as.character(tx$gene_id),
      chr_num = suppressWarnings(as.numeric(sub("^chr", "", as.character(tx$seq_name)))),
      start   = as.numeric(tx$tx_start),
      end     = as.numeric(tx$tx_end))
    canon <- data.table(tx_id = as.character(canon$tx_id),
                        gene_id = as.character(canon$gene_id))
    cr <- tx_dt[canon, on = .(tx_id, gene_id), nomatch = NULL]
    cr <- cr[!is.na(chr_num)]
    cat(sprintf("  Transcritos canônicos (Ensembl v86): %d (genes: %d)\n",
                nrow(cr), length(unique(cr$gene_id))))
    cr[, .(gene_id, chr_num, start, end)]
  }
  canon_ranges <- build_canon_ranges(ENS_DB)

  genes_from_hits <- function(hits_file) {
    if (!file.exists(hits_file)) return(NULL)
    h <- fread(hits_file)
    setnames(h, tolower(names(h)))
    if (!all(c("chr", "pos", "symbol") %in% names(h))) return(NULL)
    h[, chr_num := suppressWarnings(as.numeric(sub("^chr", "", chr)))]
    h[, pos := as.numeric(pos)]
    h <- h[!is.na(chr_num) & !is.na(symbol) & symbol != "" & symbol != "."]

    # fast-path: se o CSV já tiver a flag canônica do VEP
    if ("canonical" %in% names(h)) {
      canon <- h[tolower(trimws(canonical)) == "yes"]
      if (nrow(canon) > 0) {
        cat("  gene canônico: coluna 'canonical' do CSV usada.\n")
        h <- canon
      } else {
        cat("  AVISO: coluna 'canonical' sem YES — tentando EnsDb por overlap.\n")
      }
    }

    n_loci <- uniqueN(h[, .(chr_num, pos)])

    # gene canônico por overlap: mantém o gene cujo transcrito canônico (EnsDb v86)
    # cobre a posição da variante. Loci sem cobertura canônica → fallback (todos).
    if (!is.null(canon_ranges) && "gene" %in% names(h)) {
      ug <- unique(h[!is.na(gene) & gene != "" & gene != ".", .(chr_num, pos, gene)])
      if (nrow(ug) > 0) {
        mg <- merge(ug, canon_ranges, by.x = "gene", by.y = "gene_id",
                    allow.cartesian = TRUE, sort = FALSE)
        mg <- mg[chr_num.x == chr_num.y & start <= pos & pos <= end]
        if (nrow(mg) > 0) {
          canon_genes <- unique(mg[, .(chr_num = chr_num.x, pos, gene)])
          syms <- unique(merge(canon_genes,
                               unique(h[, .(chr_num, pos, gene, symbol)]),
                               by = c("chr_num", "pos", "gene")))
          lab <- syms[, .(genes = paste(unique(symbol), collapse = "; ")),
                      by = .(chr_num, pos)]
          cat(sprintf("  gene canônico por overlap (EnsDb v86): %d/%d loci com canônico.\n",
                      nrow(lab), n_loci))
          if (nrow(lab) < n_loci) {
            cat(sprintf("  AVISO: %d loci sem transcrito canônico cobrindo a posição (fallback).\n",
                        n_loci - nrow(lab)))
          }
          all_lab <- h[, .(genes_all = paste(unique(symbol), collapse = "; ")),
                       by = .(chr_num, pos)]
          out <- merge(all_lab, lab, by = c("chr_num", "pos"), all.x = TRUE)
          out[, genes := ifelse(is.na(genes), genes_all, genes)]
          return(out[, .(chr_num, pos, genes)])
        }
      }
    }

    if (is.null(canon_ranges)) {
      cat("  AVISO: sem mapa canônico — usando todos os SYMBOL.\n")
    } else if (!"gene" %in% names(h)) {
      cat("  AVISO: sem coluna 'gene' no hits CSV — usando todos os SYMBOL.\n")
    } else {
      cat("  AVISO: nenhum transcrito canônico cobrindo os loci — usando todos os SYMBOL.\n")
    }
    h[, .(genes = paste(unique(symbol), collapse = "; ")), by = .(chr_num, pos)]
  }
  sig_genes <- genes_from_hits(HITS_SIG)
  sug_genes <- genes_from_hits(HITS_SUG)

  plot_manh <- function(dt, genes, out_png, title, thr_lines) {
    if (is.null(genes)) {
      dt[, genes := NA_character_]
    } else {
      dt <- merge(dt, genes, by = c("chr_num", "pos"), all.x = TRUE)
    }
    dt <- dt[order(x)]
    p <- ggplot(dt, aes(x = x, y = logp)) +
      geom_point(shape = 21, colour = "grey30", fill = "#D95F02",
                 size = 2.2, alpha = 0.9) +
      scale_x_continuous(breaks = chr_center$x, labels = chr_max$chr_num,
                         expand = expansion(mult = 0.01)) +
      xlab("Cromossomo") + ylab(expression(-log[10](p))) +
      ggtitle(title) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    for (i in seq_len(nrow(thr_lines))) {
      p <- p + geom_hline(yintercept = thr_lines$y[i], linetype = "dashed",
                          colour = thr_lines$col[i])
    }
    lbl <- dt[!is.na(genes)]
    if (nrow(lbl) > 0) {
      if (have_ggrepel) {
        p <- p + ggrepel::geom_text_repel(data = lbl, aes(label = genes),
                                          size = 3.2, max.overlaps = 30,
                                          min.segment.length = 0,
                                          segment.colour = "grey40", seed = 42)
      } else {
        p <- p + geom_text(data = lbl, aes(label = genes), size = 3.2,
                           nudge_y = 0.25, vjust = 0, angle = 30)
      }
    }
    ggsave(out_png, p, width = 12, height = 7, dpi = 300)
    cat("  ", basename(out_png), " (", nrow(dt), " SNPs, ", nrow(lbl), " rotulados)\n",
        sep = "")
  }

  # thresholds: genômica (5e-8) e sugestiva (1e-5), destacando a relevante
  thr_sig <- data.frame(y = c(-log10(PCUTOFF), -log10(SUG_CUTOFF)),
                        col = c("firebrick", "grey50"))
  thr_sug <- data.frame(y = c(-log10(PCUTOFF), -log10(SUG_CUTOFF)),
                        col = c("grey60", "firebrick"))

  # apenas hits significativos (P < 5e-8)
  sig_dt <- manh[p < PCUTOFF]
  if (nrow(sig_dt) > 0) {
    cat("  Manhattan significativos (P < 5e-8):", nrow(sig_dt), "SNPs\n")
    plot_manh(sig_dt, sig_genes,
              file.path(OUT_FIGURAS, "manhattan_hits_significativos.png"),
              "GWAS — Hits significativos (P < 5e-8)", thr_sig)
  } else {
    cat("  Sem hits significativos para plotar.\n")
  }

  # apenas sugestivos (5e-8 <= P < 1e-5)
  sug_dt <- manh[p >= PCUTOFF & p < SUG_CUTOFF]
  if (nrow(sug_dt) > 0) {
    cat("  Manhattan sugestivos (5e-8 <= P < 1e-5):", nrow(sug_dt), "SNPs\n")
    plot_manh(sug_dt, sug_genes,
              file.path(OUT_FIGURAS, "manhattan_hits_sugestivos.png"),
              "GWAS — Hits sugestivos (P < 1e-5)", thr_sug)
  } else {
    cat("  Sem hits sugestivos para plotar.\n")
  }
} else {
  cat("  AVISO: ggplot2 indisponível no env — Manhattan pulado.\n")
}

cat("\nDone.\n")
