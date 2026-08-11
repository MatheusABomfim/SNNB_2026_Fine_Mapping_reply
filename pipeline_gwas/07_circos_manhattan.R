#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# 07_circos_manhattan.R
# Gera circos plot e Manhattan em arquivos separados, com genes anotados
# (VEP, do 04.4) apenas nos hits genome-wide (P < 5e-8). O gene por locus é o
# "gene canônico" derivado do próprio VEP: transcrito de maior IMPACT
# (HIGH > MODERATE > LOW > MODIFIER), desempatado pelo transcript_index.
#
# Inputs:
#   - GWAS assoc completo : /storage4/.../gwas/gwas_prod/SRR_gwas.assoc
#   - Hits significativos : resultados/tabelas/gwas_hits_significativos.csv
#   - Hits sugestivos     : resultados/tabelas/gwas_hits_sugestivos.csv
# Output:
#   - resultados/figuras/circos_plot.pdf/.png            (circos, separado)
#   - resultados/figuras/fig_circos_manhattan.pdf/.png   (só Manhattan)
#
# Dependências R (env locuszoom, via micromamba):
#   micromamba install -n locuszoom -c conda-forge r-circlize
#   # ggplot2, ggrepel, data.table já presentes via locuszoomr
# ──────────────────────────────────────────────────────────────────────────────

# ─────────────────────────── CONFIG ───────────────────────────
GWAS_ASSOC  <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/gwas_prod/SRR_gwas.assoc"
HITS_SIG    <- "resultados/tabelas/gwas_hits_significativos.csv"
HITS_SUG    <- "resultados/tabelas/gwas_hits_sugestivos.csv"
OUT_FIGURAS <- "resultados/figuras"

PCUTOFF     <- 5e-8                    # limiar genome-wide
SUG_CUTOFF  <- 1e-5                    # limiar sugestivo
# ───────────────────────────────────────────────────────────────

# Caminhos relativos à raiz do repo (o pai do dir deste script).
abs_path <- function(p) if (grepl("^/", p)) p else file.path(REPO_ROOT, p)
args0 <- commandArgs(trailingOnly = FALSE)
fidx <- grep("^--file=", args0)
SCRIPT_DIR <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))
HITS_SIG    <- abs_path(HITS_SIG)
HITS_SUG    <- abs_path(HITS_SUG)
OUT_FIGURAS <- abs_path(OUT_FIGURAS)

suppressPackageStartupMessages(library(data.table))

if (!dir.exists(OUT_FIGURAS)) dir.create(OUT_FIGURAS, recursive = TRUE)
cat("== 07_circos_manhattan.R ==\n")
cat("GWAS assoc :", GWAS_ASSOC, "\n")
cat("Thresholds : P <", PCUTOFF, "(genome-wide) | P <", SUG_CUTOFF, "(sugestivo)\n")
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

manh <- gwas_raw[, .(chr, pos, id, p)]
manh[, chr_num := suppressWarnings(as.numeric(sub("^chr", "", chr)))]
manh[, pos := as.numeric(pos)]
manh[, logp := -log10(p)]
manh <- manh[!is.na(chr_num) & !is.na(pos) & is.finite(logp)]

# ── 2. Genes dos hits (VEP) por (chr,pos) ──────────────────────────────────────
# "Gene canônico" derivado só do VEP (CSV do 04.4): por locus, o gene do
# transcrito de maior IMPACT (HIGH > MODERATE > LOW > MODIFIER), desempate pelo
# transcript_index menor (VEP lista a consequência mais grave primeiro).
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
    }
  }

  # rank do IMPACT (menor = mais grave)
  if ("impact" %in% names(h)) {
    h[, impact_rank := fifelse(tolower(trimws(impact)) == "high", 1,
                       fifelse(tolower(trimws(impact)) == "moderate", 2,
                       fifelse(tolower(trimws(impact)) == "low", 3,
                       fifelse(tolower(trimws(impact)) == "modifier", 4, 99))))]
  } else {
    h[, impact_rank := 99L]
    cat("  AVISO: sem coluna 'impact' — usando só transcript_index.\n")
  }

  # desempate por transcript_index (ordem do VEP); faltando, mantém a ordem
  if (!"transcript_index" %in% names(h)) {
    h[, transcript_index := seq_len(.N), by = .(chr_num, pos)]
    cat("  AVISO: sem coluna 'transcript_index' — usando a ordem das linhas.\n")
  }
  h[, transcript_index := as.numeric(transcript_index)]

  n_loci <- uniqueN(h[, .(chr_num, pos)])
  setorder(h, chr_num, pos, impact_rank, transcript_index)
  lab <- h[, .(genes = paste(unique(symbol), collapse = "; ")),
           by = .(chr_num, pos)]
  n1 <- sum(grepl(";", lab$genes, fixed = TRUE))
  cat(sprintf("  gene canônico por IMPACT (HIGH>MODERATE>LOW>MODIFIER): %d loci; %d com 1 gene.\n",
              nrow(lab), nrow(lab) - n1))
  if (n1 > 0) cat(sprintf("  AVISO: %d loci com >1 gene no topo (mesmo IMPACT/transcript_index).\n", n1))
  lab
}
sig_genes <- genes_from_hits(HITS_SIG)
if (!is.null(sig_genes)) {
  manh <- merge(manh, sig_genes, by = c("chr_num", "pos"), all.x = TRUE)
} else {
  manh[, genes := NA_character_]
  cat("AVISO: sem genes de hits significativos (HITS_SIG ausente/vazio).\n")
}

# ── 3. Coordenadas genômicas acumuladas (chr 1-22) ─────────────────────────────
chr_max <- manh[, .(maxp = max(pos)), by = chr_num][order(chr_num)]
gap <- round(0.02 * max(chr_max$maxp))
chr_max[, off := c(0, cumsum(maxp + gap)[- .N])]
manh <- merge(manh, chr_max[, .(chr_num, off)], by = "chr_num")
manh[, x := off + pos]
chr_center <- chr_max[, .(x = off + maxp / 2)]

# ── 4. Manhattan completo à direita (ggplot2) ──────────────────────────────────
cat("[4] Manhattan completo (direita)...\n")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  manh[, level := fifelse(p < PCUTOFF, "significativo",
                   fifelse(p < SUG_CUTOFF, "sugestivo",
                     fifelse(chr_num %% 2 == 1, "chr_impar", "chr_par")))]

  p_manh <- ggplot(manh, aes(x = x, y = logp, colour = level)) +
    geom_point(size = 1.3, alpha = 0.8) +
    scale_colour_manual(values = c(
        "significativo" = "#E31A1C",
        "sugestivo"     = "#FF7F00",
        "chr_impar"     = "#4477AA",
        "chr_par"       = "#BBCCEE"),
      breaks = c("significativo", "sugestivo"), name = NULL) +
    scale_x_continuous(breaks = chr_center$x, labels = chr_max$chr_num,
                       expand = expansion(mult = 0.01)) +
    geom_hline(yintercept = -log10(PCUTOFF), linetype = "dashed",
               colour = "#E31A1C", linewidth = 0.5) +
    geom_hline(yintercept = -log10(SUG_CUTOFF), linetype = "dashed",
               colour = "#FF7F00", linewidth = 0.5) +
    annotate("text", x = chr_center$x[1], y = -log10(PCUTOFF),
             label = "P = 5e-8", vjust = -0.5, hjust = 0, size = 3,
             colour = "#E31A1C") +
    annotate("text", x = chr_center$x[1], y = -log10(SUG_CUTOFF),
             label = "P = 1e-5", vjust = -0.5, hjust = 0, size = 3,
             colour = "#FF7F00") +
    xlab("Cromossomo") + ylab(expression(-log[10](p))) +
    ggtitle("GWAS Manhattan — CA Mama (todos os SNPs)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top")

  lbl <- manh[p < PCUTOFF & !is.na(genes)]
  if (nrow(lbl) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p_manh <- p_manh +
      ggrepel::geom_text_repel(data = lbl, aes(label = genes), size = 3,
                               max.overlaps = 30, min.segment.length = 0,
                               segment.colour = "grey40", seed = 42,
                               colour = "black")
  }
} else {
  cat("  AVISO: ggplot2 indisponível — sem painel Manhattan.\n")
  p_manh <- NULL
}

# ── 5. Circos (separado) + Manhattan (separado) ────────────────────────────────
cat("[5] Circos plot...\n")
have_circos <- requireNamespace("circlize", quietly = TRUE)

draw_circos <- function(manh, chr_max, sig_hits) {
  suppressPackageStartupMessages(library(circlize))
  chrlens <- chr_max[, .(chr_num, maxp)]
  maxlogp <- max(manh$logp)

  circos.par(gap.degree = 1, start.degree = 90, clock.wise = FALSE,
             track.margin = c(0.01, 0.01), cell.padding = c(0.02, 1, 0.02, 1))
  circos.initialize(factors = chrlens$chr_num,
                    xlim = cbind(rep(1, nrow(chrlens)), chrlens$maxp))

  # Track 1: ideograma com rótulos
  circos.track(ylim = c(0, 1), track.height = 0.08, bg.border = NA,
               panel.fun = function(x, y) {
                 i <- get.cell.meta.data("sector.numeric.index")
                 col <- if (i %% 2 == 1) "#2C3E50" else "#95A5A6"
                 circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
                             CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
                             col = col, border = NA)
                 circos.text(CELL_META$xcenter, CELL_META$ycenter,
                             get.cell.meta.data("sector.index"),
                             cex = 0.6, font = 2, col = "white",
                             facing = "downward")
               })

  # Track 2: -log10(p) de todos os SNPs
  circos.track(ylim = c(0, maxlogp * 1.05), track.height = 0.45,
               bg.border = NA, bg.col = "grey97",
               panel.fun = function(x, y) {
                 sid <- as.numeric(get.cell.meta.data("sector.index"))
                 df <- manh[chr_num == sid][order(pos)]
                 if (nrow(df) == 0) return(NULL)
                 col <- fifelse(df$p < PCUTOFF, "#E31A1C",
                         fifelse(df$p < SUG_CUTOFF, "#FF7F00",
                           fifelse(df$chr_num %% 2 == 1, "#4477AA", "#BBCCEE")))
                 circos.points(df$pos, df$logp, pch = 16, cex = 0.3, col = col)
                 circos.segments(CELL_META$cell.xlim[1], -log10(PCUTOFF),
                                 CELL_META$cell.xlim[2], -log10(PCUTOFF),
                                 col = "#E31A1C", lty = 2, lwd = 0.6)
               })

  # Track 3: hits genome-wide + genes
  circos.track(ylim = c(0, 1), track.height = 0.12, bg.border = NA,
               panel.fun = function(x, y) {
                 sid <- as.numeric(get.cell.meta.data("sector.index"))
                 df <- sig_hits[chr_num == sid]
                 if (nrow(df) == 0) return(NULL)
                 circos.points(df$pos, rep(0.4, nrow(df)), pch = 17,
                               cex = 0.9, col = "#E31A1C")
                 circos.text(df$pos, rep(1.05, nrow(df)), df$genes,
                             cex = 0.55, facing = "clockwise",
                             adj = c(0, 0), col = "#8B0000")
               })
}

sig_hits <- manh[p < PCUTOFF & !is.na(genes)]

if (have_circos) {
  # Circos sempre renderizado com device explícito (evita Rplots.pdf).
  circos_ok <- tryCatch({
    pdf(file.path(OUT_FIGURAS, "circos_plot.pdf"), width = 7, height = 7)
    draw_circos(manh, chr_max, sig_hits)
    dev.off()
    circos.clear()
    png(file.path(OUT_FIGURAS, "circos_plot.png"),
        width = 2100, height = 2100, res = 300)
    draw_circos(manh, chr_max, sig_hits)
    dev.off()
    circos.clear()
    cat("  circos_plot.pdf + circos_plot.png salvos.\n")
    TRUE
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    try(circos.clear(), silent = TRUE)
    cat("  AVISO: circos falhou (", conditionMessage(e), ").\n", sep = "")
    FALSE
  })
} else {
  cat("  AVISO: circlize indisponível — circos não gerado.\n")
}

# Manhattan sozinho (fig_circos_manhattan.*)
if (!is.null(p_manh)) {
  ggsave(file.path(OUT_FIGURAS, "fig_circos_manhattan.pdf"), p_manh,
         width = 10, height = 6)
  ggsave(file.path(OUT_FIGURAS, "fig_circos_manhattan.png"), p_manh,
         width = 10, height = 6, dpi = 300)
  cat("  fig_circos_manhattan.pdf + .png (Manhattan) salvos.\n")
} else {
  cat("  AVISO: sem Manhattan (ggplot2 indisponível).\n")
}

# limpeza: nunca deixar Rplots.pdf (default device) no working dir
if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

cat("\nDone.\n")
