#!/usr/bin/env Rscript
# 05_locus_plots.R
# Gera locus plots (manhattan regional) dos top loci fine-mapped do MVP,
# sobrepostos às summary statistics do exoma.
#
# Input :
#   - resultados/tabelas/lookup_completo.csv   (output do 02_lookup_exoma.R)
#   - dados/exoma/gwas_dbgap_+vep.tsv          (sumstats exoma, para o locus)
# Output: resultados/figuras/locus_plot_*.png

# ─────────────────────────── CONFIG ───────────────────────────
LOOKUP_CSV <- "resultados/tabelas/lookup_completo.csv"
EXOMA_SUMSTATS <- "dados/exoma/gwas_dbgap_sem_filtragem+joint_call+vep.tsv"
OUT_FIGURAS <- "resultados/figuras"

# Nº de top loci a plotar
N_TOP_LOCI <- 5
# Janela (bp) ao redor do SNP index para o locus plot
LODUS_WINDOW_BP <- 500000

# Colunas do lookup
COL_RSID <- "RSID"
COL_PIP  <- "Overall PIP"
COL_LOCUS <- "Locus"     # ex.: "chr1:1000000-2000000" (se disponível)
COL_CHR  <- "CHR"        # do exoma
COL_BP   <- "BP"         # do exoma
COL_PVAL <- "PVAL_exoma"

# Colunas do exoma (sumstats)
EX_RSID <- "rsid"; EX_CHR <- "chr"; EX_BP <- "bp"; EX_PVAL <- "pval"
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
if (!dir.exists(OUT_FIGURAS)) dir.create(OUT_FIGURAS, recursive = TRUE)

lookup <- read.csv(LOOKUP_CSV, stringsAsFactors = FALSE, check.names = FALSE)
cat("== 05_locus_plots.R ==\n")

exoma <- read.delim(EXOMA_SUMSTATS, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "#")
rename_map <- c(EX_RSID = "RSID", EX_CHR = "CHR", EX_BP = "BP", EX_PVAL = "PVAL")
exoma <- exoma %>% rename(any_of(rename_map))
exoma$PVAL <- suppressWarnings(as.numeric(exoma$PVAL))
exoma <- exoma %>% filter(!is.na(PVAL), PVAL > 0) %>% mutate(mlog10p = -log10(PVAL))

# Selecionar top loci por PIP (genotipados ou não — prioriza PIP alto)
top_loci <- lookup %>%
  filter(suppressWarnings(as.numeric(.data[[COL_PIP]])) > 0.8) %>%
  arrange(desc(suppressWarnings(as.numeric(.data[[COL_PIP]])))) %>%
  head(N_TOP_LOCI)

if (nrow(top_loci) == 0) {
  message("Nenhum locus com PIP > 0.8 encontrado. Usando top por PIP geral.")
  top_loci <- lookup %>%
    arrange(desc(suppressWarnings(as.numeric(.data[[COL_PIP]])))) %>%
    head(N_TOP_LOCI)
}

cat(sprintf("Gerando locus plots para %d loci...\n", nrow(top_loci)))

for (i in seq_len(nrow(top_loci))) {
  rsid <- top_loci[[COL_RSID]][i]

  # Definir região: usar COL_CHR/COL_BP do próprio lookup (se houver) ou Locus
  chr_i <- if (COL_CHR %in% colnames(top_loci)) top_loci[[COL_CHR]][i] else NA
  bp_i  <- if (COL_BP  %in% colnames(top_loci)) top_loci[[COL_BP]][i] else NA

  # Fallback: extrair de COL_LOCUS "chr1:1000000-2000000"
  if ((is.na(chr_i) || is.na(bp_i)) && COL_LOCUS %in% colnames(top_loci)) {
    m <- regmatches(top_loci[[COL_LOCUS]][i],
                    regexec("chr?([0-9XY]+):([0-9]+)-([0-9]+)",
                            top_loci[[COL_LOCUS]][i]))[[1]]
    if (length(m) == 4) {
      chr_i <- m[2]
      # usar o meio do intervalo
      bp_i <- (as.numeric(m[3]) + as.numeric(m[4])) / 2
    }
  }
  if (is.na(chr_i) || is.na(bp_i)) {
    warning(sprintf("Locus de %s não identificado — pulando.", rsid))
    next
  }

  chr_num <- sub("^chr", "", chr_i)
  start <- as.numeric(bp_i) - LODUS_WINDOW_BP
  end   <- as.numeric(bp_i) + LODUS_WINDOW_BP

  region <- exoma %>%
    filter(CHR == chr_num, BP >= start, BP <= end)

  if (nrow(region) == 0) {
    message(sprintf("Nenhum SNP do exoma na região %s:%s-%s (rsid=%s) — pulando.",
                    chr_i, start, end, rsid))
    next
  }

  gene_name <- "GENE"  # substituir pelo gene real (ex.: ATM, BRCA2)

  p <- ggplot(region, aes(x = BP, y = mlog10p)) +
    geom_point(alpha = 0.6, color = "#2c7bb6", size = 1.2) +
    geom_point(data = region %>% filter(RSID == rsid),
               aes(x = BP, y = mlog10p), color = "#d7191c", size = 3) +
    labs(title = paste0(gene_name, " (", rsid, ")"),
         subtitle = sprintf("Região: chr%s:%s-%s", chr_i, start, end),
         x = sprintf("Posição (chr%s)", chr_i),
         y = expression(-log[10](p))) +
    theme_minimal()

  ggsave(file.path(OUT_FIGURAS, paste0("locus_plot_", gene_name, ".png")),
         p, width = 8, height = 5, dpi = 300)
  cat(sprintf("  Salvo: locus_plot_%s.png\n", gene_name))
}

cat("OK. Locus plots em ", OUT_FIGURAS, "\n", sep = "")
