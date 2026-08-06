#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# 04.4_hits_gwas.R
# Filtra hits significativos e sugestivos do GWAS original e anota com VEP.
#
# Inputs:
#   args[1] → arquivo GWAS (.assoc): SRR_gwas.assoc
#   args[2] → VCF anotado com VEP (opcional): joint_germline.filtered.norm.AlphaMissense.vep.vcf.gz
#   args[3] → diretório de saída (default: resultados/tabelas)
#
# Outputs:
#   resultados/tabelas/gwas_hits_significativos.csv
#   resultados/tabelas/gwas_hits_sugestivos.csv
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))
options(width = 120)

# ── Args ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Uso: Rscript 04.4_hits_gwas.R <gwas_assoc> [vep_vcf] [out_dir]")
}
gwas_file <- args[1]
vep_vcf   <- if (length(args) >= 2 && file.exists(args[2])) args[2] else NULL
out_dir   <- if (length(args) >= 3) args[3] else "resultados/tabelas"

# ── Helper: caminho absoluto relativo ao repo ─────────────────────────────────
script_dir <- tryCatch({
  args0 <- commandArgs(trailingOnly = FALSE)
  fidx  <- grep("^--file=", args0)
  if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
}, error = function(e) getwd())

repo_root <- normalizePath(file.path(script_dir, ".."))
abs_path <- function(p) if (!is.null(p) && grepl("^/", p)) p else file.path(repo_root, p)

gwas_file <- abs_path(gwas_file)
out_dir   <- abs_path(out_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("== 04.4_hits_gwas.R ==\n")
cat("GWAS file:", gwas_file, "\n")
cat("VEP VCF :", if (is.null(vep_vcf)) "IGNORADO" else vep_vcf, "\n")
cat("Output  :", out_dir, "\n\n")

# ── 1. Carregar sumstats GWAS ──────────────────────────────────────────────────
cat("[1/4] Carregando sumstats GWAS...\n")
dt <- fread(gwas_file)
setnames(dt, tolower(names(dt)))  # normaliza nomes das colunas

# Mapear colunas alternativas (PLINK assoc usa 'or'/'log(or)_se' em vez de 'beta'/'se')
if (!"beta" %in% names(dt) && "or" %in% names(dt)) {
  setnames(dt, "or", "beta")
}
if (!"se" %in% names(dt) && "log(or)_se" %in% names(dt)) {
  setnames(dt, "log(or)_se", "se")
}

# Normalizar contigs para formato VCF (adicionar prefixo "chr" se ausente)
# O VCF do Sarek usa contigs com prefixo chr (ex: chr1), enquanto o .assoc pode usar apenas números
dt[, chr := ifelse(grepl("^chr|^GL|^KI|^JH|^A|un_", chr), chr, paste0("chr", chr))]

required_cols <- c("chr", "pos", "id", "ref", "alt", "beta", "p")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Colunas obrigatórias ausentes no GWAS: ", paste(missing_cols, collapse = ", "))
}

cat(sprintf("Total de SNPs carregados: %d\n", nrow(dt)))

# ── 2. Filtrar hits significativos e sugestivos ────────────────────────────────
cat("[2/4] Aplicando thresholds...\n")
dt[, is_sig := fifelse(p < 5e-8, TRUE, FALSE)]
dt[, is_sug := fifelse(p < 1e-5 & !is_sig, TRUE, FALSE)]

sig_hits <- dt[dt$is_sig]
sug_hits <- dt[dt$is_sug]

cat(sprintf("Hits significativos (P < 5e-8): %d\n", nrow(sig_hits)))
cat(sprintf("Hits sugestivos     (P < 1e-5): %d\n", nrow(sug_hits)))

# ── 3. Anotar com VEP (se disponível) ──────────────────────────────────────────
if (!is.null(vep_vcf)) {
  cat("[3/4] Anotando com VEP...\n")
  
  # Verificar se VCF e índice existem
  vep_tbi <- paste0(vep_vcf, ".tbi")
  if (!file.exists(vep_vcf)) {
    cat("  ERRO: VCF VEP não encontrado:", vep_vcf, "\n")
  } else if (!file.exists(vep_tbi)) {
    cat("  AVISO: Índice .tbi não encontrado:", vep_tbi, "\n")
    cat("  Criando índice...\n")
    system2("bcftools", args = c("index", vep_vcf), stdout = FALSE, stderr = FALSE)
  }

  # Extrair apenas regiões de interesse do VCF via bcftools
  regions_file <- tempfile(fileext = ".bed")
  all_hits <- rbind(sig_hits[, .(chr, pos, id)],
                    sug_hits[, .(chr, pos, id)],
                    fill = TRUE)
  # Converter para formato BED (0-based start, 1-based end)
  bed_dt <- data.table(
    chr   = all_hits$chr,
    start = all_hits$pos - 1,
    end   = all_hits$pos
  )
  fwrite(bed_dt, regions_file, sep = "\t", col.names = FALSE)
  cat("  Regiões BED geradas (primeiras 5):\n")
  head_lines <- readLines(regions_file, n = 5)
  for (hl in head_lines) cat("    ", hl, "\n", sep = "")

  vep_tmp <- tempfile(fileext = ".tsv")
  
  # Estratégia 1: bcftools query com -R (BED)
  bcftools_cmd <- sprintf(
    "bcftools query -R '%s' -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t%%CSQ\\n' '%s' > '%s' 2>&1",
    regions_file, vep_vcf, vep_tmp
  )
  cat("  Comando VEP (query -R): ", bcftools_cmd, "\n", sep = "")
  system(bcftools_cmd, intern = FALSE)
  
  vep_dt <- tryCatch(fread(vep_tmp, sep = "\t",
                           col.names = c("chr", "pos", "id", "ref", "alt", "csq")),
                     error = function(e) NULL)

  # Estratégia 2: fallback com bcftools view + grep (mais lento mas mais confiável)
  if (is.null(vep_dt) || nrow(vep_dt) == 0) {
    cat("  Fallback: usando bcftools view + awk...\n")
    vep_tmp2 <- tempfile(fileext = ".tsv")
    
    # Concatenar posições para -r
    regions_str <- paste(sprintf("%s:%d-%d", bed_dt$chr, bed_dt$start + 1, bed_dt$end), collapse = ",")
    bcftools_cmd2 <- sprintf(
      "bcftools view -r '%s' '%s' 2>/dev/null | bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t%%CSQ\\n' > '%s' 2>&1",
      regions_str, vep_vcf, vep_tmp2
    )
    cat("  Comando VEP (view -r): ", bcftools_cmd2, "\n", sep = "")
    system(bcftools_cmd2, intern = FALSE)
    
    vep_dt <- tryCatch(fread(vep_tmp2, sep = "\t",
                             col.names = c("chr", "pos", "id", "ref", "alt", "csq")),
                       error = function(e) NULL)
    
    if (!is.null(vep_dt)) unlink(vep_tmp2)
  }

  if (!is.null(vep_dt) && nrow(vep_dt) > 0) {
    # Extrair consequence, symbol, gene, impact do CSQ (primeiro transcript)
    parse_csq <- function(csq_str) {
      if (is.na(csq_str) || csq_str == "") return(list("", "", "", ""))
      parts <- unlist(strsplit(csq_str, ",", fixed = TRUE))
      first <- unlist(strsplit(parts[1], "|", fixed = TRUE))
      result <- list(
        consequence = if (length(first) >= 2) first[2] else "",
        symbol      = if (length(first) >= 4) first[4] else "",
        gene        = if (length(first) >= 5) first[5] else "",
        impact      = if (length(first) >= 3) first[3] else ""
      )
      return(result)
    }

    # Parsear CSQ para cada linha
    parsed_list <- lapply(vep_dt$csq, parse_csq)
    vep_dt[, consequence := sapply(parsed_list, `[[`, "consequence")]
    vep_dt[, symbol       := sapply(parsed_list, `[[`, "symbol")]
    vep_dt[, gene         := sapply(parsed_list, `[[`, "gene")]
    vep_dt[, impact       := sapply(parsed_list, `[[`, "impact")]

    vep_clean <- unique(vep_dt[, .(chr, pos, id, consequence, symbol, gene, impact)])
    dt <- merge(dt, vep_clean, by = c("chr", "pos", "id"), all.x = TRUE)
    sig_hits <- merge(sig_hits, vep_clean, by = c("chr", "pos", "id"), all.x = TRUE)
    sug_hits <- merge(sug_hits, vep_clean, by = c("chr", "pos", "id"), all.x = TRUE)

    unlink(vep_tmp); unlink(regions_file)
    if (exists("vep_tmp2")) unlink(vep_tmp2)
  } else {
    cat("⚠ AVISO: Ainda sem anotação VEP — verifique no cluster:\n")
    cat("  - VCF existe? ", file.exists(vep_vcf), "\n")
    cat("  - Índice .tbi existe? ", file.exists(vep_tbi), "\n")
    cat("  - Contig naming (chr vs number): OK (normalizado)\n")
    cat("  - bcftools query -R funcionando?\n")
    cat(sprintf("    bcftools query -R <%s> -f '%%CHROM\\t%%POS\\t%%CSQ\\n' '%s'\n",
                regions_file, vep_vcf))
  }
}

# ── 4. Exportar CSVs ───────────────────────────────────────────────────────────
cat("[4/4] Salvando resultados...\n")

fwrite(sig_hits, file.path(out_dir, "gwas_hits_significativos.csv"))
fwrite(sug_hits, file.path(out_dir, "gwas_hits_sugestivos.csv"))

cat(sprintf("\n✔ Salvo: %s/gwas_hits_significativos.csv (%d hits)\n", out_dir, nrow(sig_hits)))
cat(sprintf("✔ Salvo: %s/gwas_hits_sugestivos.csv (%d hits)\n", out_dir, nrow(sug_hits)))
cat("\nDone.\n")
