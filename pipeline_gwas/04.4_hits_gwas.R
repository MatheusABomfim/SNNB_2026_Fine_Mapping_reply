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

   # Extrair definição dos campos CSQ do header do VCF
   csq_fields <- tryCatch({
     header_lines <- system2("bcftools", args = c("view", "-h", vep_vcf), stdout = TRUE, stderr = TRUE)
     csq_line <- grep("##INFO=<ID=CSQ", header_lines, value = TRUE)
     if (length(csq_line) > 0) {
       # Extract Format: Allele|Consequence|...
        fmt_match <- regmatches(csq_line, regexec("Format: ([^\\\"]*)", csq_line))[[1]][2]
       if (length(fmt_match) > 0) {
         fields <- strsplit(sub("Format:", "", fmt_match), "\\|")[[1]]
         fields <- trimws(fields)
         cat("  Campos CSQ detectados:", length(fields), "\n")
         cat("  ", paste(fields, collapse = ", "), "\n")
         fields
       } else {
         stop("Format field not found in CSQ header")
       }
     } else {
       stop("CSQ header not found in VCF")
     }
   }, error = function(e) {
     cat("  AVISO: Não foi possível ler header CSQ do VCF:", conditionMessage(e), "\n")
     cat("  Usando lista padrão de campos.\n")
     c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
       "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
       "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
       "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
       "SYMBOL_SOURCE", "HGNC_ID", "am_class", "am_genome",
       "am_pathogenicity", "am_protein_variant", "am_transcript_id",
       "am_uniprot_id", "am_gene_name")
   })


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
    # Parsear CSQ expandindo todos os campos e todos os transcripts
    parse_csq_all <- function(csq_str, fields) {
      if (is.na(csq_str) || csq_str == "" || csq_str == ".") return(NULL)
      transcripts <- unlist(strsplit(csq_str, ",", fixed = TRUE))
      result_list <- list()
      for (i in seq_along(transcripts)) {
        vals <- unlist(strsplit(transcripts[i], "|", fixed = TRUE))
        # Pad with empty strings if transcript has fewer fields than header
        if (length(vals) < length(fields)) {
          vals <- c(vals, rep("", length(fields) - length(vals)))
        }
        row <- setNames(as.list(vals), fields)
        row$transcript_index <- i
        result_list[[i]] <- row
      }
      return(result_list)
    }

    # Parsear CSQ para cada linha do VCF query
    cat("  Parseando campos CSQ de todos os transcripts...\n")
    parsed_all <- list()
    for (i in 1:nrow(vep_dt)) {
      tx_list <- parse_csq_all(vep_dt$csq[i], csq_fields)
      if (!is.null(tx_list)) {
        for (j in seq_along(tx_list)) {
          tx_list[[j]]$chr <- vep_dt$chr[i]
          tx_list[[j]]$pos <- vep_dt$pos[i]
          # NÃO incluir id/ref/alt do VCF para evitar colunas .x/.y duplicadas
          parsed_all[[length(parsed_all) + 1]] <- tx_list[[j]]
        }
      }
    }

    if (length(parsed_all) > 0) {
      vep_dt_full <- rbindlist(parsed_all, fill = TRUE)
      cat(sprintf("  Total de anotações VEP geradas (variantes x transcripts): %d\n", nrow(vep_dt_full)))
      
      # Criar coluna ID combinada para preservar unicidade por transcript
      vep_dt_full[, vep_id := sprintf("%s:%d_transcript_%d", chr, pos, transcript_index)]
       
      # Todos os CSQ fields + transcript_index (exclui chr, pos, vep_id)
      vep_cols <- setdiff(names(vep_dt_full), c("chr", "pos", "vep_id"))
      
      # MERGE COM TODOS OS TRANSCRIPTS (não apenas primeiro)
      dt <- merge(dt, vep_dt_full[, c("chr", "pos", vep_cols), with = FALSE], 
                  by = c("chr", "pos"), all.x = TRUE, suffixes = c("", ".vep"))
      sig_hits <- merge(sig_hits, vep_dt_full[, c("chr", "pos", vep_cols), with = FALSE], 
                        by = c("chr", "pos"), all.x = TRUE, suffixes = c("", ".vep"))
      sug_hits <- merge(sug_hits, vep_dt_full[, c("chr", "pos", vep_cols), with = FALSE], 
                        by = c("chr", "pos"), all.x = TRUE, suffixes = c("", ".vep"))
      
      n_anotados <- sum(!is.na(sig_hits$Consequence) & sig_hits$Consequence != "")
      cat(sprintf("  Variantes anotadas com sucesso: %d/%d\n", n_anotados, nrow(sig_hits)))
      cat(sprintf("  Transcriptos únicos processados: %d\n", length(unique(vep_dt_full$vep_id))))
    } else {
      cat("AVISO: Nenhum transcript VEP parseado com sucesso.\n")
    }

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
