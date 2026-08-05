#!/usr/bin/env Rscript
# 05.5_prepara_exoma.R
# Converte o output GWAS-SSF (05_gwas_ssf.R) para o formato de summary
# statistics esperado pelo workflow MVP x exoma (scripts/02_lookup_exoma.R).
#
# Input :
#   - gwas_ssf.tsv (output do 05_gwas_ssf.R, colunas GWAS-SSF v1.1)
#   - VCF joint-calling dos cases (Sarek) — usado só para anotar rsID
# Output: dados/exoma/exoma_sumstats.txt
#         colunas: rsid chr bp a1 a2 maf beta se pval n

# ─────────────────────────── CONFIG ───────────────────────────
SSF_TSV <- "gwas_ssf.tsv"                       # relativo ao OUTDIR
OUTDIR  <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/gwas_prod"
OUT_FILE <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/gwas_prod/exoma_sumstats.txt"
# Copiar também para dados/exoma/ do repo (descomente se quiser)
# OUT_FILE_REPO <- "dados/exoma/exoma_sumstats.txt"

# VCF para anotar rsID (ID column). "" = pular anotação rsID.
CASE_VCF <- "/storage4/matheusbomfim/programas/sarek/results_CA_Mama_DBGAP_New/variant_calling/haplotypecaller/joint_variant_calling/filtered/annotated/joint_germline.filtered.norm.AlphaMissense.vep.vcf.gz"
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))

cat("== 05.5_prepara_exoma.R ==\n")
ssf_file <- if (dirname(SSF_TSV) == ".") file.path(OUTDIR, SSF_TSV) else SSF_TSV
if (!file.exists(ssf_file)) {
  # Fallback: procurar em OUTDIR
  alt <- file.path(OUTDIR, "gwas_ssf.tsv")
  if (file.exists(alt)) ssf_file <- alt else stop("gwas_ssf.tsv não encontrado em ", OUTDIR)
}
cat("Lendo:", ssf_file, "\n")
ssf <- fread(ssf_file, na.strings = c("#NA", "NA", ""))

required <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
              "beta", "standard_error", "effect_allele_frequency", "p_value")
missing <- setdiff(required, names(ssf))
if (length(missing) > 0) stop("Colunas SSF ausentes: ", paste(missing, collapse = ", "))

# ── Mapear para o formato do 02_lookup_exoma.R ──
exoma <- ssf[, .(
  CHR  = as.character(chromosome),
  BP   = as.integer(base_pair_location),
  A1   = effect_allele,
  A2   = other_allele,
  MAF  = fifelse(effect_allele_frequency <= 0.5,
                 effect_allele_frequency, 1 - effect_allele_frequency),
  BETA = as.numeric(beta),
  SE   = as.numeric(standard_error),
  PVAL = as.numeric(p_value),
  N    = if ("n" %in% names(ssf)) n else NA_integer_
)]

# ── Anotar rsID a partir do ID column do VCF ──
exoma[, RSID := ""]
if (CASE_VCF != "" && file.exists(CASE_VCF)) {
  cat("Anotando rsID a partir do VCF...\n")
  # bcftools query extrai CHROM, POS, REF, ALT, ID em uma linha por variante.
  # Para variantes multi-alélicas a ID se aplica à variante; usamos CHROM:POS:REF:ALT
  # para casar com variant_id do SSF (mesma normalização do quali-workflow).
  vcf_tmp <- tempfile(fileext = ".tsv")
  system2("bcftools", args = c("query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\n",
                               CASE_VCF), stdout = vcf_tmp)
  vcf_ann <- fread(vcf_tmp, header = FALSE,
                   col.names = c("chrom", "pos", "ref", "alt", "id"))
  # normalizar alelos (toupper, como no quali-workflow) e remover "chr" do cromossomo
  vcf_ann[, chrom := gsub("^chr", "", chrom, ignore.case = TRUE)]
  vcf_ann[, ref := toupper(ref)]
  vcf_ann[, alt := toupper(alt)]
  vcf_ann[, varid := paste(chrom, pos, ref, alt, sep = ":")]
  # pegar a primeira ID (rs) por variante
  vcf_ann <- vcf_ann[id != ".", .(varid, id = unlist(tstrsplit(id, "[,;]", keep = 1))[1]), by = varid]
  vcf_ann <- unique(vcf_ann, by = "varid")

  exoma[, varid := paste(CHR, BP, A1, A2, sep = ":")]
  exoma[vcf_ann, RSID := i.id, on = "varid"]
  exoma[, varid := NULL]
  unlink(vcf_tmp)
  cat("  rsID anotados:", sum(exoma$RSID != ""), "/", nrow(exoma), "\n")
} else {
  warning("CASE_VCF não encontrado — coluna RSID vazia. Use match por posição no 02_lookup_exoma.R ou anote depois.")
}

# ── Ordenar colunas e salvar ──
setcolorder(exoma, c("RSID", "CHR", "BP", "A1", "A2", "MAF", "BETA", "SE", "PVAL", "N"))
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
fwrite(exoma, OUT_FILE, sep = "\t", quote = FALSE)
cat("Salvo:", OUT_FILE, "(", nrow(exoma), "linhas)\n")

if (exists("OUT_FILE_REPO") && OUT_FILE_REPO != "") {
  dir.create(dirname(OUT_FILE_REPO), showWarnings = FALSE, recursive = TRUE)
  file.copy(OUT_FILE, OUT_FILE_REPO, overwrite = TRUE)
  cat("Copiado para:", OUT_FILE_REPO, "\n")
}

cat("\nColunas (formato do 02_lookup_exoma.R):\n")
print(names(exoma))
cat("\nOK. Agora rode: Rscript scripts/02_lookup_exoma.R\n")
