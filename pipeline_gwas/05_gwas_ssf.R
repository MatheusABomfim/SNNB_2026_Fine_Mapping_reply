library(data.table)

# GWAS-SSF (GWAS Summary Statistics Format) v1.1.0
# Spec: https://github.com/EBISPOT/gwas-summary-statistics-standard
# Converts PLINK2 --glm firth output to the standard TSV + YAML metadata format.
# EAF (effect allele frequency) is computed from 1000 Genomes Phase 3 via plink2 --freq.
#
# Mandatory SSF columns (in order):
#   chromosome, base_pair_location, effect_allele, other_allele,
#   odds_ratio (or beta or hazard_ratio), standard_error,
#   effect_allele_frequency, p_value
# Optional: variant_id, n

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# Production run: 459 cases (breast cancer), 459 1KG controls (1:1)
# Matching: PCAmatchR 1:1, 50 PCs, weighted Mahalanobis
# Final params: --firth --pca-match --match-k 1 --n-pcs 3
# Lambda GC: 1.2417, SNPs tested: 23458
# ── CONFIG — coorte phs000822 ──
BASE_DIR <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping"
KG_PREFIX <- "/storage4/matheusbomfim/quali/1kg/1kg_hg38_exome"

# ── CLI flags ──
#   --assoc <path>       PLINK2 assoc file (SRR_gwas.assoc)
#   --outdir <dir>       Output directory (default: dirname of assoc)
#   --kg-prefix <path>   1KG PLINK BED prefix
#   --amr-only           Use AMR super-population for EAF
#   --all (default)      Use all 1KG populations for EAF
args <- commandArgs(trailingOnly = TRUE)
assoc_file <- NULL
outdir <- NULL
kg_prefix <- KG_PREFIX
amr_only <- FALSE

i <- 1
while (i <= length(args)) {
  if (args[i] == "--assoc") {
    assoc_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--outdir") {
    outdir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--kg-prefix") {
    kg_prefix <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--amr-only") {
amr_only <- FALSE
    i <- i + 1
  } else if (args[i] == "--all") {
    amr_only <- FALSE
    i <- i + 1
  } else {
    i <- i + 1
  }
}

if (is.null(assoc_file) || !file.exists(assoc_file)) {
  stop("Usage: Rscript 05_gwas_ssf.R --assoc <assoc.tsv> [--outdir <dir>] [--kg-prefix <prefix>] [--amr-only|--all]")
}

if (is.null(outdir)) outdir <- dirname(assoc_file)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Case/control counts and sex from the association dataset fam ──
# target.{bed,bim,fam} is saved by 04_gwas_assoc.R in the same outdir.
# PLINK fam coding: pheno 2 = case, 1 = control; sex 1 = male, 2 = female, 0 = unknown.
fam_file <- file.path(outdir, "target.fam")
if (file.exists(fam_file)) {
  fam_meta <- fread(fam_file, col.names = c("FID", "IID", "PID", "MID", "Sex", "pheno"))
  n_cases <- sum(fam_meta$pheno == 2)
  n_controls <- sum(fam_meta$pheno == 1)
  if (all(fam_meta$Sex == 2)) {
    sex <- "female"
  } else if (all(fam_meta$Sex == 1)) {
    sex <- "male"
  } else if (any(fam_meta$Sex %in% c(1, 2))) {
    sex <- "combined"
  } else {
    # Sexo ausente no fam (PLINK2 nao carrega sexo de VCF): desenho do estudo e mulheres.
    sex <- "female"
  }
  cat(sprintf("  Samples from %s: %d cases, %d controls, sex = %s\n",
    fam_file, n_cases, n_controls, sex))
} else {
  # Fallback para o run de producao (target.fam ausente).
  cat("  WARNING: target.fam not found at", fam_file, "- using study defaults\n")
  n_cases <- 455L
  n_controls <- 910L
  sex <- "female"
}

cat("==========================================\n")
cat("GWAS-SSF v1.1.0 Summary Statistics Format\n")
cat("Start time:", format(Sys.time()), "\n")
cat("Assoc file:", assoc_file, "\n")
cat("Output dir:", outdir, "\n")
cat("1KG prefix:", kg_prefix, "\n")
cat("Population:", if (amr_only) "AMR only" else "ALL (default)", "\n")
cat("==========================================\n\n")

# ── 1. Read assoc file ──
# The cleaned assoc TSV contains: chr, pos, id, ref, alt, a1, test,
# obs_ct, beta, se, or, p, logp. Already filtered to autosomes 1-22,
# OR (0.01-100), valid P-values.
cat("[1/6] Reading assoc file...\n")
assoc <- fread(assoc_file, na.strings = c("NA", "nan", "NaN", ""))
cat("  Variants:", nrow(assoc), "\n")
cat("  Columns:", paste(names(assoc), collapse = ", "), "\n")

# ── 2. Write SNP IDs for --freq extraction ──
# The 1KG uses the same ID format (CHR:BP:REF:ALT) as the case VCF,
# normalized in the pipeline. plink2 --extract filters the 1KG PLINK
# BED to only the SNPs that passed QC, so --freq computes EAF on the
# exact same variant set used in the association.
cat("[2/6] Writing SNP IDs for --freq...\n")
snp_ids_file <- file.path(tempdir(), "gwas_ssf_snp_ids.txt")
fwrite(assoc[, .(id)], snp_ids_file, col.names = FALSE, quote = FALSE)
cat("  SNPs to extract:", nrow(assoc), "\n")

# ── 3. Run plink2 --freq on 1KG ──
# Effect allele frequency (EAF) is not output by PLINK2 --glm by default.
# We compute it post-hoc from the 1KG reference panel using --freq.
# Default (--all) uses all 2573 1KG samples; --amr-only restricts to AMR.
cat("[3/6] Running plink2 --freq on 1KG...\n")
freq_prefix <- file.path(tempdir(), "gwas_ssf_freq")
freq_args <- c("--bfile", kg_prefix, "--extract", snp_ids_file, "--freq", "--out", freq_prefix)
if (amr_only) {
  # Read PSAM to extract AMR sample IDs for --keep
  psam_file <- file.path(dirname(kg_prefix), "all_hg38.psam")
  if (file.exists(psam_file)) {
    psam <- fread(psam_file, col.names = c("fid", "pat", "mat", "sex", "superpop", "pop"))
    amr_ids <- psam[superpop == "AMR", .(fid)]
    amr_keep_file <- file.path(tempdir(), "gwas_ssf_amr_keep.txt")
    fwrite(amr_ids, amr_keep_file, sep = "\t", col.names = FALSE, quote = FALSE)
    freq_args <- c(freq_args, "--keep", amr_keep_file)
    cat("  Keeping AMR samples:", nrow(amr_ids), "\n")
  } else {
    cat("  Warning: PSAM not found at", psam_file, "- using ALL 1KG\n")
  }
}
system2("plink2", args = freq_args)

freq_file <- paste0(freq_prefix, ".afreq")
if (!file.exists(freq_file)) {
  stop("plink2 --freq failed: ", freq_file, " not found")
}
freq <- fread(freq_file)
cat("  Freq variants:", nrow(freq), "\n")

# ── 4. Merge freq with assoc ──
# plink2 --freq outputs ALT_FREQS (frequency of the ALT allele in the
# 1KG panel). We merge by variant ID and then map to the effect allele
# (a1, the tested allele in PLINK2 --glm).
cat("[4/6] Merging frequencies with assoc...\n")
setnames(freq, "ID", "id")
if (!"ALT_FREQS" %in% names(freq)) {
  stop("Expected ALT_FREQS column in --freq output. Got: ", paste(names(freq), collapse = ", "))
}
setnames(freq, "ALT_FREQS", "alt_freq")
freq <- freq[, .(id, ALT, alt_freq)]

ssf <- merge(assoc, freq, by = "id", all.x = TRUE, sort = FALSE)

# Map alt_freq to effect_allele_frequency:
#   If a1 == ALT (1KG alt allele): EAF = alt_freq
#   If a1 == ref:                  EAF = 1 - alt_freq
#   Otherwise (strand issue etc.): EAF = NA (will be written as #NA)
cat("[4/6] Mapping effect_allele_frequency...\n")
ssf[, effect_allele_frequency := fcase(
  a1 == ALT, alt_freq,
  a1 == ref, 1 - alt_freq,
  default = NA_real_
)]
cat("  Freq mapped:", sum(!is.na(ssf$effect_allele_frequency)), "/", nrow(ssf), "\n")

# ── 5. Map to GWAS-SSF v1.1 columns ──
# Required columns in order: chromosome, base_pair_location,
# effect_allele, other_allele, odds_ratio (or beta), standard_error,
# effect_allele_frequency, p_value. We add variant_id and n as extras.
cat("[5/6] Mapping to GWAS-SSF v1.1 columns...\n")

ssf[, chromosome := as.character(chr)]
ssf[, base_pair_location := pos]

ssf[, effect_allele := a1]
# other_allele is the non-tested allele (ref if a1==alt, alt if a1==ref)
ssf[, other_allele := fifelse(a1 == alt, ref, alt)]
ssf[, odds_ratio := or]
# Associação PLINK2 --glm firth reporta OR; beta = log(OR)
ssf[, "beta" := log(or)]
ssf[, standard_error := `log(or)_se`]
ssf[, p_value := p]
ssf[, variant_id := id]
# Ref allele (fixed reference allele from the genome assembly)
ssf[, ref_allele := fifelse(ref == a1, "EA", "OA")]
ssf[, n := obs_ct]

# Mandatory + optional SSF columns, in spec order
ssf_cols <- c(
  "chromosome", "base_pair_location", "effect_allele", "other_allele",
  "odds_ratio", "beta", "standard_error", "effect_allele_frequency",
  "p_value", "variant_id", "n", "ref_allele"
)
ssf_out <- ssf[, ..ssf_cols]

ssf_file <- file.path(outdir, "gwas_ssf.tsv")
# #NA is the GWAS-SSF standard placeholder for missing values
fwrite(ssf_out, ssf_file, sep = "\t", na = "#NA", quote = FALSE)
cat("  Saved:", ssf_file, "\n")

# ── 6. Generate YAML metadata ──
# The metadata file describes the study, trait, genotyping, samples,
# analysis parameters, and file integrity. It validates against the
# yamale schema at https://github.com/EBISPOT/gwas-summary-statistics-standard
cat("[6/6] Generating YAML metadata...\n")

trait_efo <- "EFO_0000305"
trait_desc <- "breast carcinoma"
samples_list <- list(
  list(
    # Case ancestry (nearest PC1-PC2 to 1KG): 18 AFR, 10 AMR, 427 EUR
    # Controls (910/2573 selected via PCAmatchR 1:2 from all 1KG)
    sample_ancestry_category = c("African", "Ad Mixed American", "European"),
    sample_size = n_cases + n_controls,
    case_control_study = TRUE,
    case_count = n_cases,
    control_count = n_controls,
    ancestry_method = c("genetically determined")
  )
)

file_md5 <- tools::md5sum(ssf_file)
meta <- list(
  gwas_id = NA_character_,
  author_notes = paste0(
    "GWAS of breast cancer in women (dbGaP phs000822.v1.p1) with 1000 Genomes Phase 3 external controls. ",
    "Cases are predominantly European by genetic ancestry (nearest-neighbor 1KG PC centroid); ",
    "small African and Ad Mixed American fractions. ",
    "Association: PLINK2 --glm firth --pca-match --match-k 1 --n-pcs 3. ",
    "QC: MAF 0.01, GENO 0.2, HWE 1e-6 (controls only), KING 0.088 (controls only)."
  ),
  trait_description = list(trait_desc),
  ontology_mapping = list(trait_efo),
  genome_assembly = "GRCh38",
  coordinate_system = "1-based",
  genotyping_technology = list(
    "Exome sequencing (cases)",
    "Whole-genome sequencing (1000 Genomes)"
  ),
  imputation_panel = "None (exome sequencing and WGS callset)",
  imputation_software = "None",
  samples = samples_list,
  sex = sex,
  data_file_name = basename(ssf_file),
  file_type = "GWAS-SSF v1.1",
  data_file_md5sum = file_md5,
  analysis_software = "PLINK2 --glm firth",
  adjusted_covariates = list("PC1", "PC2", "PC3"),
  # MAF 0.01 applied pre-association (Marees 2018; exome-adjusted)
  minor_allele_freq_lower_limit = 0.01,
  is_harmonised = FALSE,
  is_sorted = FALSE
)

# Write YAML metadata (manual serialization to avoid dependency on r-yaml)
meta_file <- paste0(ssf_file, "-meta.yaml")
yaml_lines <- c(
  "---",
  "# Study meta-data",
  paste0("gwas_id: ", meta$gwas_id %||% "null"),
  paste0("author_notes: \"", meta$author_notes, "\""),
  "",
  "# Trait Information",
  paste0("trait_description:"),
  paste0("  - \"", meta$trait_description[[1]], "\""),
  paste0("ontology_mapping:"),
  paste0("  - ", meta$ontology_mapping[[1]]),
  "",
  "# Genotyping Information",
  paste0("genome_assembly: ", meta$genome_assembly),
  paste0("coordinate_system: ", meta$coordinate_system),
  paste0("genotyping_technology:"),
  paste0("  - \"", meta$genotyping_technology[[1]], "\""),
  paste0("  - \"", meta$genotyping_technology[[2]], "\""),
  paste0("imputation_panel: \"", meta$imputation_panel, "\""),
  paste0("imputation_software: \"", meta$imputation_software, "\""),
  "",
  "# Sample Information",
  "samples:",
  "  - sample_ancestry_category:",
  paste0("      - \"", meta$samples[[1]]$sample_ancestry_category[1], "\""),
  paste0("      - \"", meta$samples[[1]]$sample_ancestry_category[2], "\""),
  paste0("      - \"", meta$samples[[1]]$sample_ancestry_category[3], "\""),
  paste0("    sample_size: ", meta$samples[[1]]$sample_size),
  paste0("    case_control_study: ", tolower(as.character(meta$samples[[1]]$case_control_study))),
  paste0("    case_count: ", meta$samples[[1]]$case_count),
  paste0("    control_count: ", meta$samples[[1]]$control_count),
  paste0("    ancestry_method:"),
  paste0("      - \"", meta$samples[[1]]$ancestry_method[[1]], "\""),
  paste0("sex: ", meta$sex),
  "",
  "# Summary Statistic information",
  paste0("data_file_name: ", meta$data_file_name),
  paste0("file_type: \"", meta$file_type, "\""),
  paste0("data_file_md5sum: ", meta$data_file_md5sum),
  paste0("analysis_software: \"", meta$analysis_software, "\""),
  paste0("adjusted_covariates:"),
  paste0("  - ", meta$adjusted_covariates[[1]]),
  paste0("  - ", meta$adjusted_covariates[[2]]),
  paste0("minor_allele_freq_lower_limit: ", meta$minor_allele_freq_lower_limit),
  "",
  "# Harmonization status",
  paste0("is_harmonised: ", tolower(as.character(meta$is_harmonised))),
  paste0("is_sorted: ", tolower(as.character(meta$is_sorted)))
)
writeLines(yaml_lines, meta_file)
cat("  Saved:", meta_file, "\n\n")

cat("=== GWAS-SSF Summary ===\n")
cat("Variants:", nrow(ssf_out), "\n")
cat("Missing EAF:", sum(is.na(ssf_out$effect_allele_frequency)), "/", nrow(ssf_out), "\n")
cat("MD5:", file_md5, "\n")
cat("==========================================\n")
