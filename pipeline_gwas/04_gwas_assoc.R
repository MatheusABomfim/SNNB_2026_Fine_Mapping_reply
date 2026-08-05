library(data.table)
library(ggplot2)
set.seed(42)

# SNPRelate for KING kinship estimation
has_snprelate <- requireNamespace("SNPRelate", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
debug <- "--debug" %in% args
genesis <- "--genesis" %in% args || "--saige" %in% args
amr_only <- "--amr-only" %in% args
all_pops <- "--all" %in% args
if (amr_only && all_pops) stop("Cannot specify both --amr-only and --all")
if (!amr_only && !all_pops) amr_only <- TRUE
firth <- "--firth" %in% args
logistic <- "--logistic" %in% args
glm_hyb <- "--glm" %in% args
n_glm <- sum(firth, logistic, glm_hyb)
if (n_glm > 1) stop("Cannot specify more than one of --firth, --logistic, --glm")
if (n_glm == 0) glm_hyb <- TRUE
n_pcs <- 4
n_pcs_idx <- which(args == "--n-pcs")
if (length(n_pcs_idx) > 0) {
  if (n_pcs_idx == length(args)) stop("--n-pcs requires a numeric argument")
  n_pcs <- as.integer(args[n_pcs_idx + 1])
  if (is.na(n_pcs) || n_pcs < 1L || n_pcs > 20L)
    stop("--n-pcs must be between 1 and 20")
}
pca_match <- "--pca-match" %in% args
match_k <- 3
match_k_idx <- which(args == "--match-k")
if (length(match_k_idx) > 0) {
  if (match_k_idx == length(args)) stop("--match-k requires a numeric argument")
  match_k <- as.integer(args[match_k_idx + 1])
  if (is.na(match_k) || match_k < 1L || match_k > 10L)
    stop("--match-k must be between 1 and 10")
}
n_pca_pcs <- if (pca_match) 50L else max(10L, n_pcs)
pc_cols_assoc <- paste0("PC", 1:n_pcs)
covar_col_range <- if (n_pcs == 1L) "3" else paste0("3-", 2L + n_pcs)
covar_str <- paste0("PC1-PC", n_pcs)
has_genesis <- requireNamespace("GENESIS", quietly = TRUE) &&
               requireNamespace("GWASTools", quietly = TRUE)

# ── CONFIG — coorte phs000822 (adaptado do quali-workflow) ──
project <- "SRR"
# VCF joint-calling dos cases (Sarek + VEP/AlphaMissense)
vcf_input <- "/storage4/matheusbomfim/programas/sarek/results_CA_Mama_DBGAP_New/variant_calling/haplotypecaller/joint_variant_calling/filtered/annotated/joint_germline.filtered.norm.AlphaMissense.vep.vcf.gz"
# Base de outputs do projeto
base_dir <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping"
# Painel 1KG hg38 (já configurado pelo quali-workflow)
kg_panel <- "/storage4/matheusbomfim/quali/1kg/1kg_hg38_exome"

if (debug) {
  debug_dir <- file.path(base_dir, "gwas", "debug")
  vcf_file <- file.path(debug_dir, "cases_debug.vcf.gz")
  kg_prefix <- file.path(debug_dir, "1kg_debug")
  output_dir <- file.path(debug_dir, "out")
  cat("=== DEBUG MODE ===\n")
  cat("  vcf:", vcf_file, "\n")
  cat("  1kg:", paste0(kg_prefix, ".bed"), "\n")
  cat("  out:", output_dir, "\n\n")
} else {
  vcf_file <- vcf_input
  kg_prefix <- kg_panel
  output_dir <- file.path(base_dir, "gwas", format(Sys.time(), "%Y%m%d_%H%M%S"))
}

kg_dir <- dirname(kg_prefix)
kg_psam_dir <- if (debug) file.path(base_dir, "1kg") else kg_dir
tmp_dir <- file.path(output_dir, "tmp")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== GWAS Association - PLINK2 + 1000 Genomes ===\n")
cat("Project:", project, "\n\n")

# ── 1. Convert VCF to PLINK BED ──
cat("[1/13] Converting VCF to PLINK BED...\n")
system2("plink2", args = c("--vcf", vcf_file,
  "--chr", "1-22",
  "--max-alleles", "2",
  "--rm-dup", "exclude-mismatch",
  "--make-bed", "--out", file.path(tmp_dir, "cases_raw"), "--silent"))

# Debug: verify files and counts
stopifnot(file.exists(paste0(file.path(tmp_dir, "cases_raw"), ".bed")))
stopifnot(file.exists(paste0(file.path(tmp_dir, "cases_raw"), ".bim")))
stopifnot(file.exists(paste0(file.path(tmp_dir, "cases_raw"), ".fam")))
n_init_vars <- nrow(fread(paste0(file.path(tmp_dir, "cases_raw"), ".bim")))
n_init_samps <- nrow(fread(paste0(file.path(tmp_dir, "cases_raw"), ".fam")))
cat("  Variants:", n_init_vars, "\n")
cat("  Samples :", n_init_samps, "\n")

# Normalize BIM: strip chr prefix, set ID = CHR:BP:REF:ALT
bim <- fread(paste0(file.path(tmp_dir, "cases_raw"), ".bim"))
bim[, V1 := gsub("^chr", "", V1, ignore.case = TRUE)]
bim[, V2 := paste(V1, V4, toupper(V5), toupper(V6), sep = ":")]
fwrite(bim, paste0(file.path(tmp_dir, "cases_raw"), ".bim"),
  sep = "\t", col.names = FALSE, quote = FALSE)
cat("  IDs set to CHR:BP:REF:ALT\n\n")

# ── 2. Sample QC on cases ──
cat("[2/13] Sample QC on cases...\n")

cat("  Heterozygosity outliers (±3SD F_inbreeding)...\n")
system2("plink2", args = c("--bfile", file.path(tmp_dir, "cases_raw"),
  "--indep-pairwise", "50", "5", "0.2",
  "--out", file.path(tmp_dir, "cases_prune"), "--silent"))
system2("plink2", args = c("--bfile", file.path(tmp_dir, "cases_raw"),
  "--extract", file.path(tmp_dir, "cases_prune.prune.in"),
  "--het", "--out", file.path(tmp_dir, "cases_het"), "--silent"))
het <- fread(file.path(tmp_dir, "cases_het.het"))
het[, f_dev := abs(F - mean(F)) / sd(F)]
het_outliers <- het[f_dev > 3, IID]
n_het_removed <- length(het_outliers)
cat("    Het outliers removed:", n_het_removed, "(",
  round(100 * n_het_removed / n_init_samps, 2), "%)\n")
fwrite(data.table(FID = 0, IID = het_outliers),
  file.path(tmp_dir, "het_outliers.txt"),
  sep = "\t", quote = FALSE, col.names = FALSE)

# Create het-clean BED (cases_qc)
system2("plink2", args = c("--bfile", file.path(tmp_dir, "cases_raw"),
  "--remove", file.path(tmp_dir, "het_outliers.txt"),
  "--make-bed", "--out", file.path(tmp_dir, "cases_qc"), "--silent"))
rel_ok <- FALSE
rel_removed <- 0

case_bed <- file.path(tmp_dir, "cases_qc")
n_vars <- nrow(fread(paste0(case_bed, ".bim")))
n_samps <- nrow(fread(paste0(case_bed, ".fam")))
cat("  Cases after sample QC:", n_samps, "samples x", n_vars, "variants\n\n")

# ── 3. Intersect SNPs with 1000G + merge ──
cat("[3/13] Intersecting SNPs with 1000 Genomes...\n")
case_bim <- fread(paste0(case_bed, ".bim"))
kg_bim <- fread(paste0(kg_prefix, ".bim"))
# Reconstruct ID as CHR:BP:REF:ALT for consistent matching
case_bim[, V2 := paste(V1, V4, toupper(V5), toupper(V6), sep = ":")]
kg_bim[, V2 := paste(V1, V4, toupper(V5), toupper(V6), sep = ":")]
common_snps <- intersect(case_bim$V2, kg_bim$V2)
n_common_snps <- length(common_snps)
cat("  Common SNPs:", n_common_snps, "/", nrow(case_bim),
  "(", round(100 * length(common_snps) / nrow(case_bim), 2), "% )\n")
if (n_common_snps == 0) {
  stop("No common SNPs between cases and 1KG. Verify both are hg38 with matching chr/pos/alleles.")
}
writeLines(common_snps, file.path(tmp_dir, "common_snps.txt"))

cat("  Filtering to common SNPs...\n")
n_case_total <- nrow(case_bim)
system2("plink2", args = c("--bfile", case_bed,
  "--extract", file.path(tmp_dir, "common_snps.txt"),
  "--make-bed", "--out", file.path(tmp_dir, "cases_common"), "--silent"))
n_case_keep <- nrow(fread(paste0(file.path(tmp_dir, "cases_common"), ".bim")))
cat("  Cases: kept", n_case_keep, "/", n_case_total,
  "variants (", round(100 * n_case_keep / n_case_total, 2), "% )\n")

cat("  Loading 1000G population info...\n")
kg_psam <- fread(file.path(kg_psam_dir, "all_hg38.psam"),
  col.names = c("IID", "PAT", "MAT", "SEX", "SuperPop", "Population"))
kg_pops <- unique(kg_psam[, .(IID, SuperPop)])
amr_iids <- kg_psam[SuperPop == "AMR", IID]
cat("  Extracting common SNPs from 1KG (matching pool uses all populations)\n")
n_kg_vars_init <- nrow(kg_bim)
n_kg_samps_init <- nrow(fread(paste0(kg_prefix, ".fam")))
system2("plink2", args = c("--bfile", kg_prefix,
  "--extract", file.path(tmp_dir, "common_snps.txt"),
  "--make-bed", "--out", file.path(tmp_dir, "kg_common"), "--silent"))
kg_keep_bim <- fread(paste0(file.path(tmp_dir, "kg_common"), ".bim"))
kg_keep_fam <- fread(paste0(file.path(tmp_dir, "kg_common"), ".fam"))
n_kg_vars_common <- nrow(kg_keep_bim)
n_kg_samps_1kg <- nrow(kg_keep_fam)
cat("  1KG variants: kept", n_kg_vars_common, "/", n_kg_vars_init,
  "(", round(100 * n_kg_vars_common / n_kg_vars_init, 2), "% )\n")
cat("  1KG samples: kept", n_kg_samps_1kg, "/", n_kg_samps_init, "\n")

cat("  Applying HWE 1e-6 on 1KG controls...\n")
n_kg_vars_before <- nrow(kg_keep_bim)
system2("plink2", args = c("--bfile", file.path(tmp_dir, "kg_common"),
  "--hwe", "1e-6",
  "--make-bed", "--out", file.path(tmp_dir, "kg_hwe"), "--silent"))
unlink(paste0(file.path(tmp_dir, "kg_common"), c(".bed", ".bim", ".fam")))
file.rename(paste0(file.path(tmp_dir, "kg_hwe"), ".bed"),
            paste0(file.path(tmp_dir, "kg_common"), ".bed"))
file.rename(paste0(file.path(tmp_dir, "kg_hwe"), ".bim"),
            paste0(file.path(tmp_dir, "kg_common"), ".bim"))
file.rename(paste0(file.path(tmp_dir, "kg_hwe"), ".fam"),
            paste0(file.path(tmp_dir, "kg_common"), ".fam"))
kg_hwe_bim <- fread(paste0(file.path(tmp_dir, "kg_common"), ".bim"))
n_kg_hwe_removed <- n_kg_vars_before - nrow(kg_hwe_bim)
cat("  1KG after HWE:", nrow(kg_hwe_bim), "variants (removed",
  n_kg_hwe_removed, ")\n")

cat("  Merging datasets via bcftools...\n")
system2("plink2", args = c("--bfile", file.path(tmp_dir, "kg_common"),
  "--export", "vcf", "bgz", "--out", file.path(tmp_dir, "kg_vcf"), "--silent"))
system2("bcftools", args = c("index", "--tbi",
  file.path(tmp_dir, "kg_vcf.vcf.gz")))
system2("plink2", args = c("--bfile", file.path(tmp_dir, "cases_common"),
  "--export", "vcf", "bgz", "--out", file.path(tmp_dir, "cases_vcf"), "--silent"))
system2("bcftools", args = c("index", "--tbi",
  file.path(tmp_dir, "cases_vcf.vcf.gz")))
system2("bcftools", args = c("merge",
  file.path(tmp_dir, "kg_vcf.vcf.gz"),
  file.path(tmp_dir, "cases_vcf.vcf.gz"),
  "-Oz", "-o", file.path(tmp_dir, "merged.vcf.gz")))
system2("plink2", args = c("--vcf", file.path(tmp_dir, "merged.vcf.gz"),
  "--make-bed", "--out", file.path(tmp_dir, "merged"), "--silent"))

n_case_samps_common <- nrow(fread(paste0(file.path(tmp_dir, "cases_common"), ".fam")))
n_kg_samps_common   <- nrow(fread(paste0(file.path(tmp_dir, "kg_common"), ".fam")))
n_merged_vars <- nrow(fread(paste0(file.path(tmp_dir, "merged"), ".bim")))
n_merged_samps <- nrow(fread(paste0(file.path(tmp_dir, "merged"), ".fam")))
cat("  Merged: cases", n_case_samps_common, "+ 1KG", n_kg_samps_common,
  "=", n_merged_samps, "samples x", n_merged_vars, "variants\n")

# ── 4. SNP QC on merged (cases + 1000G) ──
cat("[4/13] SNP QC on merged: MAF 0.01, GENO 0.2, SNPs-only...\n")
system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged"),
  "--maf", "0.01", "--geno", "0.2", "--snps-only",
  "--make-bed", "--out", file.path(tmp_dir, "merged_qc"), "--silent"))
n_qc_vars <- nrow(fread(paste0(file.path(tmp_dir, "merged_qc"), ".bim")))
n_qc_samps <- nrow(fread(paste0(file.path(tmp_dir, "merged_qc"), ".fam")))
n_snpqc_vars <- n_qc_vars
n_snpqc_removed <- n_merged_vars - n_qc_vars
cat("  After SNP QC:", n_snpqc_vars, "variants (removed", n_snpqc_removed, ",",
  round(100 * n_snpqc_removed / n_merged_vars, 2), "%) |", n_qc_samps, "samples\n")

# ── 5. LD pruning for PCA ──
cat("[5/13] LD pruning (50 5 0.2)...\n")
system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_qc"),
  "--indep-pairwise", "50", "5", "0.2",
  if (debug) "--bad-ld",
  "--out", file.path(tmp_dir, "merged_pruned"), "--silent"))
pruned <- fread(file.path(tmp_dir, "merged_pruned.prune.in"), header = FALSE)$V1
cat("  Pruned SNPs:", length(pruned), "\n")

# ── 6. PCA (10 components) ──
cat("[6/13] PCA (", n_pca_pcs, " components)...\n", sep = "")
system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_qc"),
  "--extract", file.path(tmp_dir, "merged_pruned.prune.in"),
  "--pca", as.character(n_pca_pcs), "--out", file.path(tmp_dir, "pca"), "--silent"))

pcs <- fread(file.path(tmp_dir, "pca.eigenvec"))
setnames(pcs, c("FID", "IID", paste0("PC", 1:n_pca_pcs)))

pcs[, source := fifelse(IID %in% kg_pops$IID, "1KG", "case")]
n_case_pca <- sum(pcs$source == "case")
n_1kg_pca <- sum(pcs$source == "1KG")
cat("  Case samples:", n_case_pca, "| 1KG samples:", n_1kg_pca, "\n")

# ── 7. PCA outlier removal (±6SD PC1/PC2) ──
cat("[7/13] PCA outlier removal (±6SD PC1/PC2)...\n")
pc1_mean <- mean(pcs$PC1); pc1_sd <- sd(pcs$PC1)
pc2_mean <- mean(pcs$PC2); pc2_sd <- sd(pcs$PC2)
pcs[, outlier := abs(PC1 - pc1_mean) > 6 * pc1_sd |
                abs(PC2 - pc2_mean) > 6 * pc2_sd]
n_outliers <- sum(pcs$outlier)
cat("  Samples removed:", n_outliers, "| remaining:", nrow(pcs), "\n")
if (n_outliers > 0) {
  fwrite(pcs[outlier == TRUE, .(FID, IID)],
    file.path(tmp_dir, "pca_outliers.txt"),
    sep = "\t", quote = FALSE, col.names = FALSE)
}
pcs <- pcs[outlier == FALSE]
pcs[, outlier := NULL]

# ── 8. KING relatedness on merged (after PCA cleanup) ──
cat("[8/13] Cryptic relatedness (KING > 0.088, removing only controls)...\n")

# Remove PCA outliers from merged_qc
if (n_outliers > 0) {
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_qc"),
    "--remove", file.path(tmp_dir, "pca_outliers.txt"),
    "--make-bed", "--out", file.path(tmp_dir, "merged_king_in"), "--silent"))
} else {
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_qc"),
    "--make-bed", "--out", file.path(tmp_dir, "merged_king_in"), "--silent"))
}

# Sample missingness for tie-breaking
system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_king_in"),
  "--missing", "--out", file.path(tmp_dir, "merged_miss"), "--silent"))
miss <- fread(file.path(tmp_dir, "merged_miss.smiss"))
setnames(miss, c("FID", "IID", "nmiss_ctrl", "nmiss_case", "miss_rate"))

n_king_pairs <- 0
n_king_removed <- 0

if (has_snprelate) {
  library(SNPRelate)

  gds_file <- file.path(tmp_dir, "merged_king.gds")
  snpgdsBED2GDS(
    bed.fn = paste0(file.path(tmp_dir, "merged_king_in"), ".bed"),
    fam.fn = paste0(file.path(tmp_dir, "merged_king_in"), ".fam"),
    bim.fn = paste0(file.path(tmp_dir, "merged_king_in"), ".bim"),
    out.gdsfn = gds_file, snpfirstdim = FALSE, verbose = FALSE
  )
  gds <- snpgdsOpen(gds_file)
  on.exit({ snpgdsClose(gds); unlink(gds_file) })

  gds_snp_id <- read.gdsn(index.gdsn(gds, "snp.id"))
  pruned_match <- intersect(gds_snp_id, pruned)
  cat("    Pruned SNPs in GDS:", length(pruned_match), "\n")

  if (length(pruned_match) >= 50) {
    ibd <- snpgdsIBDKING(gds, snp.id = pruned_match,
      maf = 0.01, missing.rate = 0.1, num.thread = 1, verbose = FALSE)
    kin <- snpgdsIBDSelection(ibd)
    related_pairs <- kin[kin$kinship > 0.088, ]
    n_king_pairs <- nrow(related_pairs)
    cat("    Pairs with kin > 0.088:", n_king_pairs, "\n")

    if (n_king_pairs > 0) {
      source_map <- setNames(pcs$source, pcs$IID)
      miss_map <- setNames(miss$miss_rate, miss$IID)

      to_remove <- character()
      for (i in 1:n_king_pairs) {
        s1 <- as.character(related_pairs$ID1[i])
        s2 <- as.character(related_pairs$ID2[i])
        if (s1 %in% to_remove || s2 %in% to_remove) next

        src1 <- source_map[s1]
        src2 <- source_map[s2]

        if (src1 == "case" && src2 == "case") {
          next
        } else if (src1 == "case") {
          to_remove <- c(to_remove, s2)
        } else if (src2 == "case") {
          to_remove <- c(to_remove, s1)
        } else {
          if (miss_map[s1] >= miss_map[s2]) to_remove <- c(to_remove, s1)
          else to_remove <- c(to_remove, s2)
        }
      }

      n_king_removed <- length(unique(to_remove))
      cat("    KING: removing", n_king_removed, "samples\n")
      fwrite(data.table(FID = 0, IID = unique(to_remove)),
        file.path(tmp_dir, "rel_remove_merged.txt"),
        sep = "\t", quote = FALSE, col.names = FALSE)
    } else {
      cat("    KING: no related pairs found\n")
    }
  } else {
    cat("    Too few pruned SNPs, skipping KING\n")
  }
} else {
  if (!has_snprelate) cat("    SNPRelate not installed, skipping KING\n")
}

# Create merged_clean
n_before_king <- nrow(fread(paste0(file.path(tmp_dir, "merged_king_in"), ".fam")))
if (n_king_removed > 0) {
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_king_in"),
    "--remove", file.path(tmp_dir, "rel_remove_merged.txt"),
    "--make-bed", "--out", file.path(tmp_dir, "merged_clean"), "--silent"))
} else {
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_king_in"),
    "--make-bed", "--out", file.path(tmp_dir, "merged_clean"), "--silent"))
}
n_after_king <- nrow(fread(paste0(file.path(tmp_dir, "merged_clean"), ".fam")))
cat("    Samples before KING:", n_before_king, "| after KING:", n_after_king,
  "| total removed:", n_before_king - n_after_king, "\n")

# Filter PCs to samples that survived KING
merged_clean_ids <- fread(paste0(file.path(tmp_dir, "merged_clean"), ".fam"))$V2
pcs <- pcs[IID %in% merged_clean_ids]

unlink(paste0(file.path(tmp_dir, "merged_king_in"), c(".bed", ".bim", ".fam")))

# ── 9. Ancestry assignment ──
pcs_case <- pcs[source == "case"]
pcs_1kg <- pcs[source == "1KG"]
if (n_pcs > 0) {
n_ancestry_pcs <- min(4L, n_pcs)
cat("[9/13] Assigning ancestry to cases (nearest 1KG PC1-PC", n_ancestry_pcs, ")...\n", sep = "")
pc_cols <- paste0("PC", 1:n_ancestry_pcs)
kg_pcs_mat <- as.matrix(pcs_1kg[, ..pc_cols])
case_pcs_mat <- as.matrix(pcs_case[, ..pc_cols])
ancestry <- apply(case_pcs_mat, 1, function(case_pc) {
  d <- sqrt(rowSums(sweep(kg_pcs_mat, 2, case_pc, "-")^2))
  kg_pops$SuperPop[match(pcs_1kg$IID[which.min(d)], kg_pops$IID)]
})
pcs_case[, ancestry := ancestry]
pop_table <- table(ancestry)
pops_found <- names(pop_table)
cat("  Case ancestries:", paste(pop_table, names(pop_table), collapse = ", "), "\n")
ancestry_str <- paste(sprintf("%s: %d", names(pop_table), pop_table), collapse = ", ")
} else {
  ancestry_str <- "unavailable (n_pcs = 0)"
}

# Create list of QC-passing SNP IDs (used by both PLINK2 and GENESIS)
passing_snps <- fread(paste0(file.path(tmp_dir, "merged_qc"), ".bim"))$V2
writeLines(passing_snps, file.path(tmp_dir, "passing_snps.txt"))

# ── 10. Association ──
# Control selection: PCAmatchR (--pca-match) or population filter
if (pca_match) {
  if (!requireNamespace("PCAmatchR", quietly = TRUE))
    stop("PCAmatchR not installed. Run: install.packages('PCAmatchR')")
  if (!requireNamespace("optmatch", quietly = TRUE))
    stop("optmatch not installed. Run: install.packages('optmatch')")
  library(PCAmatchR)
  library(optmatch)
  setMaxProblemSize(size = Inf)

  cat("[10/13] PCAmatchR 1:", match_k, " (", n_pca_pcs, " PCs)...\n", sep = "")

  pc_cols_all <- paste0("PC", 1:n_pca_pcs)
  pcs_match <- as.data.frame(pcs[, c("IID", pc_cols_all), with = FALSE])
  cov_match <- as.data.frame(pcs[, .(IID, source)])
  cov_match$case <- ifelse(cov_match$source == "case", 1L, 0L)

  eigen_vals <- fread(file.path(tmp_dir, "pca.eigenval"), header = FALSE)$V1
  eigen_sum <- sum(eigen_vals)
  cat("    Eigenvalues sum:", round(eigen_sum, 1), "\n")

  n_cases_match <- nrow(pcs_case)
  n_pool <- nrow(pcs_1kg)
  max_k <- floor(n_pool / n_cases_match)
  effective_k <- if (match_k > max_k) {
    cat("    Warning: 1:", match_k, "K not possible (only", n_pool, "controls for", n_cases_match, "cases). Using 1:", max_k, "\n")
    max_k
  } else {
    match_k
  }

  matched <- match_maker(
    PC = pcs_match,
    eigen_value = eigen_vals,
    data = cov_match,
    ids = "IID",
    case_control = "case",
    num_controls = effective_k,
    eigen_sum = eigen_sum
  )

  matched_df <- matched$matches
  control_ids <- matched_df[matched_df$case == 0, "IID"]
  n_matched_cases <- matched_df[matched_df$case == 1, "IID"]
  cat("    Cases matched:", length(n_matched_cases),
      "| Controls matched:", length(control_ids), "\n")
  match_label <- paste0("PCAmatchR 1:", effective_k, " (", n_pca_pcs, " PCs)")
} else {
  control_ids <- if (amr_only) intersect(pcs_1kg$IID, amr_iids) else pcs_1kg$IID
  match_label <- if (amr_only) "none (1KG AMR only)" else "none (all 1KG samples)"
}

if (genesis) {
  cat("[10/13] Association (SAIGE SPA mixed model)...\n")
  n_cases <- nrow(pcs_case)
  n_controls <- length(control_ids)
  cat("    Cases:", n_cases, "| Controls:", n_controls, "\n")

  # a. Write pheno + covariates file for SAIGE (FID, IID, pheno, PC1..PCn)
  covar_ids <- c(pcs_case$IID, control_ids)
  saige_cols <- if (n_pcs > 0) c("FID", "IID", pc_cols_assoc) else c("FID", "IID")
  pheno_saige <- pcs[IID %in% covar_ids, ..saige_cols]
  pheno_saige[, pheno := fifelse(IID %in% pcs_case$IID, 1L, 0L)]
  fwrite(pheno_saige, file.path(tmp_dir, "saige_pheno.txt"), sep = "\t", quote = FALSE)

  # b. Create merged VCF for SAIGE
  fwrite(pcs_case[, .(FID, IID)], file.path(tmp_dir, "keep_cases.txt"),
    sep = "\t", quote = FALSE, col.names = FALSE)
  fwrite(data.table(FID = 0, IID = control_ids),
    file.path(tmp_dir, "keep_controls.txt"),
    sep = "\t", quote = FALSE, col.names = FALSE)

  system2("plink2", args = c("--bfile", file.path(tmp_dir, "cases_common"),
    "--extract", file.path(tmp_dir, "passing_snps.txt"),
    "--keep", file.path(tmp_dir, "keep_cases.txt"),
    "--make-bed", "--out", file.path(tmp_dir, "assoc_cases"), "--silent"))
  n_cases <- nrow(fread(paste0(file.path(tmp_dir, "assoc_cases"), ".fam")))

  system2("plink2", args = c("--bfile", file.path(tmp_dir, "kg_common"),
    "--extract", file.path(tmp_dir, "passing_snps.txt"),
    "--keep", file.path(tmp_dir, "keep_controls.txt"),
    "--make-bed", "--out", file.path(tmp_dir, "assoc_controls"), "--silent"))
  n_controls <- nrow(fread(paste0(file.path(tmp_dir, "assoc_controls"), ".fam")))

  cat("    Cases selected:", n_cases, "\n")
  cat("    Controls selected:", n_controls, "\n")

  fam_cases <- fread(paste0(file.path(tmp_dir, "assoc_cases"), ".fam"))
  fam_cases[, V6 := 2]
  fwrite(fam_cases, paste0(file.path(tmp_dir, "assoc_cases"), ".fam"),
    sep = " ", quote = FALSE, col.names = FALSE)

  system2("plink2", args = c("--bfile", file.path(tmp_dir, "assoc_controls"),
    "--export", "vcf", "bgz", "--out", file.path(tmp_dir, "assoc_ctrl_vcf"), "--silent"))
  system2("bcftools", args = c("index", "--csi",
    file.path(tmp_dir, "assoc_ctrl_vcf.vcf.gz")))
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "assoc_cases"),
    "--export", "vcf", "bgz", "--out", file.path(tmp_dir, "assoc_case_vcf"), "--silent"))
  system2("bcftools", args = c("index", "--csi",
    file.path(tmp_dir, "assoc_case_vcf.vcf.gz")))
  system2("bcftools", args = c("merge",
    file.path(tmp_dir, "assoc_ctrl_vcf.vcf.gz"),
    file.path(tmp_dir, "assoc_case_vcf.vcf.gz"),
    "-Oz", "-o", file.path(tmp_dir, "assoc.vcf.gz")))
  system2("bcftools", args = c("index", "--csi",
    file.path(tmp_dir, "assoc.vcf.gz")))

  n_assoc_vars <- nrow(fread(paste0(file.path(tmp_dir, "merged_qc"), ".bim")))
  n_assoc_samps <- n_cases + n_controls
  cat("    Association dataset:", n_assoc_samps, "samples x", n_assoc_vars, "variants\n")

  # c. Find SAIGE R scripts (bioconda instala no PATH)
  saige_scripts <- system("which step1_fitNULLGLMM.R", intern = TRUE, ignore.stderr = TRUE)
  if (length(saige_scripts) == 0) {
    # fallback: tenta via system.file
    saige_dir <- system.file("extdata", package = "SAIGE")
    step1_script <- file.path(saige_dir, "step1_fitNULLGLMM.R")
    step2_script <- file.path(saige_dir, "step2_SPAtests.R")
  } else {
    saige_dir <- dirname(saige_scripts[1])
    step1_script <- file.path(saige_dir, "step1_fitNULLGLMM.R")
    step2_script <- file.path(saige_dir, "step2_SPAtests.R")
  }
  if (!file.exists(step1_script)) {
    stop("SAIGE step1 script not found. Is the saige conda environment active?")
  }
  cat("    SAIGE scripts:", step1_script, "\n")

  # d. SAIGE step1: null GLMM
  cat("    SAIGE step1 (null GLMM)...\n")
  saige_step1_args <- c(
    step1_script,
    "--plinkFile", file.path(tmp_dir, "merged_clean"),
    "--phenoFile", file.path(tmp_dir, "saige_pheno.txt"),
    "--phenoCol", "pheno",
    "--sampleIDColinphenoFile", "IID",
    "--traitType", "binary",
    "--outputPrefix", file.path(tmp_dir, "saige_null"),
    "--nThreads", "4",
    "--IsOverwriteVarianceRatioFile", "TRUE"
  )
  if (n_pcs > 0) {
    saige_step1_args <- c(saige_step1_args,
      "--covarColList", paste0(pc_cols_assoc, collapse = ","))
  }
  system2("Rscript", args = saige_step1_args)

  # e. SAIGE step2: single-variant tests (per-chromosome, separate output files)
  if (!file.exists(file.path(tmp_dir, "saige_null.rda"))) {
    stop("SAIGE step1 model file not found. step1 may have failed.")
  }
  bim_assoc <- fread(paste0(file.path(tmp_dir, "assoc_cases"), ".bim"),
    select = 1, header = FALSE)
  vcf_chroms <- unique(as.character(bim_assoc$V1))
  cat("    Chromosomes in assoc VCF:", paste(vcf_chroms, collapse = ","), "\n")
  cat("    SAIGE step2 (single-variant tests per chromosome)...\n")
  step2_base <- c(
    step2_script,
    "--vcfFile", file.path(tmp_dir, "assoc.vcf.gz"),
    "--vcfFileIndex", file.path(tmp_dir, "assoc.vcf.gz.csi"),
    "--vcfField", "GT",
    "--GMMATmodelFile", file.path(tmp_dir, "saige_null.rda"),
    "--varianceRatioFile", file.path(tmp_dir, "saige_null.varianceRatio.txt"),
    "--LOCO", "FALSE",
    "--is_overwrite_output", "TRUE")
  saige_chr_files <- character(0)
  for (chr in vcf_chroms) {
    cat("      Chromosome", chr, "...\n")
    chr_file <- file.path(output_dir, paste0("saige_results_chr", chr, ".txt"))
    chr_args <- c(step2_base, "--chrom", chr, "--SAIGEOutputFile", chr_file)
    system2("Rscript", args = chr_args)
    if (file.exists(chr_file)) saige_chr_files <- c(saige_chr_files, chr_file)
  }

  # f. Merge per-chromosome results and map to common format
  if (length(saige_chr_files) == 0) {
    stop("No SAIGE step2 output files found (expected per-chromosome files)")
  }
  assoc <- rbindlist(lapply(saige_chr_files, fread,
    na.strings = c("NA", "nan", "NaN", "")))
  setnames(assoc, tolower(names(assoc)))
  if ("markerid" %in% names(assoc)) setnames(assoc, "markerid", "id")
  if ("snpid" %in% names(assoc)) setnames(assoc, "snpid", "id")
  if ("snp" %in% names(assoc) && !"id" %in% names(assoc)) setnames(assoc, "snp", "id")
  if ("p.value.spa" %in% names(assoc)) setnames(assoc, "p.value.spa", "p")
  if ("p.value" %in% names(assoc) && !"p" %in% names(assoc)) setnames(assoc, "p.value", "p")
  if (!"p" %in% names(assoc)) assoc[, p := NA_real_]
  if (!"or" %in% names(assoc) && "beta" %in% names(assoc)) {
    assoc[, or := exp(as.numeric(beta))]
  }
  if (!"or" %in% names(assoc)) assoc[, or := NA_real_]
  cols_num <- intersect(c("p", "or", "beta", "pos", "chr"), names(assoc))
  assoc[, (cols_num) := lapply(.SD, as.numeric), .SDcols = cols_num]

} else {
  # ── 10. Association via PLINK2 ──
  if (firth)       { glm_mod <- "firth";    glm_label <- "Firth" }
  else if (logistic) { glm_mod <- "no-firth"; glm_label <- "logistic" }
  else               { glm_mod <- "";         glm_label <- "hybrid (Firth fallback)" }
  cat("[10/13] Association (PLINK2 --glm", if (nchar(glm_mod) > 0) paste0(" ", glm_mod, ")...\n") else ")...\n")

  n_cases <- nrow(pcs_case)
  n_controls <- length(control_ids)
  cat("  Controls:", match_label, "|", n_cases, "cases,", n_controls, "controls\n")

  fwrite(pcs_case[, .(FID, IID)], file.path(tmp_dir, "keep_cases.txt"),
    sep = "\t", quote = FALSE, col.names = FALSE)
  fwrite(data.table(FID = 0, IID = control_ids),
    file.path(tmp_dir, "keep_controls.txt"),
    sep = "\t", quote = FALSE, col.names = FALSE)

  system2("plink2", args = c("--bfile", file.path(tmp_dir, "cases_common"),
    "--extract", file.path(tmp_dir, "passing_snps.txt"),
    "--keep", file.path(tmp_dir, "keep_cases.txt"),
    "--make-bed", "--out", file.path(tmp_dir, "assoc_cases"), "--silent"))
  n_cases <- nrow(fread(paste0(file.path(tmp_dir, "assoc_cases"), ".fam")))

  system2("plink2", args = c("--bfile", file.path(tmp_dir, "kg_common"),
    "--extract", file.path(tmp_dir, "passing_snps.txt"),
    "--keep", file.path(tmp_dir, "keep_controls.txt"),
    "--make-bed", "--out", file.path(tmp_dir, "assoc_controls"), "--silent"))
  n_controls <- nrow(fread(paste0(file.path(tmp_dir, "assoc_controls"), ".fam")))

  cat("    Cases selected:", n_cases, "\n")
  cat("    Controls selected:", n_controls, "\n")

  fam_cases <- fread(paste0(file.path(tmp_dir, "assoc_cases"), ".fam"))
  fam_cases[, V6 := 2]
  fwrite(fam_cases, paste0(file.path(tmp_dir, "assoc_cases"), ".fam"),
    sep = " ", quote = FALSE, col.names = FALSE)

  system2("plink2", args = c("--bfile", file.path(tmp_dir, "assoc_controls"),
    "--export", "vcf", "bgz", "--out", file.path(tmp_dir, "assoc_ctrl_vcf"), "--silent"))
  system2("bcftools", args = c("index", "--tbi",
    file.path(tmp_dir, "assoc_ctrl_vcf.vcf.gz")))
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "assoc_cases"),
    "--export", "vcf", "bgz", "--out", file.path(tmp_dir, "assoc_case_vcf"), "--silent"))
  system2("bcftools", args = c("index", "--tbi",
    file.path(tmp_dir, "assoc_case_vcf.vcf.gz")))
  system2("bcftools", args = c("merge",
    file.path(tmp_dir, "assoc_ctrl_vcf.vcf.gz"),
    file.path(tmp_dir, "assoc_case_vcf.vcf.gz"),
    "-Oz", "-o", file.path(tmp_dir, "assoc.vcf.gz")))
  system2("plink2", args = c("--vcf", file.path(tmp_dir, "assoc.vcf.gz"),
    "--make-bed", "--out", file.path(tmp_dir, "assoc"), "--silent"))
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "assoc"),
    "--extract", file.path(tmp_dir, "passing_snps.txt"),
    "--make-bed", "--out", file.path(tmp_dir, "assoc_snpqc"), "--silent"))

  keep_cases <- fread(file.path(tmp_dir, "keep_cases.txt"), header = FALSE)$V2
  fam_assoc <- fread(paste0(file.path(tmp_dir, "assoc_snpqc"), ".fam"))
  fam_assoc[, V6 := fifelse(V2 %in% keep_cases, 2, 1)]
  fwrite(fam_assoc, paste0(file.path(tmp_dir, "assoc_snpqc"), ".fam"),
    sep = " ", quote = FALSE, col.names = FALSE)

  n_assoc_vars <- nrow(fread(paste0(file.path(tmp_dir, "assoc_snpqc"), ".bim")))
  n_assoc_samps <- nrow(fam_assoc)
  cat("    Association dataset:", n_assoc_samps, "samples x", n_assoc_vars, "variants\n")

  glm_args <- if (nchar(glm_mod) > 0) c("--glm", glm_mod) else "--glm"
  covar_ids <- c(pcs_case$IID, control_ids)
  covar <- pcs[IID %in% covar_ids, c("FID", "IID", pc_cols_assoc), with = FALSE]
  fwrite(covar, file.path(tmp_dir, "covar.txt"), sep = " ", quote = FALSE)
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "assoc_snpqc"),
    glm_args,
    "--covar", file.path(tmp_dir, "covar.txt"),
    "--covar-col-nums", covar_col_range,
    "--out", file.path(output_dir, project),
    "--silent"))
}

# ── 11. Parse results + post-association filters + plots ──
cat("[11/13] Parsing results, filtering, and plotting...\n")

if (!genesis) {
  method_short <- "PLINK2 + 1000G"
  # PLINK2: detect and parse output file
  assoc_exts <- if (firth)
    c(".glm.firth", ".glm.logistic", ".glm.logistic.hybrid", ".glm")
  else if (logistic)
    c(".glm.logistic", ".glm.logistic.hybrid", ".glm.firth", ".glm")
  else
    c(".glm.logistic.hybrid", ".glm.logistic", ".glm.firth", ".glm")
  assoc_file <- ""
  for (ext in assoc_exts) {
    f <- file.path(output_dir, paste0(project, ".PHENO1", ext))
    if (file.exists(f)) {
      assoc_file <- f
      break
    }
  }

  if (!file.exists(assoc_file)) {
    cat("  No PLINK2 results found.\n")
    assoc <- data.table()
  } else {
    hdr <- readLines(assoc_file, n = 1)
    hdr <- sub("^#", "", hdr)
    cn <- strsplit(hdr, "\t")[[1]]
    assoc <- fread(assoc_file, skip = 1, header = FALSE)
    setnames(assoc, cn)

    if ("TEST" %in% names(assoc)) {
      assoc <- assoc[TEST == "ADD"]
    }

    old <- names(assoc)
    new <- tolower(gsub("[#\\?]", "", old))
    setnames(assoc, old, new)

    if ("chrom" %in% names(assoc) && !"chr" %in% names(assoc)) {
      setnames(assoc, "chrom", "chr")
    }

    if (!"or" %in% names(assoc) && "beta" %in% names(assoc)) {
      assoc[, or := exp(beta)]
    }

    if ("or" %in% names(assoc)) {
      assoc[, or := as.numeric(or)]
    }
  }
} else {
  cat("  Using SAIGE results from step 10...\n")
  method_short <- "SAIGE + 1000G"
  model_str <- "SAIGE (SPA mixed model)"
}

# Post-association filters (common to both engines)
if (nrow(assoc) > 0) {
  n_before <- nrow(assoc)

  assoc[, chr := as.character(chr)]
  assoc[grepl("^chr", chr, ignore.case = TRUE), chr := sub("^chr", "", chr, ignore.case = TRUE)]
  assoc <- assoc[!grepl("[^0-9]", chr)]
  assoc[, chr := as.numeric(chr)]
  assoc <- assoc[!is.na(chr) & chr >= 1 & chr <= 22]
  cat("  After CHR filter (autosomes 1-22):", nrow(assoc), "(removed", n_before - nrow(assoc), ")\n")

  if (nrow(assoc) > 0) {
    n_before <- nrow(assoc)
    assoc <- assoc[(or > 0.01 & or < 100) | is.na(or)]
    cat("  After OR filter (0.01-100):", nrow(assoc), "(removed", n_before - nrow(assoc), ")\n")
  }

  if (nrow(assoc) > 0) {
    assoc <- assoc[!is.na(p) & is.finite(p) & p > 0 & p <= 1]
    assoc[, logp := -log10(p)]
  }

  n_results <- nrow(assoc)
  n_sig <- if (n_results > 0) sum(assoc$p < 5e-8, na.rm = TRUE) else 0
  n_sugg <- if (n_results > 0) sum(assoc$p < 1e-5, na.rm = TRUE) else 0
  min_p <- if (n_results > 0) min(assoc$p, na.rm = TRUE) else NA_real_
  top <- if (n_results > 0) assoc[which.min(p)] else data.table()
  cat("\n  Results:", n_results, "variants\n")
  cat("  Significant (p < 5e-8):", n_sig, "\n")
  cat("  Suggestive (p < 1e-5):", n_sugg, "\n\n")

  # Manhattan plot
  if (n_results > 0) {
    if (!requireNamespace("qqman", quietly = TRUE)) {
      stop("qqman not installed")
    }
    library(qqman)
    pdf(file.path(output_dir, "manhattan_plot.pdf"), width = 12, height = 6)
    manhattan(assoc, chr = "chr", bp = "pos", p = "p", snp = "id",
      suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8),
      col = c("#08306B", "#2171B5", "#6BAED6", "#9ECAE1"),
      main = paste0("GWAS - ", project, " (", method_short, ")"))
    dev.off()
    cat("  Saved: manhattan_plot.pdf\n")

    # QQ-plot with lambda
    lambda <- round(median(qchisq(1 - assoc$p, 1), na.rm = TRUE) /
                    median(qchisq(1 - ppoints(length(assoc$p)), 1), na.rm = TRUE), 4)
    pdf(file.path(output_dir, "qq_plot.pdf"), width = 6, height = 6)
    qq(assoc$p, main = paste0("QQ-Plot - ", project, " (", method_short, ")\nlambda = ", lambda))
    dev.off()
    cat("  Saved: qq_plot.pdf\n")
    cat("  Lambda GC:", lambda, "\n")
  }

  if (!exists("lambda")) lambda <- NA_real_

  # Save results
  cat("\nSaving results...\n")
  fwrite(assoc, file.path(output_dir, paste0(project, "_gwas.assoc")),
    sep = "\t", quote = FALSE)
  saveRDS(assoc, file.path(output_dir, "gwas_results.rds"))
  cat("  GWAS .assoc:", file.path(output_dir, paste0(project, "_gwas.assoc")), "\n")
  cat("  GWAS .rds:  ", file.path(output_dir, "gwas_results.rds"), "\n\n")

  # ── 12. LD Score regression ──
  ldsc_intercept <- NA_real_
  ldsc_ratio <- NA_real_
  # LDSC (opcional): aponta para instalação/LD scores do quali-workflow.
  # Se indisponível, o step 12 é pulado graciosamente.
  ldsc_py <- "/storage4/matheusbomfim/tools/ldsc/ldsc.py"
  ldscores_prefix <- "/storage4/matheusbomfim/quali/1kg_ld_scores/exome_chr"
  ldsc_ready <- file.exists(ldsc_py) && file.exists(paste0(ldscores_prefix, "1.l2.ldscore.gz"))

  if (ldsc_ready && !is.null(assoc) && nrow(assoc) > 0) {
    cat("[12/13] LD Score regression...\n")
    cat("  [WARN] LDSC was designed for genome-wide data (>200k SNPs).\n")
    cat("  Exome panels (~", nrow(assoc), " SNPs) violate LDSC assumptions —\n", sep = "")
    cat("  intercept and attenuation ratio should be interpreted with\n")
    cat("  caution or ignored entirely.\n")
    ldsc_file <- file.path(output_dir, "ldsc_sumstats.txt")

    if ("z_stat" %in% names(assoc)) {
      ldsc_df <- assoc[, .(
        SNP = id, CHR = chr, BP = pos,
        A1 = a1, A2 = ifelse(a1 == ref, alt, ref),
        Z = z_stat, N = obs_ct, P = p
      )]
      fwrite(ldsc_df, ldsc_file, sep = "\t", quote = FALSE)

      system2("python", args = c(ldsc_py,
        "--h2", ldsc_file,
        "--ref-ld-chr", ldscores_prefix,
        "--w-ld-chr", ldscores_prefix,
        "--out", file.path(output_dir, "ldsc")
      ))

      ldsc_log_file <- file.path(output_dir, "ldsc.log")
      if (file.exists(ldsc_log_file)) {
        log_lines <- readLines(ldsc_log_file)
        iline <- grep("Intercept", log_lines, value = TRUE)
        rline <- grep("Ratio", log_lines, value = TRUE)
        if (length(iline)) ldsc_intercept <- as.numeric(gsub(".*Intercept[^0-9.]*([0-9.]+).*", "\\1", iline[1]))
        if (length(rline)) ldsc_ratio <- as.numeric(gsub(".*Ratio[^0-9.]*([0-9.]+).*", "\\1", rline[1]))
      }
      cat("  LDSC intercept:", ldsc_intercept, "\n")
      cat("  LDSC atten.ratio:", ldsc_ratio, "\n")
    } else {
      cat("  LDSC skipped: no Z_STAT column\n")
    }
  } else {
    if (!ldsc_ready) cat("  LDSC skipped: ldsc.py or LD scores not found\n")
  }

  # ── Generate report ──
  top_id <- if (n_results > 0) top$id else NA
  top_chr <- if (n_results > 0) top$chr else NA
  top_pos <- if (n_results > 0) top$pos else NA
  top_or <- if (n_results > 0) round(top$or, 4) else NA
  top_p <- if (n_results > 0) format(top$p, scientific = TRUE, digits = 3) else NA

  report_file <- file.path(output_dir, "gwas_report.txt")
  report <- c(
    paste0("Project: ", project),
    "",
    "=== Pipeline Summary ===",
    "QC filters:             het +/-3SD, KING (kin > 0.088), hwe 1e-6, maf 0.01, geno 0.2",
    paste0("Matching:               ", match_label),
    paste0("Association:            ", if (genesis) model_str else paste0("PLINK2 --glm", if (nchar(glm_mod) > 0) paste0(" ", glm_mod) else "", " (", glm_label, ")")),
    paste0("Covariates:             ", covar_str),
    "Post-assoc filters:     OR 0.01-100, CHR 1-22, P valid",
    "",
    "=== Sample QC (Cases) ===",
    paste0("Initial samples:                        ", n_init_samps),
    paste0("After het +/-3SD:                         ", n_init_samps - n_het_removed,
      "  (removed ", n_het_removed, ", ", round(100 * n_het_removed / n_init_samps, 2), "%)"),
    paste0("After KING (kin > 0.088):               Skipped (KING on merged)"),
    paste0("Cases after sample QC:                   ", n_samps),
    "",
    "=== Sample QC (1KG Controls) ===",
    paste0("Initial samples:                         ", n_kg_samps_init),
    paste0("Samples retained:                         ", n_kg_samps_1kg),
    paste0("Variants removed by HWE 1e-6:             ", n_kg_hwe_removed),
    "",
    "=== Merge + Variant QC ===",
    paste0("Initial case variants:                    ", n_init_vars),
    paste0("Common SNPs (case intersect 1KG):         ", n_common_snps,
      "  (", round(100 * n_common_snps / n_init_vars, 2), "%)"),
    paste0("After SNP QC (MAF 0.01, GENO 0.2):        ", n_snpqc_vars,
      "  (removed ", n_snpqc_removed, ", ", round(100 * n_snpqc_removed / n_merged_vars, 2), "%)"),
    paste0("Pruned for PCA:                           ", length(pruned)),
    paste0("KING related (kin > 0.088, ctrl only):    ", n_king_removed, " removed (", n_king_pairs, " pairs)"),
    "",
    "=== PCA ===",
    paste0("PCA outliers removed:                     ", n_outliers),
    "",
    "=== Association ===",
    paste0("Cases selected:                           ", n_cases),
    paste0("Controls (1KG", if (amr_only) " AMR" else "", " samples):                   ", n_controls),
    paste0("Total:                                    ", n_cases + n_controls),
    paste0("Ancestry distribution:                    ", ancestry_str),
    "",
    "=== Results ===",
    paste0("SNPs tested:                              ", n_results),
    paste0("Significant (p < 5e-8):                   ", n_sig),
    paste0("Suggestive (p < 1e-5):                    ", n_sugg),
    paste0("Minimum p-value:                          ", min_p),
    paste0("Lambda GC:                                ", lambda),
    paste0("LDSC intercept:                           ", ldsc_intercept),
    paste0("LDSC attenuation ratio:                   ", ldsc_ratio),
    "",
    paste0("Top SNP:                                  ", top_id),
    if (n_results > 0) paste0("  CHR:                                    ", top_chr) else character(),
    if (n_results > 0) paste0("  POS:                                    ", top_pos) else character(),
    if (n_results > 0) paste0("  OR:                                     ", top_or) else character(),
    if (n_results > 0) paste0("  P-value:                                ", top_p) else character()
  )
  writeLines(report, report_file)
  cat("  Report:", report_file, "\n")
} else {
  cat("  No results found.\n")
}

# ── 13. Save target BED for PRS transferability ──
target_bed <- file.path(output_dir, "target")
if (!genesis && file.exists(paste0(file.path(tmp_dir, "assoc_snpqc"), ".bed"))) {
  cat("[13/13] Saving target BED:", target_bed, "\n")
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "assoc_snpqc"),
    "--make-bed", "--out", target_bed, "--silent"))
} else if (genesis && file.exists(paste0(file.path(tmp_dir, "merged_clean"), ".bed"))) {
  cat("[13/13] Saving target BED:", target_bed, "\n")
  system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_clean"),
    "--make-bed", "--out", target_bed, "--silent"))
}

unlink(tmp_dir, recursive = TRUE)
cat("=== GWAS Association Complete ===\n")
