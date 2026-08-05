library(data.table)
library(ggplot2)
set.seed(42)

project <- "SRR"
# ── CONFIG — coorte phs000822 ──
base_dir <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping"
vcf_file <- "/storage4/matheusbomfim/programas/sarek/results_CA_Mama_DBGAP_New/variant_calling/haplotypecaller/joint_variant_calling/filtered/annotated/joint_germline.filtered.norm.AlphaMissense.vep.vcf.gz"
kg_prefix <- "/storage4/matheusbomfim/quali/1kg/1kg_hg38_exome"
output_dir <- file.path(base_dir, "gwas", "pca_matching", format(Sys.time(), "%Y%m%d_%H%M%S"))

kg_dir <- dirname(kg_prefix)
kg_psam_dir <- kg_dir
tmp_dir <- file.path(output_dir, "tmp")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== PCA Matching - Knee Plot ===\n")
cat("Project:", project, "\n\n")

# ── 1. Convert VCF to PLINK BED ──
cat("[1/7] Converting VCF to PLINK BED...\n")
system2("plink2", args = c("--vcf", vcf_file,
  "--chr", "1-22",
  "--max-alleles", "2",
  "--rm-dup", "exclude-mismatch",
  "--make-bed", "--out", file.path(tmp_dir, "cases_raw"), "--silent"))

stopifnot(file.exists(paste0(file.path(tmp_dir, "cases_raw"), ".bed")))
stopifnot(file.exists(paste0(file.path(tmp_dir, "cases_raw"), ".bim")))
stopifnot(file.exists(paste0(file.path(tmp_dir, "cases_raw"), ".fam")))
n_init_vars <- nrow(fread(paste0(file.path(tmp_dir, "cases_raw"), ".bim")))
n_init_samps <- nrow(fread(paste0(file.path(tmp_dir, "cases_raw"), ".fam")))
cat("  Variants:", n_init_vars, "\n")
cat("  Samples :", n_init_samps, "\n")

bim <- fread(paste0(file.path(tmp_dir, "cases_raw"), ".bim"))
bim[, V1 := gsub("^chr", "", V1, ignore.case = TRUE)]
bim[, V2 := paste(V1, V4, toupper(V5), toupper(V6), sep = ":")]
fwrite(bim, paste0(file.path(tmp_dir, "cases_raw"), ".bim"),
  sep = "\t", col.names = FALSE, quote = FALSE)
cat("  IDs set to CHR:BP:REF:ALT\n\n")

# ── 2. Sample QC on cases ──
cat("[2/7] Sample QC on cases...\n")
cat("  Heterozygosity outliers (+/-3SD F_inbreeding)...\n")
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

system2("plink2", args = c("--bfile", file.path(tmp_dir, "cases_raw"),
  "--remove", file.path(tmp_dir, "het_outliers.txt"),
  "--make-bed", "--out", file.path(tmp_dir, "cases_qc"), "--silent"))

case_bed <- file.path(tmp_dir, "cases_qc")
n_vars <- nrow(fread(paste0(case_bed, ".bim")))
n_samps <- nrow(fread(paste0(case_bed, ".fam")))
cat("  Cases after sample QC:", n_samps, "samples x", n_vars, "variants\n\n")

# ── 3. Intersect SNPs with 1000G + merge ──
cat("[3/7] Intersecting SNPs with 1000 Genomes...\n")
case_bim <- fread(paste0(case_bed, ".bim"))
kg_bim <- fread(paste0(kg_prefix, ".bim"))
case_bim[, V2 := paste(V1, V4, toupper(V5), toupper(V6), sep = ":")]
kg_bim[, V2 := paste(V1, V4, toupper(V5), toupper(V6), sep = ":")]
common_snps <- intersect(case_bim$V2, kg_bim$V2)
n_common_snps <- length(common_snps)
cat("  Common SNPs:", n_common_snps, "/", nrow(case_bim),
  "(", round(100 * length(common_snps) / nrow(case_bim), 2), "% )\n")
if (n_common_snps == 0) {
  stop("No common SNPs between cases and 1KG.")
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

cat("  Filtering 1KG to common SNPs (all populations)...\n")
n_kg_vars_init <- nrow(kg_bim)
n_kg_samps_init <- nrow(fread(paste0(kg_prefix, ".fam")))
system2("plink2", args = c("--bfile", kg_prefix,
  "--extract", file.path(tmp_dir, "common_snps.txt"),
  "--make-bed", "--out", file.path(tmp_dir, "kg_common"), "--silent"))
kg_keep_bim <- fread(paste0(file.path(tmp_dir, "kg_common"), ".bim"))
kg_keep_fam <- fread(paste0(file.path(tmp_dir, "kg_common"), ".fam"))
n_kg_vars_common <- nrow(kg_keep_bim)
n_kg_samps_keep <- nrow(kg_keep_fam)
n_kg_removed <- n_kg_samps_init - n_kg_samps_keep
cat("  1KG variants: kept", n_kg_vars_common, "/", n_kg_vars_init,
  "(", round(100 * n_kg_vars_common / n_kg_vars_init, 2), "% )\n")
cat("  1KG samples: kept", n_kg_samps_keep, "/", n_kg_samps_init,
  "(removed", n_kg_removed, ",",
  round(100 * n_kg_removed / n_kg_samps_init, 2), "%)\n")

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
cat("[4/7] SNP QC on merged: MAF 0.01, GENO 0.2, SNPs-only...\n")
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
cat("[5/7] LD pruning (50 5 0.2)...\n")
system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_qc"),
  "--indep-pairwise", "50", "5", "0.2",
  "--out", file.path(tmp_dir, "merged_pruned"), "--silent"))
pruned <- fread(file.path(tmp_dir, "merged_pruned.prune.in"), header = FALSE)$V1
cat("  Pruned SNPs:", length(pruned), "\n")

# ── 6. PCA (10 components) + outlier removal ──
cat("[6/7] PCA (10 components)...\n")
system2("plink2", args = c("--bfile", file.path(tmp_dir, "merged_qc"),
  "--extract", file.path(tmp_dir, "merged_pruned.prune.in"),
  "--pca", "10", "--out", file.path(tmp_dir, "pca"), "--silent"))

pcs <- fread(file.path(tmp_dir, "pca.eigenvec"))
setnames(pcs, c("FID", "IID", paste0("PC", 1:10)))

pcs[, source := fifelse(IID %in% kg_pops$IID, "1KG", "case")]
n_case_pca <- sum(pcs$source == "case")
n_1kg_pca <- sum(pcs$source == "1KG")
cat("  Case samples:", n_case_pca, "| 1KG samples:", n_1kg_pca, "\n")

cat("  PCA outlier removal (+/-6SD PC1/PC2)...\n")
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

cat("  PCA files saved to", file.path(tmp_dir, "pca.eigenval"), "\n\n")

# ── 7. Knee plot (scree plot) ──
cat("[7/7] Generating knee plot...\n")
eigenval <- fread(file.path(tmp_dir, "pca.eigenval"), header = FALSE)$V1
var_exp <- eigenval / sum(eigenval) * 100
cum_var <- cumsum(var_exp)
n_pcs <- length(eigenval)

knee_data <- data.table(
  PC = factor(1:n_pcs, levels = 1:n_pcs),
  Variance = var_exp,
  Cumulative = cum_var
)

p <- ggplot(knee_data, aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = "#2171B5", alpha = 0.8) +
  geom_line(aes(y = Cumulative / max(Cumulative) * max(Variance), group = 1),
    color = "#D94801", linewidth = 1) +
  geom_point(aes(y = Cumulative / max(Cumulative) * max(Variance)),
    color = "#D94801", size = 2) +
  scale_y_continuous(
    name = "Variance Explained (%)",
    sec.axis = sec_axis(~ . / max(knee_data$Variance) * 100,
      name = "Cumulative Variance (%)")
  ) +
  labs(
    title = paste0("PCA Scree Plot - ", project),
    subtitle = sprintf("Total: %d samples (%d cases + %d controls)",
      nrow(pcs), n_case_pca, n_1kg_pca)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y.left = element_text(color = "#2171B5", size = 11),
    axis.title.y.right = element_text(color = "#D94801", size = 11),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

ggsave(file.path(output_dir, "pca_knee_plot.png"), p, width = 8, height = 5, dpi = 150)
cat("  Saved: pca_knee_plot.png\n\n")

# Print table
cat("  Variance explained by PC:\n")
for (i in 1:n_pcs) {
  cat(sprintf("    PC%d: %.2f%% (cumulative: %.2f%%)\n", i, var_exp[i], cum_var[i]))
}
cat("\n  Suggested K PCs (cumulative > 85%%):",
  min(which(cum_var > 85)), "\n")
cat("  Suggested K PCs (elbow, variance drop < 5%% of PC1):")
elbow <- which(var_exp[-1] < var_exp[1] * 0.05)[1]
if (is.na(elbow)) elbow <- n_pcs
cat(" ", elbow, "\n")

# Save eigenval/eigenvec copies in output_dir
file.copy(file.path(tmp_dir, "pca.eigenval"), file.path(output_dir, "pca.eigenval"), overwrite = TRUE)
file.copy(file.path(tmp_dir, "pca.eigenvec"), file.path(output_dir, "pca.eigenvec"), overwrite = TRUE)
cat("  PCA files copied to", output_dir, "\n")

unlink(tmp_dir, recursive = TRUE)
cat("=== PCA Matching Complete ===\n")
