#!/usr/bin/env Rscript
# 04.3_associa_glimpse.R
# Associação dirigida dos 33 sinais do MVP nos casos do exoma phs000822 usando
# os genótipos imputados por GLIMPSE2 (04.2_imputa_glimpse.pbs), com controles
# 1KG matched (os mesmos do GWAS). Compara a direção do efeito com o MVP.
#
# Input :
#   - dados/mvp/bc_signal.csv                  (sinais + Beta/SE/PIP do MVP)
#   - resultados/glimpse/ligate/imp_*.bcf       (imputação GLIMPSE2 dos casos)
#   - resultados/glimpse/phased/kg_chr*.vcf.gz  (painel 1KG fasedo → controles)
#   - gwas/gwas_prod/tmp/pca.eigenvec           (PCs casos+controles)
#   - gwas/gwas_prod/tmp/keep_controls.txt      (IIDs dos controles matched)
# Output: resultados/tabelas/glimpse_associacao.csv
#
# Nota: os IDs dos casos são lidos do header do BCF ligado (samples do @RG dos
# BAMs, na ordem da imputação). A ordem das dosagens em FORMAT/DS é a mesma do
# header do BCF, então basta alinhar por índice.
#
# Nota de direção: a associação é feita para o alelo ALT do painel de referência
# (mesma orientação das dosagens imputadas). A comparação com o beta do MVP
# assume alelo de efeito do MVP = ALT (use a coluna MVP_A1 se disponível no CSV
# para alinhamento explícito de alelos).

# ─────────────────────────── CONFIG ───────────────────────────
MICROMAMBA <- "/storage2/matheusbomfim/projects/micromamba/bin/micromamba"
ENV_IMP    <- "imputacao"     # env com bcftools

MVP_SIGNAL    <- "dados/mvp/bc_signal.csv"
GLIMPSE_LIG   <- "resultados/glimpse/ligate"
GLIMPSE_PHASED<- "resultados/glimpse/phased"
BAMS_TXT      <- "resultados/glimpse/bams.txt"   # só fallback (path por linha)
PCA_EIGENVEC  <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/gwas_prod/tmp/pca.eigenvec"
KEEP_CONTROLS <- "/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/gwas_prod/tmp/keep_controls.txt"
OUT_TABELAS   <- "resultados/tabelas"
N_PCS         <- 3

# Colunas do bc_signal.csv
MVP_CHR  <- "CHR"
MVP_BP   <- "BP38"
MVP_RSID <- "RSID"
MVP_BETA <- "Beta Population"
MVP_SE   <- "SE Population"
MVP_PIP  <- "Overall PIP"
MVP_A1   <- ""            # coluna do alelo de efeito do MVP; "" = não disponível
# ───────────────────────────────────────────────────────────────

abs_path <- function(p) {
  if (grepl("^/", p)) p else file.path(REPO_ROOT, p)
}
args0 <- commandArgs(trailingOnly = FALSE)
fidx <- grep("^--file=", args0)
SCRIPT_DIR <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx])) else getwd()
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
MVP_SIGNAL <- abs_path(MVP_SIGNAL)
GLIMPSE_LIG <- abs_path(GLIMPSE_LIG)
GLIMPSE_PHASED <- abs_path(GLIMPSE_PHASED)
BAMS_TXT <- abs_path(BAMS_TXT)
OUT_TABELAS <- abs_path(OUT_TABELAS)

if (!dir.exists(OUT_TABELAS)) dir.create(OUT_TABELAS, recursive = TRUE)

bcftools <- function(...) {
  res <- system2(MICROMAMBA, args = c("run", "-n", ENV_IMP, "bcftools", ...),
                 stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0) {
    stop("bcftools falhou: ", paste(res, collapse = " "))
  }
  res
}

cat("== 04.3_associa_glimpse.R ==\n")
mvp <- read.csv(MVP_SIGNAL, stringsAsFactors = FALSE, check.names = FALSE)
mvp <- mvp[!is.na(mvp[[MVP_BP]]), ]
mvp$CHR_n <- gsub("^chr", "", as.character(mvp[[MVP_CHR]]), ignore.case = TRUE)
mvp$BP_n  <- as.numeric(mvp[[MVP_BP]])
cat(sprintf("Sinais MVP: %d\n", nrow(mvp)))

lig_files <- list.files(GLIMPSE_LIG, pattern = "\\.bcf$", full.names = TRUE)
if (length(lig_files) == 0) stop("Nenhum BCF ligado em ", GLIMPSE_LIG,
                                 " — rode 04.2_imputa_glimpse.pbs primeiro")

# IDs dos casos = samples do header do primeiro BCF ligado (ordem da imputação,
# mesma do FORMAT/DS). Fallback para bams.txt se o BCF não tiver samples.
case_ids <- tryCatch(bcftools("query", "-l", lig_files[1]),
                     error = function(e) character(0))
if (length(case_ids) == 0) {
  warning("BCF sem samples; usando bams.txt (basename) como IDs dos casos")
  bam_paths <- readLines(BAMS_TXT, warn = FALSE)
  case_ids <- sub("\\.[^.]*$", "", basename(bam_paths))
}
cat(sprintf("Casos: %d\n", length(case_ids)))

# ── 1. Dosagens imputadas nos casos (FORMAT/DS do GLIMPSE2) ──
dos_mat <- matrix(NA_real_, nrow = nrow(mvp), ncol = length(case_ids),
                  dimnames = list(mvp[[MVP_RSID]], case_ids))
for (i in seq_len(nrow(mvp))) {
  chr <- mvp$CHR_n[i]; pos <- mvp$BP_n[i]
  win_file <- lig_files[grepl(paste0("chr", chr, "_"), basename(lig_files))]
  for (f in win_file) {
    q <- bcftools("query", "-r", sprintf("chr%s:%d", chr, pos), "-f",
                  "%CHROM\\t%POS\\t[%DS\\t]\\n", f)
    if (length(q) >= 1 && !grepl("^$", q[1])) {
      cols <- strsplit(q, "\t", fixed = TRUE)[[1]]
      ds <- suppressWarnings(as.numeric(cols[-c(1, 2)]))
      if (length(ds) == length(case_ids)) { dos_mat[i, ] <- ds; break }
    }
  }
  if (all(is.na(dos_mat[i, ])))
    cat(sprintf("  aviso: %s (chr%s:%d) sem dosagem no GLIMPSE\n",
                mvp[[MVP_RSID]][i], chr, pos))
}
n_dos <- rowSums(!is.na(dos_mat))
cat(sprintf("Sítios com dosagem imputada: %d/%d\n", sum(n_dos > 0), nrow(mvp)))

# ── 2. Controles: genótipos do painel 1KG fasedo (GT → dosagem) ──
control_ids <- character(0)
if (file.exists(KEEP_CONTROLS)) {
  control_ids <- read.delim(KEEP_CONTROLS, header = FALSE, stringsAsFactors = FALSE)[[2]]
}
cat(sprintf("Controles (keep_controls.txt): %d\n", length(control_ids)))

ctrl_dos <- matrix(NA_real_, nrow = nrow(mvp), ncol = length(control_ids),
                   dimnames = list(mvp[[MVP_RSID]], control_ids))
if (length(control_ids) > 0) {
  ctrl_arg <- paste(control_ids, collapse = ",")
  for (i in seq_len(nrow(mvp))) {
    chr <- mvp$CHR_n[i]; pos <- mvp$BP_n[i]
    phased <- list.files(GLIMPSE_PHASED, pattern = paste0("kg_chr", chr, ".*\\.vcf\\.gz$"),
                         full.names = TRUE)[1]
    if (is.na(phased)) { cat("  aviso: painel fasedo chr", chr, "não encontrado\n"); next }
    q <- tryCatch(bcftools("query", "-s", ctrl_arg, "-r",
                           sprintf("chr%s:%d", chr, pos),
                           "-f", "%CHROM\\t%POS\\t[%GT\\t]\\n", phased),
                  error = function(e) character(0))
    if (length(q) >= 1 && !grepl("^$", q[1])) {
      cols <- strsplit(q, "\t", fixed = TRUE)[[1]]
      gt <- cols[-c(1, 2)]
      dos <- vapply(strsplit(gt, "[|/]", fixed = TRUE), function(a) {
        a <- suppressWarnings(as.numeric(a))
        if (any(is.na(a)) || length(a) < 2) NA_real_ else sum(a)
      }, numeric(1))
      if (length(dos) == length(control_ids)) ctrl_dos[i, ] <- dos
    }
  }
}
cat(sprintf("Sítios com genótipo dos controles: %d/%d\n",
            sum(rowSums(!is.na(ctrl_dos)) > 0), nrow(mvp)))

# ── 3. PCs (casos + controles) ─────────────────────────────────
pcs <- NULL
if (file.exists(PCA_EIGENVEC)) {
  eig <- read.table(PCA_EIGENVEC, header = TRUE, stringsAsFactors = FALSE,
                    comment.char = "#")
  pcs <- eig[, c("IID", paste0("PC", seq_len(N_PCS)))]
} else {
  cat("  aviso: pca.eigenvec não encontrado — sem covariáveis de ancestralidade\n")
}

# ── 4. Regressão logística caso×controle por sítio ────────────
results <- data.frame(
  RSID = mvp[[MVP_RSID]], CHR = mvp[[MVP_CHR]], BP38 = mvp[[MVP_BP]],
  Beta_MVP = mvp[[MVP_BETA]], SE_MVP = mvp[[MVP_SE]], PIP = mvp[[MVP_PIP]],
  n_casos = n_dos,
  stringsAsFactors = FALSE
)
results$n_controles <- rowSums(!is.na(ctrl_dos))
results$MAF_imputada_casos <- rowMeans(dos_mat, na.rm = TRUE) / 2
results$BETA_exoma <- NA_real_
results$SE_exoma   <- NA_real_
results$PVAL_exoma <- NA_real_
results$classificacao <- "Sem_imputacao"

if (length(control_ids) > 0 && !is.null(pcs)) {
  pheno <- data.frame(IID = c(case_ids, control_ids),
                      pheno = c(rep(1L, length(case_ids)), rep(0L, length(control_ids))))
  dat0 <- merge(pheno, pcs, by = "IID")
  for (i in seq_len(nrow(mvp))) {
    dos <- c(dos_mat[i, ], ctrl_dos[i, ])
    d <- dat0
    d$dos <- dos[match(d$IID, c(case_ids, control_ids))]
    d <- d[!is.na(d$dos), ]
    if (sum(d$pheno == 1) < 5 || sum(d$pheno == 0) < 5 || var(d$dos) == 0) {
      results$classificacao[i] <- "Sem_imputacao"
      next
    }
    fit <- tryCatch(glm(pheno ~ dos + PC1 + PC2 + PC3,
                        data = d[, c("pheno", "dos", paste0("PC", seq_len(N_PCS)))],
                        family = binomial), error = function(e) NULL)
    if (is.null(fit)) { results$classificacao[i] <- "Sem_imputacao"; next }
    co <- summary(fit)$coefficients["dos", ]
    results$BETA_exoma[i] <- co["Estimate"]
    results$SE_exoma[i]   <- co["Std. Error"]
    results$PVAL_exoma[i] <- co["Pr(>|z|)"]
    if (is.na(mvp[[MVP_BETA]][i])) {
      results$classificacao[i] <- "Indeterminada"
    } else if (sign(co["Estimate"]) == sign(mvp[[MVP_BETA]][i])) {
      results$classificacao[i] <- "Direcao_consistente"
    } else {
      results$classificacao[i] <- "Discordante"
    }
  }
}

# ── 5. Saída ──────────────────────────────────────────────────
cat("\n--- Classificação ---\n")
print(table(results$classificacao, useNA = "ifany"))
write.csv(results, file.path(OUT_TABELAS, "glimpse_associacao.csv"), row.names = FALSE)
cat("\nOK. Salvo em ", file.path(OUT_TABELAS, "glimpse_associacao.csv"), "\n", sep = "")
