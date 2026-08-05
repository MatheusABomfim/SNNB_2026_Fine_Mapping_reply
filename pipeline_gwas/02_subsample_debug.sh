#!/bin/bash
# Gera subsamples pequenos para debug rápido do pipeline GWAS
# Salva em $debug_dir (fora do git, ignorado pelo .gitignore)
set -euo pipefail

debug_dir="/storage4/matheusbomfim/SNNB_2026_Fine_mapping/gwas/debug"
cases_vcf="/storage4/matheusbomfim/programas/sarek/results_CA_Mama_DBGAP_New/variant_calling/haplotypecaller/joint_variant_calling/filtered/annotated/joint_germline.filtered.norm.AlphaMissense.vep.vcf.gz"
kg_bed="/storage4/matheusbomfim/quali/1kg/1kg_hg38_exome"
n_case_samps=50
n_kg_samps=50
tmpdir="${debug_dir}/_tmp"

mkdir -p "$debug_dir" "$tmpdir"
echo "=== Debug Subsample ==="
echo "  Case samples: $n_case_samps"
echo "  1KG samples: $n_kg_samps"
echo "  Output: $debug_dir"

# Detectar formato chr
chr_detect=$(bcftools query -f '%CHROM\n' "$cases_vcf" | head -1000 | grep -i "^chr" | head -1 || true)
if [ -n "$chr_detect" ]; then
  region_chr="chr22"
  region_num="22"
  echo "  Chromosome format: chr-prefixed"
else
  region_chr="22"
  region_num="22"
  echo "  Chromosome format: numeric"
fi
echo "  Region: $region_chr"

if [ ! -f "${debug_dir}/cases_debug.vcf.gz" ]; then
echo "--- Cases ---"
  n_total=$(bcftools query -l "$cases_vcf" | wc -l)
  bcftools query -l "$cases_vcf" | shuf -n "$n_case_samps" > "${debug_dir}/case_samps.txt"
  echo "  Total samples: $n_total, selected: $n_case_samps"
  bcftools view -r "$region_chr" -S "${debug_dir}/case_samps.txt" \
    "$cases_vcf" -Oz -o "${debug_dir}/cases_debug.vcf.gz"
  bcftools index --tbi "${debug_dir}/cases_debug.vcf.gz"
  n_vars=$(bcftools view -H "${debug_dir}/cases_debug.vcf.gz" | wc -l)
  n_samps=$(bcftools query -l "${debug_dir}/cases_debug.vcf.gz" | wc -l)
  echo "  Cases: $n_vars variants, $n_samps samples"
else
  echo "--- Cases (exists, skipping) ---"
fi

# Gerar lista de IDs dos cases (CHR:BP:REF:ALT) para extrair do 1KG
if [ ! -f "${debug_dir}/1kg_debug.bed" ]; then
echo "--- 1000 Genomes ---"
  # Converter cases VCF pra BED (mesmo processo do step 1) e extrair IDs
  if [ ! -f "${tmpdir}/cases_debug.bim" ]; then
    plink2 --vcf "${debug_dir}/cases_debug.vcf.gz" \
      --chr 1-22 --max-alleles 2 --rm-dup exclude-mismatch \
      --make-bed --out "${tmpdir}/cases_debug" --silent
    # Normalizar BIM (strip chr, ID = CHR:BP:REF:ALT, uppercase)
    awk -v OFS='\t' '{
      chr = $1; if (chr ~ /^chr/) chr = substr(chr, 4);
      $1 = chr; $2 = chr ":" $4 ":" $5 ":" $6; print
    }' "${tmpdir}/cases_debug.bim" | tr '[:lower:]' '[:upper:]' \
      > "${tmpdir}/cases_norm.bim"
    mv "${tmpdir}/cases_norm.bim" "${tmpdir}/cases_debug.bim"
  fi
  # Extrair IDs únicos (col2) do BIM
  awk '{print $2}' "${tmpdir}/cases_debug.bim" > "${debug_dir}/case_var_ids.txt"
  n_ids=$(wc -l < "${debug_dir}/case_var_ids.txt")
  echo "  Case variant IDs for extraction: $n_ids"

  # Extrair variantes correspondentes do 1KG
  echo "  Extracting matching variants from 1KG chr$region_num..."
  plink2 --bfile "$kg_bed" --chr "$region_num" \
    --extract "${debug_dir}/case_var_ids.txt" \
    --make-bed --out "${tmpdir}/kg_exact" --silent
  n_kg_vars=$(wc -l < "${tmpdir}/kg_exact.bim")
  echo "  1KG matching variants: $n_kg_vars"

  # Se ainda tiver muitas variantes, thinnar
  if [ "$n_kg_vars" -gt 5000 ]; then
    plink2 --bfile "${tmpdir}/kg_exact" \
      --thin-count 5000 --make-bed --out "${tmpdir}/kg_thin" --silent
    mv "${tmpdir}/kg_thin.bed" "${tmpdir}/kg_exact.bed"
    mv "${tmpdir}/kg_thin.bim" "${tmpdir}/kg_exact.bim"
    mv "${tmpdir}/kg_thin.fam" "${tmpdir}/kg_exact.fam"
    echo "  Thinned to 5000 variants"
  fi

  # Amostrar 50 samples
  awk '{print $1, $2}' "${tmpdir}/kg_exact.fam" | shuf -n "$n_kg_samps" \
    > "${debug_dir}/kg_samps.txt"
  plink2 --bfile "${tmpdir}/kg_exact" \
    --keep "${debug_dir}/kg_samps.txt" \
    --make-bed --out "${debug_dir}/1kg_debug" --silent
  n_vars=$(wc -l < "${debug_dir}/1kg_debug.bim")
  n_samps=$(wc -l < "${debug_dir}/1kg_debug.fam")
  echo "  1KG final: $n_vars variants, $n_samps samples"
else
  echo "--- 1000 Genomes (exists, skipping) ---"
fi

rm -rf "$tmpdir"
echo "=== Done ==="
echo "Run: qsub 04_gwas_assoc_debug.pbs"
