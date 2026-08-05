#!/bin/bash
# setup_1kg.sh — Baixa 1000 Genomes Phase 3 (hg38) e processa autossomos
# Uso: bash scripts/setup_1kg.sh
# Pré-requisito: micromamba run -n quali (plink2, zstd, wget no PATH)
# NOTA (SNNB_2026): o painel 1KG já existe em /storage4/matheusbomfim/quali/1kg
# e é reutilizado (01_setup_1kg e 01.1_validate_setup podem ser SKIP).
set -euo pipefail

BASE="/storage4/matheusbomfim/quali"
mkdir -p "$BASE/1kg"
cd "$BASE/1kg"

# 1. 1000 Genomes Phase 3 (hg38) em PGEN
if [ ! -f all_hg38.pgen ]; then
  echo "[1/5] Downloading 1000 Genomes Phase 3 (hg38)..."
  wget -q -O all_hg38.pgen.zst "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1"
  wget -q -O all_hg38_rs.pvar.zst "https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1"
  wget -q -O all_hg38.psam "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1"
  wget -q -O deg2_hg38.king.cutoff.out.id "https://www.dropbox.com/s/4zhmxpk5oclfplp/deg2_hg38.king.cutoff.out.id?dl=1"
  echo "  Decompressing..."
  plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen
  zstd -d -f all_hg38_rs.pvar.zst -o all_hg38.pvar
  rm -f all_hg38.pgen.zst all_hg38_rs.pvar.zst
  echo "  OK"
else
  echo "[1/5] 1000 Genomes already exists"
fi

# 2. Processar autossomos, bialélicos, sem parentais
echo "[2/5] Processing autosomes (biallelic, no relatives)..."
plink2 \
  --pfile all_hg38 \
  --chr 1-22 \
  --max-alleles 2 \
  --rm-dup exclude-mismatch \
  --remove deg2_hg38.king.cutoff.out.id \
  --set-all-var-ids '@:#:$1:$2' \
  --new-id-max-allele-len 487 \
  --make-bed \
  --out 1kg_hg38_autosomes

# 3. Baixar BED do Agilent SureSelect Human All Exon v2 via UCSC (hg38)
echo "[3/5] Downloading SureSelect v2 exome regions from UCSC..."
if [ ! -f Agilent_SureSelect_V2_Regions.bed ] || [ "$(wc -l < Agilent_SureSelect_V2_Regions.bed 2>/dev/null || echo 0)" -lt 1000 ]; then
  if [ -f Agilent_SureSelect_V2_Regions.bed ]; then
    echo "  Existing BED has <1000 lines (corrupt). Re-downloading..."
    rm -f Agilent_SureSelect_V2_Regions.bed
  fi
  if ! command -v bigBedToBed &>/dev/null; then
    echo "  bigBedToBed not found. Installing via micromamba..."
    micromamba install -y -c bioconda ucsc-bigbedtobed 2>/dev/null || \
      micromamba run -n base micromamba install -y -c bioconda ucsc-bigbedtobed 2>/dev/null || true
    hash -r
  fi
  if ! command -v bigBedToBed &>/dev/null; then
    echo "  ERROR: bigBedToBed not available. Install it with:"
    echo "    micromamba install -c bioconda ucsc-bigbedtobed"
    echo "  Then re-run this script."
    exit 1
  fi
  set +e
  wget -q "https://hgdownload.soe.ucsc.edu/gbdb/hg38/exomeProbesets/S30409818_Regions.bb" \
    -O S30409818_Regions.bb 2>/dev/null
  WGET_OK=$?
  if [ $WGET_OK -ne 0 ] || [ ! -s S30409818_Regions.bb ]; then
    set -e
    echo "  ERROR: download failed (exit code $WGET_OK)"
    echo "  Try manually: wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/exomeProbesets/S30409818_Regions.bb"
    rm -f S30409818_Regions.bb
    exit 1
  fi
  bigBedToBed S30409818_Regions.bb Agilent_SureSelect_V2_Regions.bed 2>/dev/null
  BB_OK=$?
  rm -f S30409818_Regions.bb
  set -e
  if [ $BB_OK -ne 0 ] || [ ! -f Agilent_SureSelect_V2_Regions.bed ]; then
    echo "  ERROR: bigBedToBed conversion failed"
    echo "  Try manually: bigBedToBed S30409818_Regions.bb Agilent_SureSelect_V2_Regions.bed"
    exit 1
  fi
  echo "  OK ($(wc -l < Agilent_SureSelect_V2_Regions.bed) regions)"
else
  echo "  Already exists"
fi

# 4. Filtrar 1KG para regiões exônicas
echo "[4/5] Filtering 1KG to exome regions..."
plink2 \
  --bfile 1kg_hg38_autosomes \
  --extract bed1 Agilent_SureSelect_V2_Regions.bed \
  --make-bed \
  --out 1kg_hg38_exome

# 5. Limpar PGEN brutos (só mantém BED processado)
echo "[5/5] Cleaning up raw PGEN..."
rm -f all_hg38.pgen all_hg38.pvar
echo "  OK"

echo ""
echo "=== Setup complete ==="
echo "1000G autosomes: 1kg/1kg_hg38_autosomes.{bed,bim,fam} ($(wc -l < 1kg_hg38_autosomes.bim) variants)"
echo "1000G exome:     1kg/1kg_hg38_exome.{bed,bim,fam} ($(wc -l < 1kg_hg38_exome.bim) variants)"
echo "Samples:          $(wc -l < 1kg_hg38_exome.fam)"
echo "Populations:      column 5 in 1kg/all_hg38.psam (SuperPopulation)"
