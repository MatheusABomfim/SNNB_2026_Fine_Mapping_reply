#!/bin/bash
# 01.1_validate_setup.sh — Valida saída do 01_setup_1kg.sh
# Uso: bash scripts/01.1_validate_setup.sh
set -euo pipefail

BASE="/storage4/matheusbomfim/quali/1kg"
PASS=0
FAIL=0

check() {
  local desc="$1"
  shift
  if "$@" 2>/dev/null; then
    echo "  PASS: $desc"
    PASS=$((PASS + 1))
  else
    echo "  FAIL: $desc"
    FAIL=$((FAIL + 1))
  fi
}

echo "=== Validating 1KG setup ==="
echo ""

echo "[Exome files (used in downstream analyses)]"
check "1kg_hg38_exome.bed exists"      test -f "$BASE/1kg_hg38_exome.bed"
check "1kg_hg38_exome.bim exists"      test -f "$BASE/1kg_hg38_exome.bim"
check "1kg_hg38_exome.fam exists"      test -f "$BASE/1kg_hg38_exome.fam"
check "Agilent_SureSelect_V2_Regions.bed exists" test -f "$BASE/Agilent_SureSelect_V2_Regions.bed"
check "all_hg38.psam exists"           test -f "$BASE/all_hg38.psam"

echo ""
echo "[Cleanup (must NOT exist)]"
check "all_hg38.pgen removed"          test ! -f "$BASE/all_hg38.pgen"
check "all_hg38.pvar removed"          test ! -f "$BASE/all_hg38.pvar"

echo ""
echo "[Counts]"
N_WGS=$(wc -l < "$BASE/1kg_hg38_autosomes.bim" 2>/dev/null || echo 0)
N_EXOME=$(wc -l < "$BASE/1kg_hg38_exome.bim" 2>/dev/null || echo 0)
N_SAM=$(wc -l < "$BASE/1kg_hg38_exome.fam" 2>/dev/null || echo 0)
N_REMOVED=$((N_WGS - N_EXOME))
PCT_REMOVED=$(awk "BEGIN {printf \"%.2f\", ($N_REMOVED / $N_WGS) * 100}")
echo "  WGS variants:    $N_WGS"
echo "  Exome variants:  $N_EXOME"
echo "  Removed:         $N_REMOVED ($PCT_REMOVED%)"
echo "  Samples:         $N_SAM"
if [ "$N_EXOME" -gt 300000 ] && [ "$N_EXOME" -lt 2000000 ]; then
  echo "  PASS: exome variant count in expected range (300K-2M)"
  PASS=$((PASS + 1))
else
  echo "  FAIL: exome variant count ($N_EXOME) outside expected range (300K-2M)"
  FAIL=$((FAIL + 1))
fi
if [ "$N_SAM" -gt 2400 ] && [ "$N_SAM" -lt 2600 ]; then
  echo "  PASS: sample count in expected range (2400-2600)"
  PASS=$((PASS + 1))
else
  echo "  FAIL: sample count ($N_SAM) outside expected range (2400-2600)"
  FAIL=$((FAIL + 1))
fi

echo ""
echo "[Populations]"
cut -f5 "$BASE/all_hg38.psam" 2>/dev/null | sort | uniq -c | sort -rn | while read -r n pop; do
  echo "  $pop: $n"
done

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
if [ "$FAIL" -gt 0 ]; then
  exit 1
fi
