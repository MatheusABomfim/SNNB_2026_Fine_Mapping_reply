#!/bin/bash
# setup_imputacao_env.sh
# Cria o env micromamba 'imputacao' com as ferramentas de imputação de baixa
# cobertura (GLIMPSE2 + SHAPEIT5 + utilitários) e baixa os mapas genéticos b38.
#
# Uso: bash scripts/setup_imputacao_env.sh
# Pré-requisito: micromamba no PATH (ou use o caminho absoluto abaixo).
set -euo pipefail

MICROMAMBA="/storage2/matheusbomfim/projects/micromamba/bin/micromamba"
ENV_NAME="imputacao"
MAPS_DIR="/storage4/matheusbomfim/quali/1kg/maps_b38"

echo "=== Setup env '$ENV_NAME' ==="

# 1. Criar o env (se não existir)
if "$MICROMAMBA" env list | grep -q "^\s*$ENV_NAME\s"; then
  echo "[1/3] Env '$ENV_NAME' já existe — pulando criação."
else
  echo "[1/3] Criando env '$ENV_NAME'..."
  "$MICROMAMBA" create -y -n "$ENV_NAME"
fi

# 2. Instalar ferramentas (conda-forge + bioconda)
echo "[2/3] Instalando shapeit5, glimpse-bio, bcftools, samtools, plink2, htslib..."
"$MICROMAMBA" install -y -n "$ENV_NAME" \
  -c conda-forge -c bioconda \
  shapeit5 \
  glimpse-bio \
  bcftools \
  samtools \
  plink2 \
  htslib

# 3. Baixar mapas genéticos GRCh38 (usados por SHAPEIT5 e GLIMPSE2)
echo "[3/3] Baixando mapas genéticos b38 para $MAPS_DIR"
mkdir -p "$MAPS_DIR"
for CHR in $(seq 1 22); do
  OUT="$MAPS_DIR/chr${CHR}.b38.gmap.gz"
  if [ -s "$OUT" ]; then
    echo "  chr${CHR}.b38.gmap.gz já existe"
  else
    echo "  baixando chr${CHR}.b38.gmap.gz..."
    wget -q -O "$OUT" \
      "https://raw.githubusercontent.com/odelaneau/shapeit5/master/resources/maps/genetic_maps.b38/chr${CHR}.b38.gmap.gz" \
      || echo "    AVISO: falha no download chr${CHR} (baixe manualmente para $OUT)"
  fi
done

echo ""
echo "=== OK ==="
echo "Env:      micromamba run -n $ENV_NAME <cmd>"
echo "Mapas:    $MAPS_DIR/chr<CHR>.b38.gmap.gz"
echo ""
echo "Se o download dos mapas falhou, baixe manualmente de:"
echo "  https://github.com/odelaneau/shapeit5/tree/master/resources/maps/genetic_maps.b38"
