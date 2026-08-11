# pipeline_gwas — Geração do GWAS Summary do exoma phs000822

Pipeline adaptado do **quali-workflow** (`github.com/MatheusABomfim/quali-workflow`)
para gerar as summary statistics do exoma (câncer de mama, phs000822) que
alimentam a validação das variantes fine-mapped do MVP (`scripts/` neste repo).

> **Nota:** a pasta do repo no cluster é `/storage4/matheusbomfim/SNNB_2026_Fine_Mapping_scripts`
> (renomeada de `SNNB_2026_Fine_Mapping_reply`). Os inputs/dados do MVP ficam em
> `/storage4/matheusbomfim/SNNB_2026_Fine_mapping/`.

## Fluxo

```
VCF Sarek (joint-calling, VEP+AlphaMissense)
  │
  ├─► 01_setup_1kg.sh         SKIP (1KG já configurado) / 01.1 validação
  ├─► 02_subsample_debug.sh   debug (chr22, 50 cases + 50 1KG)
  ├─► 03_pca_matching.R       PCA + knee plot (diagnóstico)
  ├─► 04_gwas_assoc.R         GWAS: QC → merge → PCA → matching → assoc PLINK2
  │      (produção: 04.3_assoc_prod.pbs → --firth --pca-match --match-k 1 --n-pcs 3)
  ├─► 04.4_hits_gwas.R        → gwas_hits_significativos.csv / _sugestivos.csv
  │                             (hits + anotação VEP por transcrito)
  ├─► 05_gwas_ssf.R           → gwas_ssf.tsv + YAML (GWAS-SSF v1.1)
  ├─► 05.5_prepara_exoma.R    → exoma_sumstats.txt (formato do 02_lookup_exoma.R)
  ├─► 06_locuszoom.R          → locuszoom_*.png (por locus) + gwas_localzoom.tsv.gz/.tbi
  │                             (para LocalZoom) + Manhattan hits/sugestivos
  └─► 07_circos_manhattan.R   → circos_plot.{pdf,png} + fig_circos_manhattan.{pdf,png}
```

## Inputs

| Input | Path |
|---|---|
| VCF cases (Sarek) | `/storage4/matheusbomfim/programas/sarek/results_CA_Mama_DBGAP_New/variant_calling/haplotypecaller/joint_variant_calling/filtered/annotated/joint_germline.filtered.norm.AlphaMissense.vep.vcf.gz` |
| Painel 1KG hg38 (pronto) | `/storage4/matheusbomfim/quali/1kg/1kg_hg38_exome` |

## Outputs

| Etapa | Output |
|---|---|
| 04 | `gwas/<timestamp>/SRR_gwas.assoc` (tabela limpa), plots, `gwas_report.txt` |
| 04.3 | `gwas/gwas_prod/SRR_gwas.assoc` + `target.{bed,bim,fam}` — escrita direta via `--outdir gwas_prod` (sem cópia manual) |
| 04.4 | `resultados/tabelas/gwas_hits_significativos.csv` + `gwas_hits_sugestivos.csv` (1 linha por variante×transcrito VEP) |
| 05 | `gwas/gwas_prod/gwas_ssf.tsv` + `-meta.yaml` |
| 05.5 | `gwas/gwas_prod/exoma_sumstats.txt` → copiar para `dados/exoma/exoma_sumstats.txt` |
| 06 | `resultados/figuras/locuszoom_chr{chr}_{pos}.png`, `manhattan_hits_significativos.png`, `manhattan_hits_sugestivos.png`; `resultados/tabelas/gwas_localzoom.tsv.gz` + `.tbi` |
| 07 | `resultados/figuras/circos_plot.{pdf,png}` (circos), `fig_circos_manhattan.{pdf,png}` (Manhattan completo) |

## Ordem de execução (cluster)

```bash
# 1KG já existe — 01/01.1 podem ser pulados
# Debug (opcional):
qsub 02_subsample_debug.pbs
qsub 04_gwas_assoc_debug.pbs   # via 04_gwas_assoc.R --debug

# PCA diagnóstico (opcional):
qsub 03_pca_matching.pbs

# Produção:
qsub 04.3_assoc_prod.pbs        # → gwas/gwas_prod/SRR_gwas.assoc
qsub 05_gwas_ssf.pbs            # → gwas_ssf.tsv + YAML

# Hits + anotação VEP (env quali):
qsub 04.4_hits_gwas.pbs         # → resultados/tabelas/gwas_hits_{significativos,sugestivos}.csv

# Figuras (env locuszoom):
qsub 06_locuszoom.pbs           # → locuszoom_*.png + LocalZoom (.tsv.gz/.tbi) + Manhattan
qsub 07_circos_manhattan.pbs    # → circos_plot.* + fig_circos_manhattan.*

# Converter para o formato do MVP lookup:
Rscript 05.5_prepara_exoma.R    # → exoma_sumstats.txt

# Em seguida, no workflow de validação:
Rscript scripts/02_lookup_exoma.R
```

## Notas

- **LDSC** (step 12 do 04) é opcional e será pulado graciosamente se as LD
  scores/`ldsc.py` não estiverem disponíveis (aponta para a instalação do
  quali-workflow em `/storage4/matheusbomfim/tools/ldsc`).
- **rsID** no `exoma_sumstats.txt`: anotado a partir do ID column do VCF Sarek
  via `bcftools query` (script 05.5). Variantes sem rs ficam com `RSID=""`.
- Ambientes Micromamba/PBS: mantidos do quali-workflow
  (`/storage2/matheusbomfim/projects/micromamba`, env `quali`).
- Parâmetros de produção finais (match grid 04.2, ponto ótimo k1 pc3):
  `--firth --pca-match --match-k 1 --n-pcs 3` (λ GC 1.2417).
- **Env `locuszoom`** (06/07): `micromamba create -n locuszoom -c conda-forge -c bioconda r-locuszoomr r-circlize r-ggplot2 r-ggrepel tabix` e instalar o `bioconductor-ensembldb`/`bioconductor-ensdb.hsapiens.v86` (necessário pelo `locuszoomr`). O 07 também precisa de `r-circlize` e (opcionalmente) `r-gridgraphics`.
- **Gene canônico dos rótulos**: o CSQ do VEP/AlphaMissense do Sarek não grava o
  campo `CANONICAL`, então o gene por locus é derivado do próprio CSV de hits:
  transcrito de **maior IMPACT** (HIGH > MODERATE > LOW > MODIFIER), desempate
  pelo `transcript_index` menor (ordem do VEP = consequência mais grave primeiro).
- **LocalZoom**: o `06_locuszoom.R` exporta `gwas_localzoom.tsv.gz` + `.tbi`
  (header simples, ordenado, `tabix --skip-lines 1`). No site
  `https://statgen.github.io/localzoom/`, arrastar os DOIS arquivos e mapear as
  colunas (chrom, pos, ref, alt, id, beta, p, logp — indicadas no log).
