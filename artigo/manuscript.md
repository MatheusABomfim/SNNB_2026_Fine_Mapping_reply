# Validation of fine-mapped coding variants for breast cancer in an independent exome cohort

**Autores**: [seu nome] et al.

## Abstract (~250 palavras)

**Background:** Genome-wide association studies have identified over 200 risk loci
for breast cancer, but the causal variants and target genes remain largely unknown.
Fine-mapping in the Million Veteran Program (MVP) identified 57,601 genetic signals
across 2,068 traits, including 6,318 with high-confidence (PIP > 0.95) single-variant
resolution. However, independent validation of coding variants within these credible
sets is lacking, particularly in admixed populations and early-onset disease.

**Methods:** We extracted fine-mapped signals for breast cancer from the MVP dataset
(GCST90479802; Verma et al. 2024) and filtered for coding variants (missense,
synonymous, stop-gained, splice-site). We compared direction of effect in an
independent whole-exome sequencing study of 466 early-onset breast cancer cases
(phs000822). Variants were classified as replicated (same direction, p<0.05),
direction-consistent (same direction, p≥0.05), or discordant (opposite direction).
Gene-level replication was assessed using MAGMA.

**Results:** Among [X] breast cancer fine-mapped signals in MVP, [Y] were coding
variants. Of these, [Z] were genotyped in the exome cohort. [W%] showed consistent
direction of effect, and [V] reached nominal significance (p<0.05). Discordant
variants comprised [U%]. Gene-level analysis identified [T] genes with p<0.01 in the
exome cohort, including known breast cancer susceptibility genes (ATM, CHEK2, BRCA2).

**Conclusion:** Independent validation of fine-mapped coding variants in an exome
cohort supports their causal role at breast cancer GWAS loci. Our results highlight
the value of cross-cohort validation in diverse and early-onset populations for
prioritizing putative causal variants.

## Introdução (~1 página)

- GWAS de câncer de mama: 200+ loci (BCAC, Michailidou 2017, Zhang 2020)
- Fine-mapping: MVP (Verma 2024) — 57,601 sinais, 6,318 com PIP>0.95
- Lacuna: maioria dos estudos foca em europeus; validação em miscigenados é rara
- Contribuição: validar coding variants do MVP em exoma independente de
  early-onset breast cancer em população americana miscigenada

## Métodos (~1.5 páginas)

- **MVP (GCST90479802)**:
  - 57,601 sinais fine-mapped com SuSiE
  - 3 populações: AFR, EUR, AMR (hispanic/latino)
  - PIP, beta, VEP annotation por variante
- **Exoma phs000822**:
  - 466 casos early-onset breast cancer (≤35 anos)
  - Whole-exome (Agilent SureSelect v2, Illumina HiSeq)
  - Summary statistics para variantes comuns (MAF>1%) geradas pelo pipeline
    quali-workflow (gwas_dbgap_sem_filtragem+joint_call+vep.tsv)
- **Análise**:
  - Filtrar MVP para câncer de mama + variantes codantes
  - Look-up por rsID no exoma
  - Classificar: replicada, direção consistente, discordante
  - MAGMA gene-based analysis no exoma

## Resultados (~2 páginas)

- Tabela 1: Coding variants replicadas (direção consistente)
- Tabela 2: Discordantes
- Tabela 3: Novos candidatos do exoma
- Tabela 4: Genes MAGMA
- Figura 1: Scatter plot direção do efeito
- Figura 2: Locus plots dos top loci

## Discussão (~1 página)

- Coding variants como causais em loci de GWAS (ex: ATM, BRCA2, CHEK2)
- Consistência entre populações AFR, AMR, EUR
- Limitações:
  - N pequeno do exoma (466) → poder limitado
  - Apenas variantes comuns (MAF>1%)
  - Possível heterogeneidade de fenótipos (MVP: geral, phs000822: early-onset)
- Implicações: fine-mapping + validação em exoma independente fortalece a
  priorização de variantes causais

## Referências (~30)

- Verma A, et al. Diversity and scale: Genetic architecture of 2068 traits in the VA Million Veteran Program. *Science*. 2024. PMID 39024449.
- Zhang H, et al. (BCAC). Genome-wide association study identifies 32 novel breast cancer susceptibility loci. *Nat Genet*. 2020.
- Michailidou K, et al. (BCAC). Association analysis identifies 65 new breast cancer risk loci. *Nature*. 2017.
- phs000822 (dbGaP) — Breast cancer in women with hereditary indication.
- Sun L, et al. (coding variant associations). *Nature*. 2022.
