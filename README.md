# SNNB_2026_Fine_Mapping

Validação de variantes codantes fine-mapped do câncer de mama no **MVP (Million Veteran Program)** em coorte independente de exoma (**phs000822**).

> **Nota (nome das pastas):** o repositório no GitHub chama-se `SNNB_2026_Fine_Mapping_reply`; a pasta de trabalho no cluster foi renomeada para `SNNB_2026_Fine_Mapping_scripts`. São o mesmo projeto.

## Visão geral

| | Discovery (fine-mapping) | Validação |
|---|---|---|
| **Coorte** | MVP — Verma et al. 2024, *Science* (PMID 39024449) | Exoma phs000822 (466 casos early-onset, ≤35 anos) |
| **Dados** | `GCST90479802.tsv.gz` + `fine_mapping_GCST90479802.xlsx` (GWAS Catalog) | `gwas_dbgap_sem_filtragem+joint_call+vep.tsv` (summary statistics com VEP) |
| **Objetivo** | Credible sets com PIP | Look-up por rsID, direção de efeito, MAGMA |

O GWAS do exoma phs000822 (summary statistics) é gerado pelo pipeline
**`pipeline_gwas/`** (adaptado do **quali-workflow**,
`github.com/MatheusABomfim/quali-workflow`): VCF Sarek → scripts 01–05 →
`gwas_ssf.tsv` → `05.5_prepara_exoma.R` → `exoma_sumstats.txt`.

## Fluxo

```
Pipeline GWAS do exoma (pipeline_gwas/)
  VCF Sarek (joint-calling, VEP+AlphaMissense)
    ├─► 01/01.1 setup+validate 1KG      (SKIP: 1KG já configurado)
    ├─► 02_subsample_debug              (debug chr22, 50+50)
    ├─► 03_pca_matching                 (PCA + knee plot)
    ├─► 04_gwas_assoc (04.3 produção)   → SRR_gwas.assoc + target
    ├─► 05_gwas_ssf                     → gwas_ssf.tsv + YAML
    └─► 05.5_prepara_exoma              → exoma_sumstats.txt
                                              │
                                              ▼
MVP fine-mapping (GCST90479802)
  ├─► 01_filtrar_mvp.R        sinais de câncer de mama (Biological Mother)
  │         └─► bc_signal.csv
  │
  ▼
02_lookup_exoma.R  ──►  match por rsID nas sumstats do exoma
  │                     classificação: Replicada / Direcao_consistente / Discordante / Nao_genotipada
  │                     └─► lookup_completo.csv
  ├─► 03_analise_ancestralidade.R   (AFR / AMR / EUR)
  ├─► 04_tabelas_finais.R           tabelas 1-4 do artigo
  ├─► 05_locus_plots.R              locus plots dos top loci
  └─► 06_scatter_plot.R             scatter direção do efeito MVP × exoma

MAGMA (gene-level)  →  colaborador (etapa futura), output em resultados/magma/
```

## Estrutura

```
pipeline_gwas/    # geração do GWAS summary do exoma (quali-workflow adaptado)
dados/
├── mvp/          # GCST90479802.tsv.gz, fine_mapping_GCST90479802.xlsx (no cluster)
├── exoma/        # exoma_sumstats.txt (output do 05.5), genes.annot
└── referencias/  # gene_loc.txt, LD reference (MAGMA futuro)
scripts/          # validação MVP×exoma: 01_filtrar_mvp.R ... 06_scatter_plot.R
resultados/
├── tabelas/      # tabela1_consistentes, tabela2_discordantes, tabela3_novos_candidatos, tabela4_magma_genes
├── figuras/      # locus_plot_*.png, scatter_direcao.png
└── magma/        # exoma.genes.out (colaborador)
artigo/           # manuscript.md
```

## Uso

### 1. Gerar o GWAS summary do exoma (no cluster)

```bash
# a) produção (1KG já configurado → pular 01/01.1):
qsub pipeline_gwas/04.3_assoc_prod.pbs   # → gwas/gwas_prod/SRR_gwas.assoc
qsub pipeline_gwas/05_gwas_ssf.pbs       # → gwas_ssf.tsv + YAML

# b) converter para o formato do MVP lookup:
Rscript pipeline_gwas/05.5_prepara_exoma.R   # → dados/exoma/exoma_sumstats.txt
```

Ver `pipeline_gwas/README.md` para detalhes.

### 2. Validação MVP × exoma

```bash
# no cluster, dentro de scripts/ (os PBS usam caminhos relativos a esse dir)
cd scripts
qsub 01_filtrar_mvp.pbs
qsub 02_lookup_exoma.pbs   # depende do 01 (bc_signal.csv)
qsub 03_analise_ancestralidade.pbs   # depende do 02
qsub 04_tabelas_finais.pbs           # depende do 02
qsub 05_locus_plots.pbs              # depende do 02
qsub 06_scatter_plot.pbs             # depende do 02
```

Ou diretamente (sem PBS), no cluster:

```bash
Rscript scripts/01_filtrar_mvp.R
Rscript scripts/02_lookup_exoma.R
Rscript scripts/03_analise_ancestralidade.R
Rscript scripts/04_tabelas_finais.R
Rscript scripts/05_locus_plots.R
Rscript scripts/06_scatter_plot.R
```

Os `.pbs` de `scripts/` usam o env `quali` e rodam a partir de `scripts/`
(`cd "$PBS_O_WORKDIR"`), então devem ser submetidos **dentro de `scripts/`**.
A ordem é importante: 02–06 consomem `resultados/tabelas/lookup_completo.csv`
(output do 02).

Os caminhos dos dados e o mapeamento de colunas ficam configurados no bloco
`CONFIG` no topo de cada script — preencha conforme os headers reais dos arquivos.

## Notas

- Dados brutos (xlsx/tsv/vcf/PLINK) não são versionados — ficam no cluster
  (`/storage4/matheusbomfim/SNNB_2026_Fine_mapping/` = dados do MVP + dbGaP;
  outputs do GWAS em `SNNB_2026_Fine_Mapping_scripts/` conforme os CONFIG dos scripts).
- A pasta do repo no cluster é `/storage4/matheusbomfim/SNNB_2026_Fine_Mapping_scripts`.
- `pipeline_gwas/` reutiliza o painel 1KG já configurado pelo quali-workflow
  (`/storage4/matheusbomfim/quali/1kg/1kg_hg38_exome`) e a instalação LDSC.
- MAGMA gene analysis será desenvolvido por colaborador; `04_tabelas_finais.R`
  já lê `resultados/magma/exoma.genes.out` quando disponível.
