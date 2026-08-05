# SNNB_2026_Fine_Mapping_reply

Validação de variantes codantes fine-mapped do câncer de mama no **MVP (Million Veteran Program)** em coorte independente de exoma (**phs000822**).

## Visão geral

| | Discovery (fine-mapping) | Validação |
|---|---|---|
| **Coorte** | MVP — Verma et al. 2024, *Science* (PMID 39024449) | Exoma phs000822 (466 casos early-onset, ≤35 anos) |
| **Dados** | `GCST90479802.tsv.gz` + `fine_mapping_GCST90479802.xlsx` (GWAS Catalog) | `gwas_dbgap_sem_filtragem+joint_call+vep.tsv` (summary statistics com VEP) |
| **Objetivo** | Credible sets com PIP | Look-up por rsID, direção de efeito, MAGMA |

O GWAS do exoma phs000822 é gerado pelo pipeline **quali-workflow** (scripts 01–05:
`SRR_gwas.assoc` / `gwas_ssf.tsv` → summary statistics), o mesmo que produz as
sumstats do exoma usadas aqui. Ver: `github.com/MatheusABomfim/quali-workflow`.

## Fluxo

```
MVP fine-mapping (GCST90479802)
  ├─► 01_filtrar_mvp.R        coding variants (VEP) de câncer de mama
  │         └─► bc_coding.csv
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
dados/
├── mvp/          # GCST90479802.tsv.gz, fine_mapping_GCST90479802.xlsx (no cluster)
├── exoma/        # gwas_dbgap_sem_filtragem+joint_call+vep.tsv, genes.annot
└── referencias/  # gene_loc.txt, LD reference (MAGMA futuro)
scripts/          # 01_filtrar_mvp.R ... 06_scatter_plot.R
resultados/
├── tabelas/      # tabela1_replicadas, tabela2_discordantes, tabela3_novos_candidatos, tabela4_magma_genes
├── figuras/      # locus_plot_*.png, scatter_direcao.png
└── magma/        # exoma.genes.out (colaborador)
artigo/           # manuscript.md
```

## Uso

```bash
# no cluster, dentro do diretório do repo
Rscript scripts/01_filtrar_mvp.R
Rscript scripts/02_lookup_exoma.R
Rscript scripts/03_analise_ancestralidade.R
Rscript scripts/04_tabelas_finais.R
Rscript scripts/05_locus_plots.R
Rscript scripts/06_scatter_plot.R
```

Os caminhos dos dados e o mapeamento de colunas ficam configurados no bloco
`CONFIG` no topo de cada script — preencha conforme os headers reais dos arquivos.

## Notas

- Dados brutos (xlsx/tsv/vcf/PLINK) não são versionados — ficam no cluster
  (`/storage4/matheusbomfim/SNNB_2026_Fine_mapping/` e outputs do quali-workflow).
- MAGMA gene analysis será desenvolvido por colaborador; `04_tabelas_finais.R`
  já lê `resultados/magma/exoma.genes.out` quando disponível.
