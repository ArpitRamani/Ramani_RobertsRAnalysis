# Ramani_RobertsRAnalysis

Differential-expression analysis of **protein nitration** (nitro-tyrosine / nitro-tryptophan)
in the **McEachin et al. C9orf72 ALS spinal-cord** TMTpro dataset — the R/statistics half of
the lab's nitration project (applying an in-house CSF nitration workflow to the McEachin spinal-cord dataset).

Upstream (mass-spec → site matrix) is done in **FragPipe → TMT-Integrator**. This repo starts
from the TMT-Integrator single-site abundance matrix and produces per-contrast DE tables.

---

## Pipeline

```
TMT-Integrator single-site abundance  (sites × 36 channels)
        │
        ▼   drop pooled reference channel 126  → 34 patient samples
   TAMPOR (noGIS)            ← batch + sample-loading normalization (Seyfried-lab median polish)
        │                       (requires a site be present in BOTH plexes → restricts to the
        │                        cross-batch-detectable core, ~124 sites)
        ▼
   missingness filter        ← keep sites ≥50% present
        │
        ▼
   missForest imputation     ← random-forest fill of remaining gaps
        │
        ▼
   limma                     ← moderated t-stats per site, 3 contrasts, BH-FDR
        │
        ▼
   DE_<contrast>.csv         ← logFC, p, FDR  (volcanoes & GO done downstream by user)
```

### Design (from `data/metadata/original_metaData.csv`)
34 patients across 2 TMTpro-18 plexes (b01/b02), **balanced across plexes**:

| Group   | n  | Meaning                                   |
|---------|----|-------------------------------------------|
| c9ALS   | 18 | C9orf72-expansion ALS, untreated          |
| ASO     | 6  | C9 ALS treated with BIIB078 antisense drug |
| Control | 10 | non-ALS controls                          |

Channel **126** in each plex is the pooled reference (excluded from patient DE).

### Contrasts
- `c9ALS_vs_Control` — disease effect
- `ASO_vs_Control` — does treatment normalize nitration?
- `c9ALS_vs_ASO` — treatment effect

---

## Two normalization modes (important)

| Mode | Sites analyzed | When it's used |
|------|----------------|----------------|
| **TAMPOR (noGIS)** | ~124 | default; rigorous cross-batch normalization, but only keeps sites detected in **both** plexes |
| **log2 + median-center (fallback)** | ~419 | auto-used if TAMPOR deps are missing; keeps the sparse one-plex sites, adds batch as a limma covariate |

The gap is real and methodological: a sparse PTM is often sampled in only one plex, and TAMPOR
can't compute a cross-batch ratio for those. TAMPOR = robust core; fallback = fuller but noisier.
Choose per question; both are defensible.

---

## Setup

```r
# R >= 4.x. One-time install:
install.packages(c("missForest","randomForest","ggplot2","ggrepel","ggpubr","pheatmap",
                   "png","patchwork","doParallel","foreach","snow","doSNOW","data.table"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","vsn"))
```
> - `snow` is required for TAMPOR's parallel backend — without it the pipeline falls back to median-centering.
> - `pheatmap` powers the heatmaps; `png`/`patchwork` assemble the network+enrichment figures.
> - `run_STRING_PPI.R` needs only internet access (it calls the STRING REST API via `download.file`).

## Run

**1. Differential expression** — pick a mode (each writes to its own results folder):

```bash
Rscript scripts/run_McEachin_DE_paper.R      # PRIMARY: manuscript-faithful (log2, no TAMPOR);
                                             #          limma + Welch concordance. 419 sites.
Rscript scripts/run_McEachin_DE_TAMPOR.R     # robustness check: TAMPOR core. 124 sites.
Rscript scripts/run_McEachin_DE_fallback.R   # log2 + median-center, batch in model. 419 sites.
```
Each writes per-contrast `DE_*.csv` (logFC, p, FDR; the paper mode adds `welch_P`/`welch_adjP`),
`DE_summary.csv`, plus `cleanDat_normalized_imputed.tsv` + `samples.csv` for the steps below.

**2. Downstream** (default to the paper results; pass another folder as arg 1):

```bash
Rscript scripts/make_plots.R           [results_dir]   # volcano (.pdf/.png) + heatmap per contrast
Rscript scripts/run_STRING_PPI.R       [results_dir]   # STRING network image + edges + enrichment (API)
Rscript scripts/make_PPI_figures.R     [results_dir]   # network + enrichment-bar figures (pathways, GO)
Rscript scripts/compare_limma_welch.R  [results_dir]   # limma vs Welch concordance table + plots
```
Outputs land in `<results_dir>/plots/`, `/PPI/`, and `/comparison/`.

**3. Slide deck** (optional — collects every figure into a PowerPoint to share):

```bash
Rscript scripts/export_slide_pngs.R  [results_dir]   # render all figures to <results_dir>/slides_png/
python3  scripts/build_pptx.py       [results_dir]   # -> <results_dir>/McEachin_nitration_figures.pptx
```
One figure per slide, title on top, an editable notes box at the bottom. Needs `python-pptx`
(`pip install python-pptx`). The `.pptx` and PNGs live under `results/` (gitignored), not the repo.

> **limma vs Welch.** The paper driver runs both tests on the same sites (`concordance=TRUE`),
> so each `DE_*.csv` carries limma (`P.Value`/`adj.P.Val`) **and** Welch (`welch_P`/`welch_adjP`).
> `compare_limma_welch.R` quantifies the difference: `limma_vs_welch_concordance.csv` (per-contrast
> Spearman/Pearson agreement, sig-count overlap, % where limma is more significant), per-contrast
> p-value scatters (`*_p_concordance.pdf`), a method-count bar (`method_significance_counts.pdf`),
> and `variance_moderation.pdf` — the empirical-Bayes shrinkage of each site's variance toward the
> prior, which is *why* limma and Welch diverge at this sample size.

> **Enrichment / GO:** STRING's enrichment (`/PPI/*_enrichment.tsv`) already covers GO Biological
> Process / Molecular Function / Cellular Component, pathways (KEGG, Reactome, WikiPathways), and HPO
> — filter by the `category` column. `make_PPI_figures.R` renders the top terms as −log10(FDR) bars
> beside each STRING network, so there is no separate GO step.

Both DE drivers share `R/nitro_DE_functions.R`; edit the `cfg` block at the top of any driver.

---

## Repo layout

```
Ramani_RobertsRAnalysis/
├── README.md
├── scripts/
│   ├── run_McEachin_DE_paper.R     # PRIMARY DE driver (manuscript-faithful)
│   ├── run_McEachin_DE_TAMPOR.R    # TAMPOR-mode driver
│   ├── run_McEachin_DE_fallback.R  # fallback-mode driver
│   ├── make_plots.R                # volcano + heatmap
│   ├── run_STRING_PPI.R            # STRING PPI networks + enrichment (API)
│   ├── make_PPI_figures.R          # network + enrichment-bar combined figures
│   └── compare_limma_welch.R       # limma vs Welch concordance + variance-moderation plots
├── R/
│   ├── nitro_DE_functions.R        # shared DE pipeline (load→normalize→impute→test)
│   ├── plotting_functions.R        # volcano, heatmap, enrichment bars, network combiner
│   ├── enrichment_functions.R      # STRING API helper
│   └── vendor/TAMPOR.R             # official TAMPOR (edammer/TAMPOR), vendored
├── data/
│   ├── McEachin_single-site_nitro_snapshot.tsv   # TMT-Integrator input (snapshot)
│   └── metadata/
│       ├── original_metaData.csv   # full clinical metadata (used for grouping)
│       └── numeric_metaData.csv    # slim Sample/Batch/Disease/Sex/Age
└── results/
    └── McEachin_DE_paper/          # DE outputs + plots/ GO/ PPI/ subfolders
```

## Provenance / upstream parameters
TMT-Integrator config that produced the input matrix: `min_snr=180`, `min_resolution=0`,
`min_site_prob=0.1`, `min_purity=0.5`, `min_pep_prob=0.9`, `ms1_int=true`, `add_Ref=1`
(matches the reference in-house run). The `min_snr=180` value is critical — `1000` collapses the result.

TAMPOR: Dammer et al., Seyfried lab — https://github.com/edammer/TAMPOR
