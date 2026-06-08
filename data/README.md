# data/ — inputs (NOT committed)

Real data is **deliberately excluded** from this repo (see `.gitignore`): it contains
identifiable patient clinical metadata and unpublished study results. Place files here locally.

## Expected inputs

```
data/
├── McEachin_single-site_nitro_snapshot.tsv   # TMT-Integrator single-site abundance matrix
│                                             #   cols: Index, Gene, ProteinID, ... + b01_126 … b02_135N
└── metadata/
    └── original_metaData.csv                 # sample sheet (see minimal columns below)
```

## Minimal metadata the pipeline actually uses

The scripts only read **four** columns — you do **not** need (and should not publish) the
identifiable clinical fields (Case Number, Age at Death, Race, Braak/Thal, brain bank, etc.):

| column        | example | purpose                          |
|---------------|---------|----------------------------------|
| `Batch`       | b01     | plex; builds the channel id       |
| `Channel`     | 127C    | TMT channel; `Batch_Channel` = matrix column |
| `DiseaseGroup`| c9ALS   | grouping variable for DE          |
| `Sex`         | F       | optional covariate                |

A de-identified sheet with just these columns is sufficient to run everything.
