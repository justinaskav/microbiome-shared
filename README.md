# Microbiome Shared Resources

R functions, Quarto assets and environment scripts shared across the microbiome projects.
Projects include this repository as a git submodule (`shared/`) and `source()` the files they
need. The Nextflow CLI scripts (`calc_alpha_diversity.R`, `calc_beta_diversity.R`, ...) live in
the `microbiome-workflows` repository, not here.

## Layout

```
R/
├── alpha_diversity/        LMM/Kruskal summaries, sensitivity analyses, panels, diagnostics,
│                           per-subject deltas, rarefaction curves
├── beta_diversity/         PERMANOVA (sample- and patient-level), dbRDA + variable selection,
│                           ordination panels, patient aggregation
├── differential_abundance/ MaAsLin3 result handling, forest/volcano/EMM plots, interaction tables
├── metadata/               metadata correlation heatmap and network
├── ml/                     mikropml tasks, metrics/ROC plots, production models (save/load/score)
├── network/                NetCoMi networks and ggraph plotting
├── taxonomic/              aggregation, summaries, composition/heatmap/box-jitter plots
└── util/                   significance stars, figure saving, table helpers
quarto/                     CSL styles, reference .docx/.pptx, logo, SCSS
env/                        conda env and package/QIIME 2 setup scripts
```

## Usage

```r
shared <- "shared/R"
source(file.path(shared, "util/classify_significance.R"))
source(file.path(shared, "taxonomic/calc_taxonomic_summary.R"))
source(file.path(shared, "alpha_diversity/calc_alpha_diversity_summary.R"))
```

Files are plain scripts, not a package, so a function's helpers must be sourced too. Each
file's header names any other shared file it needs ("Requires ..."). The most common pairs:

| Function(s)                                   | Also source                                      |
|-----------------------------------------------|--------------------------------------------------|
| `calc_alpha_diversity_summary()`              | `taxonomic/calc_taxonomic_summary.R`, `util/classify_significance.R` |
| `plot_alpha_diversity_summary_lmm()`, panels  | `alpha_diversity/calc_alpha_diversity_summary.R`, `alpha_plot_style.R` |
| `beta_permanova_split()`, `run_dbrda_analysis()`, `plot_rda_summary()` | `beta_diversity/beta_diversity.R` (`tidy_vegan_anova()`) |
| `beta_permanova_split()`, `plot_taxonomic_heatmap(aggregate_by = )` | `beta_diversity/aggregate_to_patient.R` |
| `dbrda_variable_selection()`                  | `beta_diversity/run_dbrda_analysis.R`            |
| MaAsLin3 plots, `extract_maaslin3_stats()`    | `differential_abundance/prepare_maaslin3_data.R` |
| `score_ml_model()`                            | `ml/prepare_ml_data.R`, `ml/prepare_ml_covariates.R` |

Functions use `pkg::fn()` and the native pipe, so they do not require any package to be attached.

## Requirements

R ≥ 4.4 (base `%||%`), ggplot2 ≥ 3.5 (tested with 4.0), plus the packages listed in
`env/r-microbiome.yml` and `env/setup_packages.R` (mia, vegan, lmerTest, emmeans, maaslin3,
mikropml, NetCoMi, gt, ...).
