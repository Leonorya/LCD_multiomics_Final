# Multi-omics responses to a low-calorie diet in type 2 diabetes

Analysis and plotting code for the longitudinal multi-omics study of a 24-week low-calorie diet intervention in adults with
type 2 diabetes.

## Requirements

Analyses were performed in R 4.4.2. Packages used by the scripts:

- Mixed models and statistics: lme4, lmerTest, vegan, permute, pROC, glmnet
- Trajectory clustering: navmix
- Pathway and network analysis: clusterProfiler, KEGGREST, org.Hs.eg.db, AnnotationDbi, STRINGdb, igraph
- High-dimensional mediation: HIMA, hdi, HDMT (the scripts call internal HIMA functions, so results depend on the installed
  HIMA version)
- Plotting: ggplot2, patchwork, ggrepel, ggalluvial, circlize, ggplotify, scales

Pairwise associations in the networks are partial repeated-measures correlations, an extension of repeated-measures
correlation (rmcorr) that also removes visit means (`partial_rmcorr()` in `analysis/utils.R`).

## Data

Deidentified individual participant data are available to qualified researchers upon reasonable request to the corresponding
author, subject to approval of a research proposal and a data access agreement. The scripts expect the input files under
`data/` (clinical measures, medication use, omics feature tables with annotations, pooled QC tables, species abundances, gene
counts and species-level KEGG annotations); the file names are listed at the top of each script.

## Running the code

Run the scripts from the repository root in the following order; results are written to `results/` and figures to `figures/`.

1. `analysis/00_preprocessing.R`
2. `analysis/Figure1_time_effects.R`
3. `analysis/Figure2_blood_trajectories_pathways.R`
4. `analysis/Figure3_fecal_diversity_trajectories.R`
5. `analysis/Figure4_cross_omics_network.R`
6. `analysis/Figure5_clinical_associations_mediation.R`
7. `analysis/Figure6_remission.R`
8. `analysis/Figure7_prediction.R`
9. `plotting/Figure1.R` to `plotting/Figure7.R`

Parallel steps use at most `N_CORES` workers (`analysis/utils.R`).

## Repository structure

```
analysis/
  utils.R                                      Shared settings and helper functions: transformation and baseline
                                               standardization, time model, effect-modification model, NAvMix clustering,
                                               KEGG over-representation analysis, partial repeated-measures correlation
  00_preprocessing.R                           Transformation and baseline standardization of the six omics layers;
                                               principal component analysis of study and pooled QC samples (Figure S1)
  Figure1_time_effects.R                       Baseline characteristics and clinical outcomes; time model; visit
                                               localization; effect modification by diet group, sex and pre-intervention
                                               metformin use
  Figure2_blood_trajectories_pathways.R        Trajectory clustering, KEGG analysis of clusters and STRING protein
                                               interactions (k-core) for the blood omics layers
  Figure3_fecal_diversity_trajectories.R       Alpha diversity and gene count; PERMANOVA and dispersion; distance from
                                               baseline; trajectory clustering and species-to-pathway mapping for the fecal
                                               omics layers; effect modification of trajectory-cluster scores
  Figure4_cross_omics_network.R                Partial repeated-measures correlation network of the top 50 time-responsive
                                               features per layer; Louvain modules; cross-layer connectivity; KEGG analysis
                                               of cross-layer modules
  Figure5_clinical_associations_mediation.R    Longitudinal associations with clinical traits; adiposity and glycemia
                                               feature sets (KEGG analysis and correlation networks); adjustment for
                                               time-varying BMI; mediation and serial mediation by trajectory-cluster
                                               scores; high-dimensional mediation (HIMA)
  Figure6_remission.R                          Remission GLMM; KEGG analysis of remission-associated features; agreement
                                               with the longitudinal HbA1c association; gut species by remission status
  Figure7_prediction.R                         Base and omics models with repeated nested cross-validation (LASSO);
                                               performance of participant-averaged predictions with 95% CIs, calibration
                                               and selection frequency; Williams' and DeLong tests; decision curve analysis
annotation/
  label_exclusions.csv                         Features putatively annotated as exogenous compounds, which are not labeled
                                               in the figures
plotting/
  plot_utils.R                                 Palettes, theme and shared panels (trajectories, heatmaps, KEGG, networks)
  Figure1.R                                    Figure 1 and Figure S1
  Figure2.R ... Figure7.R                      Figures 2-7
```

The plotting scripts draw the content of each panel; the final panel arrangement and styling were edited manually.
