# Sperm miRNA signatures associated with male fertility status

This repository contains the analysis code and reproducibility information for a dry-lab MSc thesis investigating sperm microRNA (miRNA) expression, machine-learning panel development, and biological interpretation in relation to male fertility status.

## Study overview

The supplied pilot dataset contained sperm small RNA-sequencing counts and accompanying metadata for 75 male participants recruited in Singapore:

- 43 Fertile samples
- 17 Borderline samples
- 15 Infertile samples

The 58 Fertile and Infertile samples formed the labelled cohort used for differential-expression analysis and supervised model development. Borderline samples were excluded from model fitting and threshold selection and were projected onto the fixed final score only after model development.

## Analysis workflow

1. Match the miRNA count matrix to sample metadata.
2. Remove sparsely detected miRNAs using a label-independent expression-prevalence filter.
3. Apply trimmed mean of M-values (TMM) normalisation and limma-voom modelling.
4. Identify Fertile-versus-Infertile candidates while adjusting the differential-expression model for age and body mass index.
5. Compare Elastic Net, random forest, ridge logistic regression, and decision tree models using repeated nested stratified five-fold cross-validation.
6. Tune the selected model and reduce the candidate panel using selection stability, leave-one-miRNA-out ablation, panel-size comparison, and repeated seed checks.
7. Calculate internal out-of-fold predictions and an exploratory 0-100 score from the locked panel.
8. Query database-supported miRNA-target interactions with multiMiR and test KEGG pathway over-representation with clusterProfiler.

## Software

The analyses were performed in R 4.5.2. Principal packages included:

- edgeR 4.6.3
- limma 3.64.3
- caret 7.0.1
- glmnet 5.0
- randomForest 4.7-1.2
- rpart 4.1.27
- pROC 1.19.0.1
- multiMiR 1.30.0
- clusterProfiler 4.16.0
- ggplot2 4.0.3

The complete computational environment is recorded in `sessionInfo.txt`.

## Random seeds and resampling

The following seeds were fixed before the corresponding analyses. Consecutive ranges are inclusive.

| Analysis stage | Fold design | Recorded seeds |
|---|---|---|
| Nested model-class comparison | 15 repeated stratified outer five-fold allocations | Outer seeds 1001-1015; inner-fold seeds were generated deterministically as `100000 + 100 x repeat index + outer-fold index` |
| Elastic Net full-data parameter comparison | 15 repeated stratified five-fold allocations | Seeds 50123-50137 |
| Elastic Net stability, ablation and panel-size comparison | 50 stratified five-fold allocations | Seeds 20260628-20260677 |
| Four-miRNA repeated out-of-fold score analysis | 50 stratified five-fold allocations | Seeds 20260724-20260773 |
| Random-forest backward-ablation sensitivity analysis | 50 stratified five-fold allocations | Seeds 8001-8050 |

Within-fold model-fitting seeds were derived deterministically from the corresponding allocation seed and fold index, as specified in the analysis scripts. Seed values are integer identifiers and should not be interpreted as calendar dates.

## Recommended repository structure

```text
.
├── README.md
├── sessionInfo.txt
├── code/
│   ├── 01_data_preparation.R
│   ├── 02_differential_expression.R
│   ├── 03_model_comparison.R
│   ├── 04_elastic_net_optimisation.R
│   ├── 05_panel_reduction.R
│   ├── 06_score_and_roc.R
│   ├── 07_target_and_pathway_analysis.R
│   └── 08_generate_figures.R
└── outputs/
    └── non-identifiable summary tables and figures
```

The script names above provide a recommended upload order. They may be adapted to match the final script filenames, provided the execution order remains clear.

## Data availability and privacy

Participant-level count data, metadata, sample identifiers, and other potentially sensitive study materials are **not included** in this repository. These data were supplied for the MSc project and should not be uploaded to a public repository without explicit permission from the data owner and project supervisors.

The repository should contain only analysis code, software information, and non-identifiable summary outputs. Users with authorised access to the source data must place the required input files locally and update the input paths before running the scripts.

## Reproducibility notes

- Random seeds and cross-validation folds should be fixed within the relevant scripts.
- Preprocessing must be estimated within training partitions whenever a model is evaluated on held-out observations.
- The Borderline group must remain excluded from supervised model development.
- Reported post-selection out-of-fold performance is an internal development estimate, not independent validation.
- Any future update to the panel, coefficients, preprocessing, or cut-off should be recorded as a new model version.

## Citation

If this repository is cited, use the final thesis citation and the permanent repository URL. No participant-level data are distributed with the code.
