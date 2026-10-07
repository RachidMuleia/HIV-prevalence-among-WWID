# HIV prevalence among WWID

R code to reproduce the tables in:

Banze AR, Muleia R, Nuvunga S, Craveiro I, Seabra SG, Baltazar CS. *Prevalence and Risk Factors of HIV, HBV, and HCV Among Women Who Inject Drugs in Five Mozambican Cities: A Cross-Sectional Biological and Behavioral Survey, 2023–2024*.

This pack is a cleaned version of the analysis used for the published tables. Individual-level survey data are **not** included.

## Manuscript tables and figures

| Paper item | What it reports | Script | Output |
|---|---|---|---|
| **Table 1** | Crude and pooled RDS-weighted prevalence of HIV, HBsAg, HCV, and co-infections | `analysis/02_table1_prevalence.R` | `output/table1_overall_prevalence.xlsx` |
| **Table 2** | Stratified RDS-weighted prevalence (HIV, HBsAg, HCV) | `analysis/03_table2_stratified_regression.R` | `output/table2_*_weighted.xlsx` |
| **Table 2** | Rao–Scott chi-square tests | same | `output/table2_rao_scott_*.xlsx` |
| **Table 2** | HIV crude prevalence ratios (modified Poisson, robust SE) | same | `output/table2_hiv_cPR.xlsx` |
| **Table 2** | HIV adjusted prevalence ratios (modified Poisson GEE) | same | `output/table2_hiv_aPR_gee.xlsx` |
| **HIV cascade / 95–95–95** | Awareness, ART, viral suppression | `analysis/04_cascade.R` | `output/figure_hiv_cascade_95_95_95.xlsx` |
| **Figure 1 (recruitment)** | Coupon/screening/enrollment flow | `figures/Figure_recruitment_flow.tex` | compile with `pdflatex` |
| **RDS diagnostics** | Seeds, waves, degree by sex (Methods) | `analysis/05_rds_diagnostics.R` | `output/rds_diagnostics.xlsx` |

Participant characteristics in the Results (median age, education, housing, etc.) are unweighted tabulations of the same covariates in `analysis/01_prepare_data.R` (`descriptive_kp` logic in the original `ANALISES_WWID_AURIA.R`).

## Estimators (as in the paper)

1. **RDS domain estimates among women.** City RDS objects are built from the mixed-gender PWID sample. Gile successive-sampling estimates are obtained with `subset` equal to women, which rescales the city PWID size \(N\) by the network-weighted female share and then recomputes SS weights in the female subsample (`RDS` package).
2. **Pooled prevalence.** Site-specific estimates are combined with a random-effects meta-analysis (REML) in `metafor`, following Lohr (2014) for multi-site RDS.
3. **Associations.** Rao–Scott chi-square tests use `survey::svychisq` with recruiter as the cluster and city as the stratum. HIV prevalence ratios use modified Poisson regression: crude models with sandwich robust variance, and the adjusted model with GEE (`geepack::geeglm`, exchangeable correlation, recruiter clusters), unweighted as described in the Methods.

## How to run

From the repository root, with the analytic file available locally:

```bash
export WWID_DATA_FILE="/path/to/DATA_PID_CLUSTER.csv"
# optional: export WWID_N_BOOT=15000
Rscript run_all.R
```

Or run scripts one table at a time:

```r
source("analysis/02_table1_prevalence.R")
source("analysis/03_table2_stratified_regression.R")
source("analysis/04_cascade.R")
```

Required R packages: `RDS`, `metafor`, `survey`, `sandwich`, `lmtest`, `geepack`, `dplyr`, `rio`.

## Data

The analytic file is the BBS-PWID II dataset restricted to the variables used here (HIV/HBV/HCV results, recruiter IDs, reported network size, and covariates). It is not posted because it contains individual-level information on a criminalized key population. Access is through the Instituto Nacional de Saúde, Mozambique.

City PWID size estimates used for Gile SS (Maputo, Beira, Tete, Quelimane, Nampula):

`990, 1029, 3219, 681, 1792`

## Repository layout

```
R/00_config.R              Paths, population sizes, bootstrap size
R/01_functions.R           RDS domain estimator, meta-analysis, Rao–Scott, PR extraction
analysis/01_prepare_data.R Load data, recodes, city RDS objects
analysis/02_table1_prevalence.R
analysis/03_table2_stratified_regression.R
analysis/04_cascade.R
analysis/05_rds_diagnostics.R
figures/Figure_recruitment_flow.tex
run_all.R
```

## Original working scripts

The exploratory local scripts (`ANALISES_WWID_AURIA.R`, `AURIA_ANALISE.R`, `MetaAnalysis.R`, `METANALYSIS_FUN.R`) contain the same estimators plus unused sensitivity runs (full-sample vs subset weights, homophily, toy RDS). This repository keeps only the code needed to reproduce the paper tables.
