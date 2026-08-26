# Cite
Xinmeng Zhao, Li Deng, Nicole D. Ford, and Sharon Saydah. Long COVID Prevalence Among U.S. Adults: A State-Level Ecological Analysis of Associations with SARS-CoV-2 Incidence, COVID-19 Hospitalization Rates, COVID-19 Vaccine Coverage, and Multimorbidity. (Preprint, PLOS ONE).

# Code Overview

This folder contains the analysis scripts for the ecological analysis of state-level Long COVID prevalence. The scripts use `Data/analysis_data.csv` as the input file and write figures or tables to a `results/` folder.


## Files

### `state_level_analysis.r`

Creates descriptive state-level figures and correlation summaries:

- Scatter plots of each factor against crude Long COVID prevalence.
- A pair plot of factors and crude Long COVID prevalence.
- A state-level crude Long COVID prevalence plot with standard error bars.
- A factor-outcome correlation bar chart.
- A Spearman correlation heatmap among factors.

Main outputs are saved under:

- `results/scatterplot/`
- `results/correlation/`
- `results/long_covid_state_level_prevalence.tiff`

### `ecological_regression.r`

Runs ecological regression analyses using crude Long COVID prevalence as the outcome:

- Descriptive summaries for factor variables.
- Univariable linear regressions for each factor.
- A multivariable crude Long COVID regression model.
- Interpretable-scale regression summaries.
- Variance inflation factor diagnostics.
- Leave-one-state-out sensitivity analyses.

Main outputs are saved under:

- `results/regression_diagnostics/factor_descriptive_summary.csv`
- `results/regression_diagnostics/long_covid_regression_summary.csv`
- `results/regression_diagnostics/long_covid_regression_summary_interpretable_scale.csv`
- `results/regression_diagnostics/vif_results.csv`
- `results/regression_diagnostics/leave_one_state_out_results.csv`
- `results/regression_diagnostics/leave_one_state_out_summary.csv`

### `uncertainty_analysis.r`

Runs an uncertainty analysis for crude Long COVID prevalence using both outcome (Y) and predictor (X) uncertainty by default. For each of 1,000 simulated datasets, the script resamples `Crude_Long_COVID` from its state-level mean and `Crude_Long_COVID_StdErr`, resamples `cumulative_booster` from its observed value and `cumulative_booster_SE`, keeps other factors without standard-error columns fixed, fits the full multivariable model, and pools coefficients with Rubin's Rules. The uncertainty source can also be limited to `X` or `Y` with a command-line argument.

Main outputs are saved under:

- `results/uncertainty_analysis.csv`
- `results/uncertainty_analysis_interpretable_scale.csv`

## Running The Scripts

Run each script from the folder so the relative path `Data/analysis_data.csv` resolves correctly.

Suggested order:

1. `state_level_analysis.r`
2. `ecological_regression.r`
3. `uncertainty_analysis.r`
