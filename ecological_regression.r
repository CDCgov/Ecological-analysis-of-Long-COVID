# Author: X.Z.
# Created: 2025-09
# Last Modified: 2026-06-15
# Cite: Xinmeng Zhao, Li Deng, Nicole D. Ford, and Sharon Saydah. Long COVID Prevalence Among U.S. Adults: A State-Level Ecological Analysis of Associations with SARS-CoV-2 Incidence, COVID-19 Hospitalization Rates, COVID-19 Vaccine Coverage, and Multimorbidity. (Preprint, PLOS ONE).

library(broom)
library(car)
library(dplyr)
library(lmtest)
library(purrr)
library(sandwich)
library(tidyr)

analysis_data <- read.csv("Data/analysis_data.csv")
rownames(analysis_data) <- analysis_data$State

diagnostics_dir <- "results/regression_diagnostics"
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Factor Configuration
# ==============================================================================
factors <- c(
  "cumulative_incidence",
  "cumulative_hospitalization",
  "cumulative_booster",
  "ThreePlusComorbidities",
  "age_65_over_share",
  "male_given_age_65_over",
  "female_given_age_65_under"
)

factor_reporting_increments <- c(
  cumulative_incidence  = 0.01,
  cumulative_hospitalization  = 0.001,
  cumulative_booster          = 0.10,
  ThreePlusComorbidities          = 0.10,
  age_65_over_share               = 0.10,
  male_given_age_65_over          = 0.10,
  female_given_age_65_under       = 0.10
)

factor_reporting_labels <- c(
  cumulative_incidence  = "per 10 reported cases/ 1,000 persons-year",
  cumulative_hospitalization  = "per 1 hospitalizations/ 1,000 persons-year",
  cumulative_booster          = "per 10 percentage-point increase",
  ThreePlusComorbidities          = "per 10 percentage-point increase",
  age_65_over_share               = "per 10 percentage-point increase",
  male_given_age_65_over          = "per 10 percentage-point increase",
  female_given_age_65_under       = "per 10 percentage-point increase"
)

# ==============================================================================
# Model Specifications
# ==============================================================================
model_specs <- data.frame(
  Model = "Crude",
  Outcome = "Crude_Long_COVID",
  Weighted = FALSE,
  stringsAsFactors = FALSE
)

# ==============================================================================
# Helper Functions
# ==============================================================================
make_model_formula <- function(outcome_var, model_factors) {
  as.formula(paste(outcome_var, "~", paste(model_factors, collapse = " + ")))
}

fit_model <- function(data, outcome_var, weighted = TRUE,
                      model_factors = factors) {
  model_formula <- make_model_formula(outcome_var, model_factors)
  if (weighted) {
    lm(model_formula, data = data, weights = total_population)
  } else {
    lm(model_formula, data = data)
  }
}

get_outcome_multiplier <- function(model_name) {
  # Crude_Long_COVID is stored as a proportion.
  rep(100, length(model_name))
}

add_interpretable_scale <- function(results) {
  results %>%
    mutate(
      factor_Increment = unname(factor_reporting_increments[factor]),
      Reported_Increment  = unname(factor_reporting_labels[factor]),
      Outcome_Unit        = "percentage-point change in Long COVID prevalence",
      Outcome_Multiplier  = get_outcome_multiplier(Model),
      Scaled_Beta         = Beta      * factor_Increment * Outcome_Multiplier,
      Scaled_CI_Lower     = CI_Lower  * factor_Increment * Outcome_Multiplier,
      Scaled_CI_Upper     = CI_Upper  * factor_Increment * Outcome_Multiplier
    )
}

# ==============================================================================
# Descriptive Summary
# ==============================================================================
summary_table <- analysis_data %>%
  select(all_of(factors)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    `Mean (SD)`      = paste0(round(mean(value, na.rm = TRUE), 3),
                               " (", round(sd(value, na.rm = TRUE), 3), ")"),
    `Median [Q1-Q3]` = paste0(round(median(value, na.rm = TRUE), 3),
                               " [", round(quantile(value, 0.25, na.rm = TRUE), 3),
                               " - ", round(quantile(value, 0.75, na.rm = TRUE), 3), "]"),
    `Range`          = paste0(round(min(value, na.rm = TRUE), 3), " - ",
                               round(max(value, na.rm = TRUE), 3)),
    .groups = "drop"
  )

write.csv(summary_table,
          file.path(diagnostics_dir, "factor_descriptive_summary.csv"),
          row.names = FALSE)

# ==============================================================================
# Univariable Models
# ==============================================================================
fit_univariate_model <- function(var) {
  lm(as.formula(paste("Crude_Long_COVID ~", var)), data = analysis_data) %>%
    tidy(conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    transmute(
      Model    = "Univariable",
      factor = term,
      Beta      = estimate,
      CI_Lower  = conf.low,
      CI_Upper  = conf.high,
      P_value   = p.value
    )
}

univariate_summary <- map_dfr(factors, fit_univariate_model)

# ==============================================================================
# Multivariable Models
# ==============================================================================
fitted_models <- pmap(
  model_specs,
  function(Model, Outcome, Weighted) fit_model(analysis_data, Outcome, Weighted)
)
names(fitted_models) <- model_specs$Model

multivariable_summary <- imap_dfr(
  fitted_models,
  function(model, model_name) {
    tidy(model, conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      transmute(
        Model     = model_name,
        factor = term,
        Beta      = estimate,
        CI_Lower  = conf.low,
        CI_Upper  = conf.high,
        P_value   = p.value
      )
  }
)


all_results <- bind_rows(univariate_summary, multivariable_summary)
write.csv(all_results,
          file.path(diagnostics_dir, "long_covid_regression_summary.csv"),
          row.names = FALSE)

scaled_regression_summary <- add_interpretable_scale(all_results)
write.csv(scaled_regression_summary,
          file.path(diagnostics_dir, "long_covid_regression_summary_interpretable_scale.csv"),
          row.names = FALSE)


# ==============================================================================
# Variance Inflation Factors
# ==============================================================================
calculate_vif <- function(model, model_name) {
  vif_values <- car::vif(model)
  if (is.matrix(vif_values)) {
    return(data.frame(
      Model     = model_name,
      factor = rownames(vif_values),
      vif_values,
      row.names    = NULL,
      check.names  = FALSE
    ))
  }
  data.frame(
    Model            = model_name,
    factor        = names(vif_values),
    VIF              = unname(vif_values),
    stringsAsFactors = FALSE
  )
}

vif_results <- imap_dfr(fitted_models, calculate_vif)
write.csv(vif_results,
          file.path(diagnostics_dir, "vif_results.csv"),
          row.names = FALSE)

# ==============================================================================
# Leave-One-State-Out Sensitivity Analysis
# ==============================================================================
calculate_leave_one_state_out <- function(model_name, outcome_var, weighted) {
  full_model <- fit_model(analysis_data, outcome_var, weighted)
  full_beta  <- coef(full_model)

  map_dfr(
    analysis_data$State,
    function(omitted_state) {
      loo_data <- analysis_data %>% filter(State != omitted_state)
      rownames(loo_data) <- loo_data$State

      loo_model <- fit_model(loo_data, outcome_var, weighted)

      tidy(loo_model, conf.int = TRUE) %>%
        filter(term != "(Intercept)") %>%
        transmute(
          Model                   = model_name,
          Omitted_State           = omitted_state,
          factor               = term,
          Beta                    = estimate,
          CI_Lower                = conf.low,
          CI_Upper                = conf.high,
          P_value                 = p.value,
          Full_Data_Beta          = unname(full_beta[term]),
          Beta_Difference         = Beta - Full_Data_Beta,
          Percent_Change_From_Full = ifelse(
            abs(Full_Data_Beta) > .Machine$double.eps,
            100 * Beta_Difference / abs(Full_Data_Beta),
            NA_real_
          )
        )
    }
  )
}

leave_one_state_out_results <- pmap_dfr(
  model_specs,
  function(Model, Outcome, Weighted) {
    calculate_leave_one_state_out(Model, Outcome, Weighted)
  }
)

write.csv(leave_one_state_out_results,
          file.path(diagnostics_dir, "leave_one_state_out_results.csv"),
          row.names = FALSE)

leave_one_state_out_summary <- leave_one_state_out_results %>%
  group_by(Model, factor) %>%
  summarise(
    Full_Data_Beta               = first(Full_Data_Beta),
    Min_Beta                     = min(Beta, na.rm = TRUE),
    Max_Beta                     = max(Beta, na.rm = TRUE),
    Max_Absolute_Percent_Change  = max(abs(Percent_Change_From_Full), na.rm = TRUE),
    State_With_Max_Change        = Omitted_State[which.max(abs(Percent_Change_From_Full))],
    Same_Sign_All_Leave_One_Out  = all(sign(Beta) == sign(Full_Data_Beta), na.rm = TRUE),
    Significant_All_Leave_One_Out = all(P_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(leave_one_state_out_summary,
          file.path(diagnostics_dir, "leave_one_state_out_summary.csv"),
          row.names = FALSE)
