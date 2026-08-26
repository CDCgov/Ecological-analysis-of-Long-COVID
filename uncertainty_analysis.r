# Author: X.Z.
# Created: 2025-09
# Last Modified: 2026-06-15
# Xinmeng Zhao, Li Deng, Nicole D. Ford, and Sharon Saydah. Long COVID Prevalence Among U.S. Adults: A State-Level Ecological Analysis of Associations with SARS-CoV-2 Incidence, COVID-19 Hospitalization Rates, COVID-19 Vaccine Coverage, and Multimorbidity. (Preprint, PLOS ONE).
library(dplyr)
library(mitools)

analysis_data <- read.csv("Data/analysis_data.csv")

# Set uncertainty_source to "Y", "X", or "both".
# Optional command-line override: Rscript uncertainty_analysis.r X
uncertainty_source <- "both"
command_args <- commandArgs(trailingOnly = TRUE)
if (length(command_args) >= 1 && nzchar(command_args[[1]])) {
  uncertainty_source <- command_args[[1]]
}

# ==============================================================================
# Helper Functions
# ==============================================================================

normalize_uncertainty_source <- function(source) {
  source <- tolower(trimws(source))
  valid_sources <- c("x", "y", "both")

  if (!source %in% valid_sources) {
    stop(
      paste(
        "uncertainty_source must be one of:",
        paste(valid_sources, collapse = ", ")
      )
    )
  }

  source
}

#' Uncertainty analysis via outcome/predictor resampling and Rubin's Rules pooling.
#'
#' @param mean_vec   Numeric vector of outcome point estimates (one per state).
#' @param se_vec     Numeric vector of outcome standard errors.
#' @param factors Character vector of factor column names.
#' @param factor_data Data frame containing factors and state identifier.
#' @param predictor_se_cols Optional named character vector mapping predictor
#'                          names to their standard-error columns in
#'                          \code{factor_data}. Each listed predictor is
#'                          resampled in each simulated dataset.
#' @param uncertainty_source Character value indicating uncertainty sources:
#'                           "Y" for outcome only, "X" for predictor only, or
#'                           "both" for outcome and predictor uncertainty.
#' @param weights    Optional numeric vector of regression weights (same length
#'                   and order as \code{factor_data}). Pass \code{NULL} for
#'                   unweighted models.
#' @param state_id   Name of the state identifier column in \code{factor_data}.
#' @param n_imp      Number of simulated outcome draws (default 1000).
#' @param seed       Random seed for reproducibility (default 123).
#'
#' @return Data frame with pooled coefficients, standard errors, CIs, and p-values.

run_uncertainty_analysis <- function(mean_vec, se_vec, factors,
                                     factor_data, predictor_se_cols = NULL,
                                     uncertainty_source = "both",
                                     weights = NULL,
                                     state_id, n_imp = 1000, seed = 123) {
  set.seed(seed)
  uncertainty_source <- normalize_uncertainty_source(uncertainty_source)

  if (uncertainty_source %in% c("y", "both")) {
    resampled_Y <- replicate(n_imp, {
      rnorm(n = length(mean_vec), mean = mean_vec, sd = se_vec)
    })
  } else {
    resampled_Y <- matrix(
      rep(mean_vec, times = n_imp),
      nrow = length(mean_vec),
      ncol = n_imp
    )
  }
  rownames(resampled_Y) <- factor_data[[state_id]]

  X <- factor_data[, factors]
  rownames(X) <- factor_data[[state_id]]

  resampled_predictors <- list()
  if (uncertainty_source %in% c("x", "both")) {
    if (is.null(predictor_se_cols) || length(predictor_se_cols) == 0) {
      stop("uncertainty_source includes X, but no predictor_se_cols were provided.")
    }

    for (predictor in names(predictor_se_cols)) {
      if (!predictor %in% factors) {
        stop(paste("Predictor", predictor, "is not in the factor list."))
      }

      se_col <- predictor_se_cols[[predictor]]
      if (!se_col %in% names(factor_data)) {
        stop(paste("Standard-error column", se_col, "not found in factor_data."))
      }

      resampled_predictors[[predictor]] <- replicate(n_imp, {
        rnorm(
          n    = nrow(factor_data),
          mean = factor_data[[predictor]],
          sd   = factor_data[[se_col]]
        )
      })
      rownames(resampled_predictors[[predictor]]) <- factor_data[[state_id]]
    }
  }

  if (!is.null(weights)) {
    weight_vec <- weights
    names(weight_vec) <- factor_data[[state_id]]
    weight_vec <- weight_vec[rownames(X)]
  }

  regression_results <- lapply(seq_len(n_imp), function(i) {
    y_sim    <- resampled_Y[, i]
    X_sim <- X

    if (length(resampled_predictors) > 0) {
      for (predictor in names(resampled_predictors)) {
        X_sim[[predictor]] <- resampled_predictors[[predictor]][rownames(X_sim), i]
      }
    }

    df_model <- cbind(X_sim, y_sim = y_sim)
    df_model <- df_model[complete.cases(df_model), ]

    if (!is.null(weights)) {
      weight_use <- weight_vec[rownames(df_model)]
      lm(y_sim ~ ., data = df_model, weights = weight_use)
    } else {
      lm(y_sim ~ ., data = df_model)
    }
  })

  betas  <- lapply(regression_results, coef)
  vars   <- lapply(regression_results, vcov)
  pooled <- MIcombine(betas, vars)

  estimates <- pooled$coefficients
  se        <- sqrt(diag(pooled$variance))
  t_vals    <- estimates / se
  df_vals   <- pooled$df
  p_vals    <- 2 * pt(-abs(t_vals), df = df_vals)
  ci_lower  <- estimates - qt(0.975, df_vals) * se
  ci_upper  <- estimates + qt(0.975, df_vals) * se

  data.frame(
    factor = names(estimates),
    Beta      = estimates,
    StdError  = se,
    t         = t_vals,
    df        = df_vals,
    CI_Lower  = ci_lower,
    CI_Upper  = ci_upper,
    p         = p_vals,
    check.names = FALSE
  )
}

# ==============================================================================
# Factor Lists
# ==============================================================================
factor_list <- c(
  "cumulative_incidence",
  "cumulative_hospitalization",
  "cumulative_booster",
  "ThreePlusComorbidities",
  "age_65_over_share",
  "male_given_age_65_over",
  "female_given_age_65_under"
)

# ==============================================================================
# Interpretable Scale Configuration
# ==============================================================================
factor_reporting_increments <- c(
  cumulative_incidence        = 0.01,
  cumulative_hospitalization  = 0.001,
  cumulative_booster      = 0.10,
  ThreePlusComorbidities      = 0.10,
  age_65_over_share           = 0.10,
  male_given_age_65_over      = 0.10,
  female_given_age_65_under   = 0.10
)

factor_reporting_labels <- c(
  cumulative_incidence        = "per 10 reported cases/ 1,000 persons-year",
  cumulative_hospitalization  = "per 1 hospitalizations/ 1,000 persons-year",
  cumulative_booster      = "per 10 percentage-point increase",
  ThreePlusComorbidities      = "per 10 percentage-point increase",
  age_65_over_share           = "per 10 percentage-point increase",
  male_given_age_65_over      = "per 10 percentage-point increase",
  female_given_age_65_under   = "per 10 percentage-point increase"
)

get_outcome_multiplier <- function(model_name) {
  rep(100, length(model_name))
}

add_interpretable_scale <- function(results) {
  results %>%
    mutate(
      factor_Increment = unname(factor_reporting_increments[factor]),
      Reported_Increment  = unname(factor_reporting_labels[factor]),
      Outcome_Unit        = "percentage-point change in Long COVID prevalence",
      Outcome_Multiplier  = get_outcome_multiplier(Model),
      Scaled_Beta         = Beta     * factor_Increment * Outcome_Multiplier,
      Scaled_CI_Lower     = CI_Lower * factor_Increment * Outcome_Multiplier,
      Scaled_CI_Upper     = CI_Upper * factor_Increment * Outcome_Multiplier
    )
}

# ==============================================================================
# Run Uncertainty Analyses
# ==============================================================================

uncertainty_source <- normalize_uncertainty_source(uncertainty_source)
message("Uncertainty source: ", uncertainty_source)

# Scenario: Crude
summary_1 <- run_uncertainty_analysis(
  mean_vec       = analysis_data$Crude_Long_COVID,
  se_vec         = analysis_data$Crude_Long_COVID_StdErr,
  factors     = factor_list,
  factor_data = analysis_data,
  predictor_se_cols = c(cumulative_booster = "cumulative_booster_SE"),
  uncertainty_source = uncertainty_source,
  weights        = NULL,
  state_id       = "State"
)

# ==============================================================================
# Prepare and Export
# ==============================================================================
all_results <- summary_1 %>%
  mutate(Model = "Crude") %>%
  filter(factor != "(Intercept)") %>%
  select(Model, factor, Beta, CI_Lower, CI_Upper, p)

write.csv(all_results, "results/uncertainty_analysis.csv", row.names = FALSE)

scaled_results <- add_interpretable_scale(all_results)
write.csv(scaled_results, "results/uncertainty_analysis_interpretable_scale.csv", row.names = FALSE)
