# Author: X.Z.
# Created: 2025-09
# Last Modified: 2026-06-15
# Xinmeng Zhao, Li Deng, Nicole D. Ford, and Sharon Saydah. Long COVID Prevalence Among U.S. Adults: A State-Level Ecological Analysis of Associations with SARS-CoV-2 Incidence, COVID-19 Hospitalization Rates, COVID-19 Vaccine Coverage, and Multimorbidity. (Preprint, PLOS ONE).
library(dplyr)
library(ggplot2)
library(ggrepel)
library(GGally)
library(reshape2)

analysis_data <- read.csv("Data/analysis_data.csv")

# ==============================================================================
# Variable Configuration
# ==============================================================================
factor_vars <- c(
  "cumulative_incidence",
  "cumulative_hospitalization",
  "cumulative_booster",
  "ThreePlusComorbidities",
  "age_65_over_share",
  "male_given_age_65_over",
  "female_given_age_65_under"
)

factor_labels <- c(
  cumulative_incidence  = "SARS-CoV-2 incidence per 1,000 persons-year",
  cumulative_hospitalization  = "COVID-19 hospitalizations per 1,000 persons-year",
  cumulative_booster          = "COVID-19 vaccine coverage",
  ThreePlusComorbidities          = "Prevalence of Multimorbidity",
  age_65_over_share               = "Proportion of Adult >=65 yrs",
  male_given_age_65_over          = "Proportion of male >=65 yrs",
  female_given_age_65_under       = "Proportion of female <65 yrs"
)

correlation_labels <- c(
  cumulative_incidence  = "SARS-CoV-2\nincidence",
  cumulative_hospitalization  = "COVID-19\nhospitalizations",
  cumulative_booster          = "COVID-19 Vaccine\ncoverage",
  ThreePlusComorbidities          = "Prevalence of Multimorbidity",
  age_65_over_share               = "Proportion of Adults\n>=65 yrs",
  male_given_age_65_over          = "Proportion of male\n>=65 yrs",
  female_given_age_65_under       = "Proportion of female\n<65 yrs"
)

scatter_file_names <- c(
  cumulative_incidence  = "scatterplot_SARSCoV2_incidence_longCOVID.tiff",
  cumulative_hospitalization  = "scatterplot_COVID19_hospitalization_longCOVID.tiff",
  cumulative_booster          = "scatterplot_CV_longCOVID.tiff",
  ThreePlusComorbidities          = "scatterplot_multimorbidity_longCOVID.tiff",
  age_65_over_share               = "scatterplot_age65_longCOVID.tiff",
  male_given_age_65_over          = "scatterplot_male_age65over_longCOVID.tiff",
  female_given_age_65_under       = "scatterplot_female_age65under_longCOVID.tiff"
)

percent_vars <- c(
  "cumulative_booster",
  "ThreePlusComorbidities",
  "age_65_over_share",
  "male_given_age_65_over",
  "female_given_age_65_under"
)

rate_per_1000_vars <- c(
  "cumulative_incidence",
  "cumulative_hospitalization"
)

# ==============================================================================
# Helper Functions
# ==============================================================================
save_high_res_tiff <- function(plot, filename, width = 7, height = 5) {
  ggsave(
    filename    = filename,
    plot        = plot,
    width       = width,
    height      = height,
    units       = "in",
    dpi         = 600,
    device      = "tiff",
    compression = "lzw"
  )
}

# ==============================================================================
# Scatter Plots
# ==============================================================================
make_scatter_plot <- function(var) {
  plot_data <- analysis_data %>%
    mutate(plot_x = if (var %in% rate_per_1000_vars) .data[[var]] * 1000 else .data[[var]])

  p <- ggplot(plot_data, aes(x = plot_x, y = Crude_Long_COVID)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, level = 0.95,
                color = "#1f4e79", fill = "#9bbce0") +
    geom_text_repel(aes(label = State), max.overlaps = 20, size = 3) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = factor_labels[[var]],
      y = "Long COVID prevalence"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin      = margin(8, 12, 8, 8)
    )

  if (var %in% percent_vars) {
    p <- p + scale_x_continuous(labels = scales::percent_format(accuracy = 1))
  } else if (var %in% rate_per_1000_vars) {
    p <- p + scale_x_continuous(labels = scales::number_format(accuracy = 0.1))
  } else {
    p <- p + scale_x_continuous(labels = scales::number_format(accuracy = 0.001))
  }

  p
}

scatter_dir <- "results/scatterplot"
dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)

for (var in factor_vars) {
  p <- make_scatter_plot(var)
  print(p)
  save_high_res_tiff(p, file.path(scatter_dir, scatter_file_names[[var]]))
}

# ==============================================================================
# Pair Plots
# ==============================================================================
GGally::ggpairs(analysis_data[, c(factor_vars, "Crude_Long_COVID")])

# ==============================================================================
# Long COVID State-Level Prevalence Plot
# ==============================================================================
long_covid_state_plot <- ggplot(
  analysis_data,
  aes(x = Crude_Long_COVID * 100, y = reorder(State, Crude_Long_COVID))
) +
  geom_errorbarh(
    aes(
      xmin = (Crude_Long_COVID - Crude_Long_COVID_StdErr) * 100,
      xmax = (Crude_Long_COVID + Crude_Long_COVID_StdErr) * 100
    ),
    height    = 0.2,
    linewidth = 0.4,
    color     = "#5A5A5A"
  ) +
  geom_point(size = 2.4, color = "#1F4E79") +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = "Crude Long COVID Prevalence (%)",
    y = "State"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y      = element_text(size = 9.5),
    axis.text.x      = element_text(size = 10.5),
    axis.title       = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(long_covid_state_plot)

results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
save_high_res_tiff(
  long_covid_state_plot,
  file.path(results_dir, "long_covid_state_level_prevalence.tiff"),
  width  = 7,
  height = 10
)

# ==============================================================================
# factor–Outcome Correlation Bar Chart
# ==============================================================================
correlations <- data.frame(
  Variable    = unname(factor_labels[factor_vars]),
  Correlation = sapply(
    factor_vars,
    function(var) cor(analysis_data[[var]], analysis_data$Crude_Long_COVID,
                      use = "complete.obs")
  )
)

ggplot(correlations,
       aes(x = reorder(Variable, Correlation), y = Correlation,
           fill = Correlation > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("Negative", "Positive")) +
  labs(x = "Potential Risk Factors Associated with Long COVID",
       y = "Correlation Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

# ==============================================================================
# Correlation Matrix (Spearman)
# ==============================================================================
subset_data <- analysis_data %>% select(all_of(factor_vars))
colnames(subset_data) <- unname(correlation_labels[factor_vars])

cor_matrix <- cor(subset_data, use = "complete.obs", method = "spearman")
cor_melted  <- reshape2::melt(cor_matrix)

correlation_matrix_plot <- ggplot(cor_melted,
                                   aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "#F2F2F2", linewidth = 0.25) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(
    low      = "#B66A6A",
    mid      = "#F7F7F4",
    high     = "#1F4E79",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Spearman\ncorrelation"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  labs(title = "", x = "", y = "") +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1,
                                  size = 9, color = "#3F3F3F"),
    axis.text.y   = element_text(size = 9, color = "#3F3F3F"),
    panel.grid    = element_blank(),
    legend.position = "right",
    legend.title  = element_text(size = 9),
    legend.text   = element_text(size = 8),
    plot.margin   = margin(6, 6, 6, 6)
  )

print(correlation_matrix_plot)

correlation_dir <- "results/correlation"
dir.create(correlation_dir, recursive = TRUE, showWarnings = FALSE)

save_high_res_tiff(
  correlation_matrix_plot,
  file.path(correlation_dir, "correlation_between_factors.tiff"),
  width  = 6.5,
  height = 5.5
)
ggsave(
  filename = file.path(correlation_dir, "correlation_between_factors.png"),
  plot     = correlation_matrix_plot,
  width    = 6.5,
  height   = 5.5,
  units    = "in",
  dpi      = 600
)
