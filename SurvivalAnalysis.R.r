# Cancer Drug Survival Analysis using Pharmaverse
# Install required packages (run once)
# install.packages(c("admiral", "dplyr", "survival", "survminer", "ggplot2"))

library(admiral)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# **************************************************************************
# CREATE TOY SDTM DATASETS (CDISC standards)
# **************************************************************************
# CDISC - Clinical data interchange standards consortium
# ADSL - Subject Level Analysis Dataset
# ADSL Subject - Level Analysis Dataset. It serves as source of STDM
# TRT01P = planned treatment
# TRT01A = treatment actually received
# SAFFL, ITTFL, EFFFL Safety, Intent to treat and efficacy pop flag

adsl <- tibble(
  STUDYID = "ONCO-001",
  USUBJID = paste0("SUBJ-", sprintf("%03d", 1:200)),
  SUBJID = sprintf("%03d", 1:200),
  ARM = sample(c("Placebo", "Drug 50mg", "Drug 100mg"), 200, replace = TRUE),
  ARMCD = case_when(
    ARM == "Placebo" ~ "PBO",
    ARM == "Drug 50mg" ~ "D50",
    ARM == "Drug 100mg" ~ "D100"
  ),
  TRT01P = ARM,
  TRT01A = ARM,
  AGE = round(rnorm(200, mean = 62, sd = 10)),
  SEX = sample(c("M", "F"), 200, replace = TRUE, prob = c(0.55, 0.45)),
  RACE = sample(c("WHITE", "BLACK OR AFRICAN AMERICAN", "ASIAN", "OTHER"), 
                200, replace = TRUE, prob = c(0.7, 0.15, 0.1, 0.05)),
  ETHNIC = sample(c("HISPANIC OR LATINO", "NOT HISPANIC OR LATINO"), 
                  200, replace = TRUE, prob = c(0.2, 0.8)),
  COUNTRY = sample(c("USA", "CANADA", "UK", "GERMANY"), 200, replace = TRUE),
  RANDDT = as.Date("2023-01-01") + sample(0:180, 200, replace = TRUE)
) %>%
  mutate(
    SAFFL = "Y",
    ITTFL = "Y",
    EFFFL = "Y"
  )

# Create survival outcomes based on treatment (better outcomes for higher doses)
survival_data <- adsl %>%
  mutate(
    # Simulate survival times with treatment effect
    baseline_hazard = case_when(
      ARMCD == "PBO" ~ 1.0,
      ARMCD == "D50" ~ 0.7,   # 30% hazard reduction
      ARMCD == "D100" ~ 0.5   # 50% hazard reduction
    ),
    # Add age effect (older = higher hazard)
    # λ(t∣age)=λ0(t)eβ×age  (here =λ0(t) is the  hazard  function,
    # beta is the regression coefficient representing the age effect )
    age_effect = exp((AGE - 62) * 0.02),
    total_hazard = baseline_hazard * age_effect,
    
    # Generate survival times from exponential distribution
    SURVTIME = rexp(n(), rate = total_hazard / 12),  # median ~12 months for placebo
   
     # Generate censoring times
    CENSORTIME = runif(n(), min = 6, max = 24),
    # Observed time is minimum of survival and censoring
    AVAL = pmin(SURVTIME, CENSORTIME),
    # Event indicator (1 = death, 0 = censored)
    CNSR = ifelse(SURVTIME <= CENSORTIME, 0, 1),
    EVENT = 1 - CNSR
  )

# ****************************************************************************
# CREATE ADTTE (Time-to-Event Analysis Dataset) using Admiral
# ****************************************************************************
adtte <- survival_data %>%
  mutate(
    PARAMCD = "OS",
    PARAM = "Overall Survival",
    AVALU = "MONTHS",
    STARTDT = RANDDT,
    ADT = RANDDT + round(AVAL * 30.4375),  # Convert months to days
    EVNTDESC = ifelse(EVENT == 1, "Death", "Censored")
  ) %>%
  select(STUDYID, USUBJID, SUBJID, PARAMCD, PARAM, 
         ARM, ARMCD, TRT01A, AGE, SEX,
         STARTDT, ADT, AVAL, AVALU, CNSR, EVENT, EVNTDESC)


# ****************************************************************************
# SURVIVAL ANALYSIS
# ****************************************************************************

# Create survival object
# Refere Biostatistics text book for interpretation of Kaplan-Meier Estimate
# Refer Biostatistics text book to know how confidence interval of Kaplan-Meier
# Estimate confidence interval is calculated

surv_obj <- Surv(time = adtte$AVAL, event = adtte$EVENT)

# Kaplan-Meier survival curves by treatment
km_fit <- survfit(surv_obj ~ ARM, data = adtte)

# Print summary statistics
print("Kaplan-Meier Summary by Treatment Arm:")
print(km_fit)

# Median survival times
print("\nMedian Survival Times (months):")
print(summary(km_fit)$table)

# ****************************************************************************
# COX PROPORTIONAL HAZARDS MODEL
# ****************************************************************************

# Uni variate Cox model - Treatment effect
# Cox proportion hazards regression model for survival data
# It models the hazard (instantaneous event rate) at time t as a product of an 
# unspecified baseline hazard and a parametric function of covariates as follows
# h(t∣X)=h0(t ) ×exp(beta1X1+beta2X2+..+betapXp)
cox_treatment <- coxph(surv_obj ~ ARM, data = adtte)
print("\nCox Model - Treatment Effect:")
print(summary(cox_treatment))

# Multivariable Cox model - Adjusted for age and sex
cox_adjusted <- coxph(surv_obj ~ ARM + AGE + SEX, data = adtte)
print("\nAdjusted Cox Model (Treatment + Age + Sex):")
print(summary(cox_adjusted))

# ****************************************************************************
# LOG-RANK TEST
# ****************************************************************************

logrank_test <- survdiff(surv_obj ~ ARM, data = adtte)
print("\nLog-Rank Test for Treatment Comparison:")
print(logrank_test)

# ***************************************************************************
# VISUALIZATIONS
# ***************************************************************************

# Kaplan-Meier Curves
km_plot <- ggsurvplot(
  km_fit,
  data = adtte,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  ggtheme = theme_bw(),
  palette = c("#E7B800", "#2E9FDF", "#00AFBB"),
  title = "Overall Survival by Treatment Arm",
  xlab = "Time (Months)",
  ylab = "Overall Survival Probability",
  legend.title = "Treatment",
  legend.labs = levels(factor(adtte$ARM))
)

print(km_plot)

# Forest plot for hazard ratios - using positional indexing
cox_coef <- coef(cox_treatment)
cox_ci <- confint(cox_treatment)

hr_data <- tibble(
  Comparison = c("Drug 50mg vs Placebo", "Drug 100mg vs Placebo"),
  HR = exp(cox_coef),  # This automatically gets both coefficients
  Lower = exp(cox_ci[, 1]),  # Lower bounds
  Upper = exp(cox_ci[, 2])   # Upper bounds
)

forest_plot <- ggplot(hr_data, aes(x = HR, y = Comparison)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(trans = "log", breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2)) +
  labs(
    title = "Hazard Ratios for Overall Survival",
    x = "Hazard Ratio (95% CI)",
    y = ""
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

print(forest_plot)


# ****************************************************************************
# SURVIVAL RATES AT KEY TIMEPOINTS
# ****************************************************************************

# Calculate survival rates at 6, 12, and 18 months
timepoints <- c(6, 12, 18)
surv_summary <- summary(km_fit, times = timepoints)

survival_table <- tibble(
  Time = rep(timepoints, each = length(unique(adtte$ARM))),
  Treatment = rep(levels(factor(adtte$ARM)), times = length(timepoints)),
  Survival = surv_summary$surv,
  Lower_CI = surv_summary$lower,
  Upper_CI = surv_summary$upper
) %>%
  mutate(
    Survival_Pct = sprintf("%.1f%% (%.1f%%-%.1f%%)", 
                           Survival * 100, Lower_CI * 100, Upper_CI * 100)
  )

print("\nSurvival Rates at Key Timepoints:")
print(survival_table)

# ****************************************************************************
# EXPORT SUMMARY STATISTICS
# ****************************************************************************

# Create summary table for reporting
summary_stats <- adtte %>%
  group_by(ARM) %>%
  summarise(
    N = n(),
    Events = sum(EVENT),
    Censored = sum(CNSR),
    `Event Rate (%)` = round(Events / N * 100, 1),
    `Median FU (months)` = round(median(AVAL), 1),
    `Max FU (months)` = round(max(AVAL), 1)
  )

print("\nSummary Statistics by Treatment Arm:")
print(summary_stats)

cat("\n=== Analysis Complete ===\n")
cat("Note: This is simulated data for demonstration purposes.\n")
cat("In real analysis, ensure proper data validation and quality checks.\n")