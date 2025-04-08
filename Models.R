###############################################################
# Peace Agreement Analysis - Cox Regression Models
# Analyzing how economic provisions affect peace agreement duration
## - Code by Elisa D'Amico
## - Uploaded on 8 April 2025
###############################################################

# VARIABLES USED IN MODELS:
# Dependent Variable:
# - duration_in_days: Duration of peace agreement in days

# Key Independent Variables (Private Goods):
# - funds_dummy_vSum22: Direct funds for ex-combatants
# - inclusive_dummy_vSum22: Explicit reintegration for ex-combatants
# - tp_all_vSum22: Training programs for ex-combatants
# - PC1: Principal component derived from private goods variables

# Key Independent Variables (Public Goods):
# - EpsFis: Fiscal federalism provision (binary)
# - Dev: Development/socio-economic reconstruction (categorical)
# - NEC: National economic plan (binary)
# - PC1: Principal component derived from public goods variables

# Control Variables:
# - log_Lgt: Log of battle deaths (conflict intensity)
# - cumulative_intensity: Cumulative conflict intensity
# - log_GDP_PPP: Log of GDP per capita (PPP)
# - log_develassistance: Log of development assistance
# - v2x_libdem: Liberal democracy index
# - PrevAgmt_Bin: Previous agreement binary indicator
# - post_1990: Period indicator (after 1990)
# - post_2005: Period indicator (after 2005)
# - ddr: Disarmament, demobilization, and reintegration provisions
# - contype_2: Conflict type 2
# - contype_3: Conflict type 3
# - year: Year of agreement
# - conflict_id: Conflict identifier
# - ifeconint_treat: Treatment variable for clustering

###############################################################
# SETUP - Load Libraries and Data
###############################################################

# Clear environment
rm(list = ls())

# Load required libraries
library(broom)      # For converting statistical objects to tidy data frames
library(coxme)      # For mixed-effects Cox models
library(lubridate)  # For date handling
library(survival)   # For survival analysis
library(dplyr)      # For data manipulation
library(tidyr)      # For data tidying
library(sandwich)   # For robust standard errors
library(lmtest)     # For hypothesis testing
library(stargazer)  # For creating regression tables
library(MatchIt)    # For propensity score matching
library(zoo)        # For time series analysis
library(mice)       # For multiple imputation
library(purrr)      # For functional programming
library(nortest)    # For normality tests
library(MASS)       # For statistical functions
library(forestplot) # For forest plots
library(coefplot)   # For coefficient plots
library(arm)        # For Bayesian data analysis
library(caret)      # For classification and regression training
library(psych)      # For psychological research analyses
library(readr)      # For reading data
library(countrycode)# For country codes
library(sf)         # For spatial data handling
library(rnaturalearth)     # For map data
library(rnaturalearthdata) # For earth data
library(readxl)

# Load data
privgoods <- read_csv("privgoods.csv")
pubgoods <- read_xlsx("pubgoods.xlsx")

###############################################################
# PART 1: PRIVATE GOODS MODELS
###############################################################

## Principal Component Analysis for Private Goods
# Select variables for PCA
selected_columns <- c("funds_dummy_vSum22", "inclusive_dummy_vSum22", "tp_all_vSum22")
covariates_for_pca <- pubgoods[, selected_columns]

# Perform PCA
pca_result <- prcomp(covariates_for_pca, scale. = TRUE)

# Plot variance explained
plot(cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2)), 
     xlab = "Number of Principal Components",
     ylab = "Cumulative Proportion of Variance Explained", 
     type = "b")

# Add reference line at 90% cumulative variance
abline(h = 0.9, col = "red", lty = 2)

# Display factor loadings
factor_loadings <- pca_result$rotation
print(factor_loadings)

# Extract principal components and combine with original data
pca_components <- predict(pca_result, newdata = covariates_for_pca)
privgoods <- cbind(pubgoods, pca_components)

## MODEL 1: Cox model with Principal Component
cox_model_with_pca <- coxph(
  Surv(duration_in_days, censored) ~ 
    PC1 + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id + 
    cluster(ifeconint_treat),
  data = privgoods
)
summary(cox_model_with_pca)

## MODEL 2: Explicit Reintegration for Ex-Combatants
cox_mod_incl <- coxph(
  Surv(duration_in_days, censored) ~
    inclusive_dummy_vSum22 + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id + 
    cluster(ifeconint_treat),
  data = pubgoods
)
summary(cox_mod_incl)

## MODEL 3: Direct Funds for Ex-Combatants
cox_mod_funds <- coxph(
  Surv(duration_in_days, censored) ~
    funds_dummy_vSum22 + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id + 
    cluster(ifeconint_treat),
  data = pubgoods
)
summary(cox_mod_funds)

## MODEL 4: Economic Provisions
cox_mod_econ <- coxph(
  Surv(duration_in_days, censored) ~
    economic_dummy_vSum22 + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id + 
    cluster(ifeconint_treat),
  data = pubgoods
)
summary(cox_mod_econ)

## MODEL 5: Training Programs for Ex-Combatants
cox_mod_full <- coxph(
  Surv(duration_in_days, censored) ~
    tp_all_vSum22 + log_Lgt + cumulative_intensity + log_GDP_PPP + 
    log_develassistance + v2x_libdem + PrevAgmt_Bin + post_1990 + 
    post_2005 + ddr + contype_2 + contype_3 + year + conflict_id + 
    cluster(ifeconint_treat),
  data = pubgoods
)
summary(cox_mod_full)

## Create list of private goods models for comparison
models_list <- list(
  cox_mod_incl,  # Model with explicit reintegration
  cox_mod_funds, # Model with direct funds
  cox_mod_full,  # Model with training programs
  cox_model_with_pca # Model with principal component
)

## Extract coefficients for visualization
coefficients <- lapply(models_list, tidy)

# Select key coefficients of interest
selected_coefficients <- lapply(coefficients, function(model_coef) {
  model_coef[model_coef$term %in% c("funds_dummy_vSum22", "inclusive_dummy_vSum22", "tp_all_vSum22", "PC1"), ]
})

# Combine and format for plotting
selected_coefficients_df <- bind_rows(selected_coefficients, .id = "Model")

# Rename variables for clearer interpretation
selected_coefficients_df$Variable <- factor(
  selected_coefficients_df$term,
  levels = c("inclusive_dummy_vSum22", "funds_dummy_vSum22", "tp_all_vSum22", "PC1"),
  labels = c("Explicit Reintegration for Ex-Combatants", 
             "Direct Funds for Ex-Combatants", 
             "Training Program for Ex-Combatants", 
             "Economic Reintegration Principal Component")
)

## Create coefficient plot
plot_full <- ggplot(selected_coefficients_df, 
                    aes(x = Model, y = estimate, color = Variable, shape = Variable)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  labs(title = "", x = "Private Goods Models", y = "Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(family = "Times New Roman")) +
  scale_y_continuous(limits = c(-1, .5))

# Display plot
print(plot_full)

# Create detailed regression table
stargazer(models_list, type = "text", title = "Cox Regression Models - Private Goods", header = FALSE)

###############################################################
# PART 2: PUBLIC GOODS MODELS
###############################################################

# Clear environment for public goods analysis
rm(list = ls())

# Reload data
privgoods <- read_csv("privgoods.csv")

# Remove any existing principal components
privgoods <- privgoods[, !(names(privgoods) %in% c("PC1", "PC2", "PC3", "PC4", "PC5"))]

# Handle missing values in public goods variables
privgoods <- privgoods %>%
  mutate(across(c("EpsFis", "Dev", "NEC"), ~ ifelse(is.na(.), 0, .)))

## Principal Component Analysis for Public Goods
selected_columns <- c("EpsFis", "Dev", "NEC")
covariates_for_pca <- privgoods[, selected_columns]

# Perform PCA
pca_result <- prcomp(covariates_for_pca, scale. = TRUE)

# Plot variance explained
plot(cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2)), 
     xlab = "Number of Principal Components",
     ylab = "Cumulative Proportion of Variance Explained", 
     type = "b")

# Add reference line at 90% cumulative variance
abline(h = 0.9, col = "red", lty = 2)

# Display factor loadings
factor_loadings <- pca_result$rotation
print(factor_loadings)

# Extract principal components and combine with original data
pca_components <- predict(pca_result, newdata = covariates_for_pca)
privgoods <- cbind(privgoods, pca_components)

## MODEL 1: Public Goods Principal Component
cox_pub_pca <- coxph(
  Surv(duration_in_days, censored) ~
    PC1 + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id +
    cluster(ifeconint_treat),
  data = privgoods
)
summary(cox_pub_pca)

## MODEL 2: Fiscal Federalism
cox_EpsFis <- coxph(
  Surv(duration_in_days, censored) ~
    EpsFis + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id +
    cluster(ifeconint_treat),
  data = privgoods
)
summary(cox_EpsFis)

## MODEL 3: National Economic Plan
cox_NEC <- coxph(
  Surv(duration_in_days, censored) ~
    NEC + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id +
    cluster(ifeconint_treat),
  data = privgoods
)
summary(cox_NEC)

## MODEL 4: Development - Social Sector
cox_DevInfra <- coxph(
  Surv(duration_in_days, censored) ~
    DevSoc + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id +
    cluster(ifeconint_treat),
  data = privgoods
)
summary(cox_DevInfra)

## Create data frame for coefficient plot
coef_data <- data.frame(
  Model = c("Fiscal Federalism", "National Economic Plan", 
            "Development and Reconstruction", "Public Goods Principal Component"),
  estimate = c(-0.4205, -0.2474, 0.09578, 0.01603),
  robust_se = c(0.0588, 0.6209, 0.1002, 0.09228),
  p_value = c(8.58e-13, 0.690341, 0.338966, 0.862114)
)

## Create coefficient plot for public goods models
plot_coef <- ggplot(coef_data, aes(x = Model, y = estimate, color = Model, shape = Model)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = estimate - robust_se, ymax = estimate + robust_se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  labs(title = "", x = "Public Goods Models", y = "Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        text = element_text(family = "Times New Roman"),
        axis.title.x = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(-1, 0.5))

# Display plot
print(plot_coef)

## Create list of public goods models for comparison
models_list_pub <- list(
  cox_EpsFis,    # Fiscal federalism model
  cox_NEC,       # National economic plan model
  cox_DevInfra,  # Development/reconstruction model
  cox_pub_pca    # Principal component model
)

# Create detailed regression table
stargazer(models_list_pub, type = "text", title = "Cox Regression Models - Public Goods", header = FALSE)
