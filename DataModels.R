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

# Load data
PAFull_pmm_pc <- read_csv("PAFull_pmm_pc.csv")
padata <- read.csv("PA_Coded.csv")

###############################################################
# DATA PREPARATION
###############################################################

# Create a binary variable for agency list
padata$agency_binary <- ifelse(
  !is.na(trimws(padata$agency_list_vSum22)) & trimws(padata$agency_list_vSum22) != "0" & trimws(padata$agency_list_vSum22) != "",
  1,
  0
)

# Filling NAs
padata <- padata %>%
  mutate(
    tp_excombatant_vSum22 = ifelse(is.na(tp_excombatant_vSum22), 0, tp_excombatant_vSum22),
    economic_dummy_vSum22 = ifelse(is.na(economic_dummy_vSum22), 0, economic_dummy_vSum22),
    funds_dummy_vSum22 = ifelse(is.na(funds_dummy_vSum22), 0, funds_dummy_vSum22),
    inclusive_dummy_vSum22 = ifelse(is.na(inclusive_dummy_vSum22), 0, inclusive_dummy_vSum22), 
    tp_security_vSum22 = ifelse(is.na(tp_security_vSum22), 0, tp_security_vSum22), 
    tp_all_vSum22 = ifelse(is.na(tp_all_vSum22), 0, tp_all_vSum22), 
    tp_categorical_vSum22 = ifelse(is.na(tp_categorical_vSum22), 0, tp_categorical_vSum22)
  )

# Subset to only data that I care about
columnkeep <- c(
  "region", "paid", "pa_name", "conflict_id", "conflict_name", 
  "actor_id", "actor_name",  "year", 
  "pa_date", "duration", "ended", "ddr", "locgov", "justice_prov", "pko", 
  "inclusive", "economic_dummy_vSum22", "funds_dummy_vSum22", 
  "inclusive_dummy_vSum22", "tp_excombatant_vSum22", "agency_binary", "tp_security_vSum22", "tp_all_vSum22", "tp_categorical_vSum22"
)
PAsubset <- padata[, columnkeep]

# Converting dates to survival format
PAsubset$pa_date <- as.Date(PAsubset$pa_date)
PAsubset$duration <- as.Date(PAsubset$duration)

PAsubset$censored <- ifelse(PAsubset$ended, 1, 0)
PAsubset$censored[is.na(PAsubset$censored)] <- 0

PAsubset$duration_in_days <- as.numeric(PAsubset$duration - PAsubset$pa_date)

distinct_conflict_ids <- unique(PAsubset$conflict_id)
print(distinct_conflict_ids)

PAsubset <- PAsubset %>%
  separate_rows(conflict_id, sep = ",\\s*")
PAsubset$conflict_id <- as.integer(PAsubset$conflict_id)

######################## - Merging Data -########################

# UCDP - https://ucdp.uu.se/downloads/ucdpprio/ucdp-prio-acd-231.pdf
# variables - conflict duration, type, and intensity
ucdp <- read.csv("ucdp.csv")
str(ucdp)

# Rename variables
names(ucdp)[names(ucdp) == "start_date"] <- "conflict_startdate"
names(ucdp)[names(ucdp) == "ep_end_date"] <- "conflict_enddate"

# Convert startdate and enddate to Date format if not already
ucdp$conflict_startdate <- as.Date(ucdp$conflict_startdate, format = "%m/%d/%Y")
ucdp$conflict_enddate <- as.Date(ucdp$conflict_enddate, format = "%m/%d/%Y")
str(ucdp)

# Create a new dataset with the earliest start dates
earliest_start_dates <- ucdp %>%
  group_by(conflict_id) %>%
  summarise(conflict_start_date = min(conflict_startdate))

# Create a new dataset with the latest end dates
latest_end_dates <- ucdp %>%
  group_by(conflict_id) %>%
  summarise(conflict_end_date = max(conflict_enddate))

latest_end_dates$totoday <- as.Date("2023-11-17")  # Replace with your actual "totoday" value

latest_end_dates <- latest_end_dates %>%
  mutate(conflict_end_date = coalesce(conflict_end_date, as.Date("2023-11-17")))  # Replace with your actual "totoday" value

conflictdates <- left_join(earliest_start_dates, latest_end_dates, by = "conflict_id")

conflictdates$duration_conflict <- as.numeric(conflictdates$conflict_end_date - conflictdates$conflict_start_date)

UCDPsubset <- ucdp[, c("conflict_id", "year", "conflict_startdate", "conflict_enddate", "intensity_level", "cumulative_intensity", "type_of_conflict")]

UCDP_last_year <- UCDPsubset %>%
  group_by(conflict_id) %>%
  filter(year == max(year))

UCDP_last_year <- UCDP_last_year[, c("conflict_id", "intensity_level", "cumulative_intensity", "type_of_conflict")]

# MERGE UCDP w/ PA
PA_UCDP <- left_join(PAsubset, conflictdates, by = "conflict_id")
PA_UCDP <- left_join(PA_UCDP, UCDP_last_year, by = "conflict_id")

# Checking missingness
missing_values <- colSums(is.na(PA_UCDP))
print(missing_values[missing_values > 0])

# World Bank data - GPD, IMF, Refugees, Aid
wb <- read.csv("WorldBank.csv")
str(wb)

# Additional PA Info - Country code identifiers, Previous Agreements, 1995/2005 dummies, Agreement Page Length
addPA <- read.csv("HazardData.csv")
str(addPA)

# Merge addPA with full PA data (PAUCDP)
PAccode <- full_join(PA_UCDP, addPA, by = c("paid", "year"))

PAccode <- PAccode %>%
  rename(ccode = ccode2)

# World Bank join
wb2 <- read.csv("wb2.csv")
PA_WB <- left_join(PAccode, wb2, by = c("ccode", "year")) 
PA_WB <- left_join(PA_WB, wb, by = c("ccode", "year")) 

missing_values <- colSums(is.na(PA_WB))
print(missing_values[missing_values > 0])

# VDEM - Pull in regime variables
vdem <- read.csv("VDEM.csv")
str(vdem)

VDEMsubset <- vdem[vdem$year > 1900, c("country_text_id", "year", "COWcode", "v2x_polyarchy", "v2x_libdem", "v2x_delibdem", "v2x_egaldem")]

VDEMsubset <- VDEMsubset %>%
  rename(ccode = country_text_id)

# VDEM Join 
PAFull <- left_join(PA_WB, VDEMsubset, by = c("ccode", "year")) 

summary_stats <- summary(PAFull)
print(summary_stats)

# Adjusting for censored data
PAFull$totoday <- as.Date("2023-11-17")  # Replace with your actual "totoday" value
difference_days <- as.numeric(PAFull$totoday - PAFull$pa_date)
na_indices <- is.na(PAFull$duration_in_days)
PAFull$duration_in_days[na_indices] <- difference_days[na_indices]

# Check missingness
missing_values <- colSums(is.na(PAFull))
print(missing_values)

# Specify the variables of interest
vars_of_interest <- c(
  "funds_dummy_vSum22", 
  "RefugeePop_origin", "GDP_PPP", "RefugeePop_asylum", "develassistance",
  "develassistance_aid", "Lgt","v2x_polyarchy","post_2005"
)

# Subset the data
PAFull_nonmissing <- PAFull[complete.cases(PAFull[, vars_of_interest]), ]

# Multiple Imputations - PMM
# Compute the correlation matrix
numeric_data <- PAFull[sapply(PAFull, is.numeric)]
cor_matrix <- cor(numeric_data, use = "complete.obs")
print(cor_matrix)

numeric_data <- PAFull[, sapply(PAFull, is.numeric)]

# Calculate the correlation matrix
cor_matrix <- cor(numeric_data)
# Set a threshold for correlation
threshold <- 0.5  # You can adjust this value based on your needs

# Find highly correlated variable pairs
highly_correlated <- which(upper.tri(cor_matrix, diag = TRUE) & abs(cor_matrix) > threshold, arr.ind = TRUE)

# Print the pairs of highly correlated variables
for (i in 1:nrow(highly_correlated)) {
  var1 <- rownames(cor_matrix)[highly_correlated[i, 1]]
  var2 <- colnames(cor_matrix)[highly_correlated[i, 2]]
  correlation <- cor_matrix[highly_correlated[i, 1], highly_correlated[i, 2]]
  
  cat(sprintf("Variables %s and %s are highly correlated with a correlation of %.2f\n", var1, var2, correlation))
}

# Columns to exclude from imputation
columns_to_exclude <- c("paid", "conflict_id", "year", "censored", "region", "censored", "duration_in_days", "COWcode")

columns_for_imputation <- setdiff(c("RefugeePop_origin", "GDP_PPP", "RefugeePop_asylum", "develassistance",
                                    "develassistance_aid", "Lgt", "v2x_polyarchy", "v2x_libdem","v2x_delibdem",
                                    "v2x_egaldem", "IMF", "corruption_percentile", "ROL_percentile", "ROL_estimate", 
                                    "voiceaccount_sources"), columns_to_exclude)

# Scale only the columns for imputation
scaled_data_for_imputation <- scale(PAFull[, columns_for_imputation])

# Combine scaled imputed data with non-imputed data
PAFull_pmm <- cbind(PAFull[, setdiff(names(PAFull), columns_for_imputation)], 
                    complete(mice(scaled_data_for_imputation, method = "pmm", m = 10)))

# Check missing values
missing_values <- colSums(is.na(PAFull_pmm))
print(missing_values)

# Check normality
check_normality <- function(variable) {
  shapiro_result <- shapiro.test(variable)
  ad_result <- ad.test(variable)
  
  return(list(
    Variable = deparse(substitute(variable)),
    ShapiroWilk_P_Value = shapiro_result$p.value,
    AndersonDarling_P_Value = ad_result$p.value
  ))
}

numeric_vars <- PAFull_pmm %>% select_if(is.numeric)
normality_results <- map(numeric_vars, check_normality)

# Filter variables with non-normal distributions
non_normal_vars <- Filter(function(x) {
  x$ShapiroWilk_P_Value < 0.05 | x$AndersonDarling_P_Value < 0.05
}, normality_results)

cat("Non-normal variables:\n")
if (length(non_normal_vars) == 0) {
  cat("None\n")
} else {
  for (result in non_normal_vars) {
    cat("Variable:", result$Variable, "\n")
    cat("  Shapiro-Wilk p-value:", result$ShapiroWilk_P_Value, "\n")
    cat("  Anderson-Darling p-value:", result$AndersonDarling_P_Value, "\n\n")
  }
}

# Log-transform variables
# Add a small constant value to handle zero values
constant_value <- 1e-6  

# Log-transform the specified variables
PAFull$log_PrevAgmt <- log(PAFull$PrevAgmt + constant_value)
PAFull$log_RefugeePop_origin <- log(PAFull$RefugeePop_origin + constant_value)
PAFull$log_GDP_PPP <- log(PAFull$GDP_PPP + constant_value)
PAFull$log_RefugeePop_asylum <- log(PAFull$RefugeePop_asylum + constant_value)
PAFull$log_develassistance <- log(PAFull$develassistance + constant_value)
PAFull$log_develassistance_aid <- log(PAFull$develassistance_aid + constant_value)
PAFull$log_Lgt <- log(PAFull$Lgt + constant_value)
PAFull$log_IMF <- log(PAFull$IMF + constant_value)
PAFull$log_v2x_polyarchy <- log(PAFull$v2x_polyarchy + constant_value)
PAFull$log_v2x_libdem <- log(PAFull$v2x_libdem + constant_value)
PAFull$log_v2x_delibdem <- log(PAFull$v2x_delibdem + constant_value)
PAFull$log_v2x_egaldem <- log(PAFull$v2x_egaldem + constant_value)

PAFull_pmm$log_PrevAgmt <- log(ifelse(PAFull_pmm$PrevAgmt <= 0, constant_value, PAFull_pmm$PrevAgmt))
PAFull_pmm$log_RefugeePop_origin <- log(ifelse(PAFull_pmm$RefugeePop_origin <= 0, constant_value, PAFull_pmm$RefugeePop_origin))
PAFull_pmm$log_GDP_PPP <- log(ifelse(PAFull_pmm$GDP_PPP <= 0, constant_value, PAFull_pmm$GDP_PPP))
PAFull_pmm$log_RefugeePop_asylum <- log(ifelse(PAFull_pmm$RefugeePop_asylum <= 0, constant_value, PAFull_pmm$RefugeePop_asylum))
PAFull_pmm$log_develassistance <- log(ifelse(PAFull_pmm$develassistance <= 0, constant_value, PAFull_pmm$develassistance))
PAFull_pmm$log_develassistance_aid <- log(ifelse(PAFull_pmm$develassistance_aid <= 0, constant_value, PAFull_pmm$develassistance_aid))
PAFull_pmm$log_Lgt <- log(ifelse(PAFull_pmm$Lgt <= 0, constant_value, PAFull_pmm$Lgt))
PAFull_pmm$log_IMF <- log(ifelse(PAFull_pmm$IMF <= 0, constant_value, PAFull_pmm$IMF))
PAFull_pmm$log_v2x_polyarchy <- log(ifelse(PAFull_pmm$v2x_polyarchy <= 0, constant_value, PAFull_pmm$v2x_polyarchy))
PAFull_pmm$log_v2x_libdem <- log(ifelse(PAFull_pmm$v2x_libdem <= 0, constant_value, PAFull_pmm$v2x_libdem))
PAFull_pmm$log_v2x_delibdem <- log(ifelse(PAFull_pmm$v2x_delibdem <= 0, constant_value, PAFull_pmm$v2x_delibdem))
PAFull_pmm$log_v2x_egaldem <- log(ifelse(PAFull_pmm$v2x_egaldem <= 0, constant_value, PAFull_pmm$v2x_egaldem))
PAFull_pmm$log_corruption_percentile <- log(ifelse(PAFull_pmm$corruption_percentile <= 0, constant_value, PAFull_pmm$corruption_percentile))
PAFull_pmm$log_ROL_percentile <- log(ifelse(PAFull_pmm$ROL_percentile <= 0, constant_value, PAFull_pmm$ROL_percentile))
PAFull_pmm$log_ROL_estimate <- log(ifelse(PAFull_pmm$ROL_estimate <= 0, constant_value, PAFull_pmm$ROL_estimate))
PAFull_pmm$log_voiceaccount_sources <- log(ifelse(PAFull_pmm$voiceaccount_sources <= 0, constant_value, PAFull_pmm$voiceaccount_sources))

PAFull_nonmissing <- PAFull_nonmissing %>%
  mutate(
    log_PrevAgmt = log(ifelse(PrevAgmt <= 0, constant_value, PrevAgmt)),
    log_RefugeePop_origin = log(ifelse(RefugeePop_origin <= 0, constant_value, RefugeePop_origin)),
    log_GDP_PPP = log(ifelse(GDP_PPP <= 0, constant_value, GDP_PPP)),
    log_RefugeePop_asylum = log(ifelse(RefugeePop_asylum <= 0, constant_value, RefugeePop_asylum)),
    log_develassistance = log(ifelse(develassistance <= 0, constant_value, develassistance)),
    log_develassistance_aid = log(ifelse(develassistance_aid <= 0, constant_value, develassistance_aid)),
    log_Lgt = log(ifelse(Lgt <= 0, constant_value, Lgt)),
    log_IMF = log(ifelse(IMF <= 0, constant_value, IMF)),
    log_v2x_polyarchy = log(ifelse(v2x_polyarchy <= 0, constant_value, v2x_polyarchy)),
    log_v2x_libdem = log(ifelse(v2x_libdem <= 0, constant_value, v2x_libdem)),
    log_v2x_delibdem = log(ifelse(v2x_delibdem <= 0, constant_value, v2x_delibdem)),
    log_v2x_egaldem = log(ifelse(v2x_egaldem <= 0, constant_value, v2x_egaldem)),
    log_corruption_percentile = log(ifelse(corruption_percentile <= 0, constant_value, corruption_percentile)),
    log_ROL_percentile = log(ifelse(ROL_percentile <= 0, constant_value, ROL_percentile)),
    log_ROL_estimate = log(ifelse(ROL_estimate <= 0, constant_value, ROL_estimate)),
    log_voiceaccount_sources = log(ifelse(voiceaccount_sources <= 0, constant_value, voiceaccount_sources))
  )

# Check missing values
missing_values <- colSums(is.na(PAFull_pmm))
print(missing_values)

# Check Main IV Correlation
cor_matrix <- cor(PAFull_pmm[, c("funds_dummy_vSum22", "inclusive_dummy_vSum22", "tp_excombatant_vSum22")])
print(cor_matrix)

# Create interaction terms
PAFull_pmm$full_int <- as.numeric(PAFull_pmm$tp_excombatant_vSum22 == 1 & PAFull_pmm$funds_dummy_vSum22 == 1)
PAFull_pmm$total_int <- rowSums(PAFull_pmm[, c("funds_dummy_vSum22", "inclusive_dummy_vSum22", "tp_excombatant_vSum22", "economic_dummy_vSum22")] == 1)

cor_matrix <- cor(PAFull_pmm[, c("funds_dummy_vSum22", "inclusive_dummy_vSum22", "tp_excombatant_vSum22")])
print(cor_matrix)

# Create binary indicators for each conflict type category
PAFull_pmm$contype_3 <- as.integer(PAFull_pmm$type_of_conflict == 3)
PAFull_pmm$contype_4 <- as.integer(PAFull_pmm$type_of_conflict == 4)
PAFull_pmm$contype_2 <- as.integer(PAFull_pmm$type_of_conflict == 2)

# Cox Stepwise - Full INT
surv_obj <- with(PAFull_pmm, Surv(duration_in_days, censored))

# Create a data frame with the variables of interest
covariates <- PAFull_pmm %>%
  dplyr::select(full_int, ddr, locgov, justice_prov, pko, inclusive, intensity_level, cumulative_intensity,
                PrevAgmt, PrevAgmt_Bin, post_1990, post_2005, v2x_polyarchy, v2x_libdem, v2x_delibdem, v2x_egaldem,
                log_PrevAgmt, log_RefugeePop_origin, log_GDP_PPP, log_RefugeePop_asylum, log_develassistance,
                log_Lgt, log_IMF, v2x_polyarchy, contype_2, contype_3, contype_4, log_corruption_percentile, 
                log_ROL_percentile, log_ROL_estimate, log_voiceaccount_sources)

covariates <- cbind(covariates, duration_in_days = PAFull_pmm$duration_in_days, censored = PAFull_pmm$censored)
cox_model <- coxph(Surv(duration_in_days, censored) ~ full_int + ., data = covariates)
summary(cox_model)

# Variable selection
surv_obj <- with(covariates, Surv(duration_in_days, censored))
step_model <- stepAIC(coxph(Surv(duration_in_days, censored) ~ full_int + ., data = covariates), direction = "both")
summary(step_model)

# Treatment variable creation
PAFull_pmm <- PAFull_pmm %>%
  mutate(ifeconint_treat = ifelse(funds_dummy_vSum22 != 0 | tp_excombatant_vSum22 != 0 | inclusive_dummy_vSum22 != 0, 1, 0))

PAFull_nonmissing <- PAFull_nonmissing %>%
  mutate(ifeconint_treat = ifelse(funds_dummy_vSum22 != 0 | tp_excombatant_vSum22 != 0 | inclusive_dummy_vSum22 != 0, 1, 0))

# Execute Matching
matching_covariates <- c("ddr", "locgov", "justice_prov", "pko", "inclusive", "intensity_level",
                         "cumulative_intensity", "PrevAgmt", "PrevAgmt_Bin", "post_1990", "post_2005",
                         "v2x_polyarchy", "v2x_libdem", "v2x_delibdem", "v2x_egaldem", "log_PrevAgmt",
                         "log_RefugeePop_origin", "log_GDP_PPP", "log_RefugeePop_asylum", "log_develassistance",
                         "log_Lgt", "log_IMF", "v2x_polyarchy", "contype_2", "contype_3", "contype_4",
                         "log_corruption_percentile", "log_ROL_percentile", "log_ROL_estimate", "log_voiceaccount_sources")

formula_matching <- as.formula(paste("full_int ~", paste(matching_covariates, collapse = " + ")))
matched_data <- matchit(formula_matching, data = PAFull_pmm, method = "nearest")

# Extract matched data
matched_data <- match.data(matched_data)

# Convert year and conflict_id to factors
PAFull_pmm$year <- as.factor(PAFull_pmm$year)
PAFull_pmm$conflict_id <- as.factor(PAFull_pmm$conflict_id)

matched_data$year <- as.factor(matched_data$year)
matched_data$conflict_id <- as.factor(matched_data$conflict_id)


###############################################################
# PART 1: PRIVATE GOODS MODELS
###############################################################

## Principal Component Analysis for Private Goods
# Select variables for PCA
selected_columns <- c("funds_dummy_vSum22", "inclusive_dummy_vSum22", "tp_all_vSum22")
covariates_for_pca <- PAFull_pmm[, selected_columns]

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
PAFull_pmm_pc <- cbind(PAFull_pmm, pca_components)

## MODEL 1: Cox model with Principal Component
cox_model_with_pca <- coxph(
  Surv(duration_in_days, censored) ~ 
    PC1 + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id + 
    cluster(ifeconint_treat),
  data = PAFull_pmm_pc
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
  data = PAFull_pmm
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
  data = PAFull_pmm
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
  data = PAFull_pmm
)
summary(cox_mod_econ)

## MODEL 5: Training Programs for Ex-Combatants
cox_mod_full <- coxph(
  Surv(duration_in_days, censored) ~
    tp_all_vSum22 + log_Lgt + cumulative_intensity + log_GDP_PPP + 
    log_develassistance + v2x_libdem + PrevAgmt_Bin + post_1990 + 
    post_2005 + ddr + contype_2 + contype_3 + year + conflict_id + 
    cluster(ifeconint_treat),
  data = PAFull_pmm
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
PAFull_pmm_pc <- read_csv("PAFull_pmm_pc.csv")

# Remove any existing principal components
PAFull_pmm_pc <- PAFull_pmm_pc[, !(names(PAFull_pmm_pc) %in% c("PC1", "PC2", "PC3", "PC4", "PC5"))]

# Handle missing values in public goods variables
PAFull_pmm_pc <- PAFull_pmm_pc %>%
  mutate(across(c("EpsFis", "Dev", "NEC"), ~ ifelse(is.na(.), 0, .)))

## Principal Component Analysis for Public Goods
selected_columns <- c("EpsFis", "Dev", "NEC")
covariates_for_pca <- PAFull_pmm_pc[, selected_columns]

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
PAFull_pmm_pc <- cbind(PAFull_pmm_pc, pca_components)

## MODEL 1: Public Goods Principal Component
cox_pub_pca <- coxph(
  Surv(duration_in_days, censored) ~
    PC1 + log_Lgt + cumulative_intensity + 
    log_GDP_PPP + log_develassistance +
    v2x_libdem + PrevAgmt_Bin + post_1990 + post_2005 + ddr +
    contype_2 + contype_3 + year + conflict_id +
    cluster(ifeconint_treat),
  data = PAFull_pmm_pc
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
  data = PAFull_pmm_pc
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
  data = PAFull_pmm_pc
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
  data = PAFull_pmm_pc
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