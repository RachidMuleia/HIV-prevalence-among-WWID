# Table 2. Stratified RDS-weighted prevalence (HIV, HBsAg, HCV),
# Rao-Scott chi-square tests, and HIV crude/adjusted prevalence ratios.

source("analysis/01_prepare_data.R") 
source("R/01_functions.R")


library(survey)
library(geepack)

covariate_levels <- function(data, vars) {
  lapply(stats::setNames(vars, vars), function(v) {
    x <- data[[v]]
    if (is.factor(x)) {
      levels(x)
    } else {
      sort(unique(as.character(x[!is.na(x)])))
    }
  })
}

# --- Weighted stratified prevalence (Gile SS among women, then REML pool) ---
run_stratified <- function(outcome, label) {
  message("Stratified RDS estimates: ", label)
  if (!exists("cov_levels", inherits = TRUE) || is.null(cov_levels)) {
    cov_levels <- covariate_levels(pid, covariates)
  }
  site_est <- rds_prov_all(
    data_list = data_list_rds,
    var_vector = covariates,
    var_outcome = outcome,
    pop_size = pop_size,
    dec = 1,
    dec.i = 1,
    site = "SITE_CITY",
    n_boot = n_boot,
    level_list = cov_levels
  )
  assign(paste0(tolower(gsub("[^A-Za-z0-9]+", "_", label)), "_site_est"),
         site_est, envir = .GlobalEnv)
  pooled <- meta_analyze_rds(site_est)
  pooled$outcome <- label
  pooled
}
print("======================HIV stratified prevalence====================================")

hiv_strat <- run_stratified("RESUL_HIV_PREV", "HIV")
print("======================HBsAg stratified prevalence====================================")
hbv_strat <- run_stratified("RESUL_HBV_CAT", "HBsAg")
hbv_strat
print("======================HCV stratified prevalence====================================")
hcv_strat  <- run_stratified("RESUL_HCV_CAT", "HCV")
hcv_strat
rio::export(hiv_strat, file.path(output_dir, "table2_hiv_weighted.xlsx"))
rio::export(hbv_strat, file.path(output_dir, "table2_hbsag_weighted.xlsx"))
rio::export(hcv_strat, file.path(output_dir, "table2_hcv_weighted.xlsx"))

# --- Rao-Scott chi-square among women (recruiter clusters, city strata) ---
wwid_df <- as.data.frame(wwid)
if (!"weights" %in% names(wwid_df)) {
  stop("Column 'weights' is required for the Rao-Scott tests.")
}

wid_design <- svydesign(
  id = ~recruiter.id,
  strata = ~SITE_CITY,
  data = pid,
  weights = ~weights,
  nest = TRUE
)
wid_design_women <- subset(wid_design, SEX_CAT == "2_FEMENINO")

rao_hiv <- cross_svy(wwid_df, covariates, "RESUL_HIV_PREV", wid_design_women)
print("======================Rao-Scott chi-square among women (recruiter clusters, city strata) for HIV====================================")
rao_hiv[,c("var", "p_value")]
rao_hbv <- cross_svy(wwid_df, covariates, "RESUL_HBV_CAT", wid_design_women)
print("======================Rao-Scott chi-square among women (recruiter clusters, city strata) for HBsAg====================================")
rao_hbv[,c("var", "p_value")]
rao_hcv <- cross_svy(wwid_df, covariates, "RESUL_HCV_CAT", wid_design_women)
print("======================Rao-Scott chi-square among women (recruiter clusters, city strata) for HCV====================================")
rao_hcv[,c("var", "p_value")]

rio::export(rao_hiv, file.path(output_dir, "table2_rao_scott_hiv.xlsx"))
rio::export(rao_hbv, file.path(output_dir, "table2_rao_scott_hbsag.xlsx"))
rio::export(rao_hcv, file.path(output_dir, "table2_rao_scott_hcv.xlsx"))

# Univariable models: Crude prevalence ratios: modified Poisson GEE, recruiter clusters
wwid_reg <- wwid_df[!is.na(wwid_df$response), ]

gee_poisson_pr <- function(dat, rhs) {
  vars <- all.vars(as.formula(paste("~", rhs)))
  dat <- droplevels(na.omit(dat[, c("response", "recruiter.id", vars)]))
  dat$recruiter.id <- droplevels(as.factor(dat$recruiter.id))
  fit <- geeglm(
    as.formula(paste("response ~", rhs)),
    family = poisson(link = "log"),
    id = recruiter.id,
    corstr = "exchangeable",
    data = dat
  )
  sm <- as.matrix(coef(summary(fit)))
  colnames(sm)[1:4] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  list(tab = extract_pr_glm(sm), fit = fit)
}

cpr <- lapply(covariates, function(v) {
  out <- gee_poisson_pr(wwid_reg, v)$tab
  out$var <- v
  out
})
cpr_tab <- do.call(rbind, cpr)
rio::export(cpr_tab, file.path(output_dir, "table2_hiv_cPR.xlsx"))



# selected covariates to be used in the multivariable models
var_wid2 <- c(
  "AGE_CAT", "SITE_CITY", "EDU_CAT1", "DEMARSTA_CAT", "DEACT",
  "MSEX6_CAT", "DRUG6_BIN", "ANOS_INJECT_CAT", "FREQ_INJ",
  "PEEREDU1_CAT", "GRAVIDEZ", "SINTOMAS_ITS", "STGXCLP_CAT"
)

# Multivariable models: Adjusted prevalence ratios: modified Poisson GEE, recruiter clusters
gee_hiv <- gee_poisson_pr(wwid_reg, paste(var_wid2, collapse = "+"))
apr_tab <- gee_hiv$tab
gee_hiv <- gee_hiv$fit
rio::export(apr_tab, file.path(output_dir, "table2_hiv_aPR_gee.xlsx"))
print(summary(gee_hiv))




# This follows the approach by Fox and Monette in the paper "Generalized Collinearity Diagnostics"
if (!requireNamespace("car", quietly = TRUE)) {
  stop("Please install.packages(\"car\") to compute VIF.")
}
gee_dat <- droplevels(na.omit(wwid_reg[, c("response", var_wid2)]))
vif_glm <- glm(
  as.formula(paste("response ~", paste(var_wid2, collapse = "+"))),
  family = poisson(link = "log"),
  data = gee_dat
)
vif_raw <- car::vif(vif_glm)
if (is.matrix(vif_raw)) {
  vif_tab <- data.frame(
    variable = rownames(vif_raw),
    GVIF = round(vif_raw[, "GVIF"], 2),
    Df = as.integer(vif_raw[, "Df"]),
    adj_GVIF = round(vif_raw[, "GVIF^(1/(2*Df))"], 2),
    stringsAsFactors = FALSE
  )
} else {
  vif_tab <- data.frame(
    variable = names(vif_raw),
    VIF = round(as.numeric(vif_raw), 2),
    stringsAsFactors = FALSE
  )
}
print(vif_tab)
rio::export(vif_tab, file.path(output_dir, "table2_hiv_vif.xlsx"))

# LRT for a recruiter random intercept vs the same mean model without it.
# 50:50 mixture of chi-square (Self and Liang 1987; Stram and Lee 1994).
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Please install.packages(\"lme4\") for the random-intercept LRT.")
}
re_dat <- droplevels(na.omit(wwid_reg[, c("response", "recruiter.id", var_wid2)]))
re_dat$recruiter.id <- droplevels(as.factor(re_dat$recruiter.id))
eq_mean <- as.formula(paste("response ~", paste(var_wid2, collapse = "+")))
m0_pois <- glm(eq_mean, family = poisson(link = "log"), data = re_dat)
m1_pois <- lme4::glmer(
  stats::update(eq_mean, . ~ . + (1 | recruiter.id)),
  family = poisson(link = "log"),
  data = re_dat
)
re_lrt <- lrt_boundary_mixture(m0_pois, m1_pois)
print(lme4::VarCorr(m1_pois))
print(re_lrt)
rio::export(re_lrt, file.path(output_dir, "table2_hiv_re_lrt_mixture.xlsx"))
