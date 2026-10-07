# Functions used for Tables 1–2 and the HIV cascade.

positive_level <- function(estimate_vec) {
  nm <- names(estimate_vec)
  if (!is.null(nm) && any(grepl("POSITIVO|1_SIM", nm))) {
    which(grepl("POSITIVO|1_SIM", nm))[1]
  } else {
    1L
  }
}

extract_positive <- function(rds_fit) {
  if (is.null(rds_fit) || (length(rds_fit) == 1 && is.na(rds_fit))) {
    return(list(est = NA_real_, se = NA_real_))
  }
  est <- rds_fit$estimate
  se_vec <- tryCatch(
    as.numeric(attr(rds_fit$interval, "bsresult")$se_estimate),
    error = function(e) as.numeric(attr(rds_fit, "bsresult")$se_estimate)
  )
  idx <- positive_level(est)
  list(est = as.numeric(est[idx]), se = as.numeric(se_vec[idx]))
}

# Original Table 2 extraction: first RDS estimate and its bootstrap SE.
extract_rds_first <- function(rds_fit) {
  if (is.null(rds_fit) || (length(rds_fit) == 1 && is.na(rds_fit))) {
    return(list(est = NA_real_, se = NA_real_))
  }
  se <- tryCatch(
    as.numeric(attr(rds_fit, "bsresult")$se_estimate[1]),
    error = function(e) NA_real_
  )
  list(est = as.numeric(rds_fit$estimate[1]), se = se)
}

is_woman <- function(data) {
  !is.na(data$SEX_CAT) & as.character(data$SEX_CAT) == "2_FEMENINO"
}

rdsat_format <- function(data, site_list, network_size, recruiter_ide, var_site, ...) {
  data_city <- lapply(site_list, function(x) data[data[[var_site]] == x, ])
  lapply(data_city, function(x) {
    RDS::as.rds.data.frame(
      x,
      recruiter.id = recruiter_ide,
      network.size = network_size,
      ...
    )
  })
}

gile_domain <- function(rds_dat, outcome, N_site, subset, n_boot) {
  fit <- tryCatch(
    RDS::RDS.bootstrap.intervals(
      rds_dat,
      outcome.variable = outcome,
      weight.type = "Gile's SS",
      uncertainty = "Gile",
      confidence.level = 0.95,
      number.of.bootstrap.samples = n_boot,
      to.factor = TRUE,
      N = N_site,
      subset = subset
    ),
    error = function(e) {
      message(e$message)
      NULL
    }
  )
  extract_positive(fit)
}

extract_lab <- function(data, variable) {
  if (!variable %in% names(data)) {
    stop(paste("Variable", variable, "not found in data."))
  }
  x <- data[[variable]]
  if (is.factor(x) && length(levels(x)) > 0) {
    lev <- levels(x)
    return(lev[!is.na(lev) & lev != "NA" & nzchar(lev)])
  }
  all_values <- as.character(x)
  unique(all_values[!is.na(all_values) & all_values != "NA"])
}

# Site-specific Gile SS prevalence among women, within one covariate category.
rds_prev <- function(data, x_var, var_outcome, population_size, dec, dec.i, n_boot,
                    cat_levels = NULL, ...) {
  categories <- if (!is.null(cat_levels)) {
    cat_levels
  } else {
    extract_lab(data = data, variable = x_var)
  }
  woman <- is_woman(data)

  prevalence <- lapply(categories, function(x) {
    in_cell <- woman & !is.na(data[[x_var]]) & as.character(data[[x_var]]) == as.character(x)
    n_cell <- sum(in_cell, na.rm = TRUE)
    if (n_cell < 2) {
      if (n_cell == 1) {
        message("Skipping ", x_var, " = ", x, " (n = ", n_cell, ")")
      }
      return(NA)
    }
    tryCatch(
      {
        rds_res <- suppressWarnings(
          RDS::RDS.bootstrap.intervals(
            data,
            outcome.variable = var_outcome,
            weight.type = "Gile's SS",
            uncertainty = "Gile",
            confidence.level = 0.95,
            number.of.bootstrap.samples = n_boot,
            subset = in_cell,
            N = population_size,
            ...
          )
        )
        if (is.null(rds_res) || (length(rds_res) == 1 && is.na(rds_res))) {
          return(NA)
        }
        rds_res
      },
      error = function(e) {
        message("Skipping ", x_var, " = ", x, ": ", e$message)
        NA
      }
    )
  })

  prev_raw <- vapply(prevalence, function(result) {
    extract_rds_first(result)$est
  }, numeric(1))
  se_raw <- vapply(prevalence, function(result) {
    extract_rds_first(result)$se
  }, numeric(1))

  data.frame(
    var = x_var,
    cat = categories,
    prevalence = round(prev_raw * 100, dec),
    se = round(se_raw * 100, dec.i),
    stringsAsFactors = FALSE
  )
}

rds_prev_all <- function(data, vec_var, var_outcome, level_list = NULL, ...) {
  do.call(rbind, lapply(vec_var, function(x) {
    rds_prev(
      data = data,
      x_var = x,
      var_outcome = var_outcome,
      cat_levels = if (is.null(level_list)) NULL else level_list[[x]],
      ...
    )
  }))
}

rds_prov_all <- function(data_list, var_vector, var_outcome, pop_size, dec, dec.i, site,
                        n_boot, level_list = NULL) {
  all_results_list <- lapply(seq_along(data_list), function(i) {
    data <- data_list[[i]]
    current_site_name <- unique(as.character(data[[site]]))[1]
    df <- rds_prev_all(
      data = data,
      vec_var = var_vector,
      var_outcome = var_outcome,
      population_size = pop_size[i],
      dec = dec,
      dec.i = dec.i,
      n_boot = n_boot,
      level_list = level_list
    )
    names(df)[names(df) %in% c("prevalence", "se")] <-
      paste(c("prevalence", "se"), current_site_name, sep = "_")
    df
  })

  combined_df <- Reduce(function(x, y) {
    dplyr::full_join(x, y, by = c("var", "cat"))
  }, all_results_list)

  combined_df$var_ord <- match(combined_df$var, var_vector)
  combined_df$cat_ord <- vapply(seq_len(nrow(combined_df)), function(i) {
    lev <- if (is.null(level_list)) NULL else level_list[[combined_df$var[i]]]
    if (is.null(lev)) {
      return(i)
    }
    m <- match(combined_df$cat[i], lev)
    if (is.na(m)) length(lev) + 1L else as.integer(m)
  }, integer(1))
  combined_df <- combined_df[order(combined_df$var_ord, combined_df$cat_ord, na.last = TRUE), ]
  combined_df$var_ord <- NULL
  combined_df$cat_ord <- NULL
  rownames(combined_df) <- NULL

  combined_df$SITE_CATEGORY <- paste(combined_df$var, combined_df$cat, sep = ":")
  combined_df
}

scalar_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1 || !is.finite(x)) NA_real_ else x
}

pool_rds_row <- function(study_data_clean, meta_model) {
  n <- nrow(study_data_clean)
  out <- list(
    pooled_prevalence = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    I2 = NA_real_
  )
  if (!is.null(meta_model)) {
    out$pooled_prevalence <- round(scalar_or_na(meta_model$b) * 100, 2)
    out$ci_lower <- round(scalar_or_na(meta_model$ci.lb) * 100, 2)
    out$ci_upper <- round(scalar_or_na(meta_model$ci.ub) * 100, 2)
    out$I2 <- round(scalar_or_na(meta_model$I2), 1)
    return(out)
  }
  if (n == 1) {
    yi <- study_data_clean$yi_prop[1]
    sei <- sqrt(max(study_data_clean$vi_prop[1], 0))
    out$pooled_prevalence <- round(yi * 100, 2)
    out$ci_lower <- round((yi - 1.96 * sei) * 100, 2)
    out$ci_upper <- round((yi + 1.96 * sei) * 100, 2)
  }
  out
}

# Random-effects (REML) pool of site-specific RDS estimates.
meta_analyze_rds <- function(rds_prov_all_output) {
  prev_cols <- names(rds_prov_all_output)[grepl("^prevalence_", names(rds_prov_all_output))]
  se_cols <- names(rds_prov_all_output)[grepl("^se_", names(rds_prov_all_output))]

  meta_results <- rds_prov_all_output %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      study_data_raw = list(
        data.frame(
          yi_prop = dplyr::c_across(dplyr::all_of(prev_cols)) / 100,
          vi_prop = (dplyr::c_across(dplyr::all_of(se_cols)) / 100)^2
        )
      ),
      study_data_clean = list(
        study_data_raw %>%
          dplyr::filter(!is.na(yi_prop) & !is.na(vi_prop))
      ),
      n_valid_studies = nrow(study_data_clean),
      meta_model = list(
        if (n_valid_studies >= 2) {
          tryCatch(
            metafor::rma(
              yi = yi_prop, vi = vi_prop,
              data = study_data_clean, method = "REML"
            ),
            error = function(e) NULL
          )
        } else {
          NULL
        }
      ),
      pooled = list(pool_rds_row(study_data_clean, meta_model)),
      pooled_prevalence = pooled$pooled_prevalence,
      ci_lower = pooled$ci_lower,
      ci_upper = pooled$ci_upper,
      I2 = pooled$I2,
      pooled_prevalence_ci = paste0(
        pooled_prevalence, " (", ci_lower, " - ", ci_upper, ")"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      SITE_CATEGORY, var, cat, pooled_prevalence_ci,
      pooled_prevalence, ci_lower, ci_upper, I2, n_valid_studies,
      dplyr::starts_with("prevalence_"), dplyr::starts_with("se_")
    ) %>%
    as.data.frame()

  meta_results
}

cross_svy <- function(data, var_vec, response, tsdesign) {
  list_for <- lapply(var_vec, function(x) {
    as.formula(paste("~", paste(x, response, sep = "+")))
  })
  prop_svy <- lapply(list_for, function(x) {
    round(prop.table(survey::svytable(x, design = tsdesign), 1) * 100, 1)
  })
  p_value <- lapply(list_for, function(x) {
    survey::svychisq(x, design = tsdesign)
  })
  raw_table <- lapply(var_vec, function(x) {
    table(data[, x], data[, response])
  })

  table_ls <- list()
  for (i in seq_along(var_vec)) {
    col1 <- paste(raw_table[[i]][, 1], " (", prop_svy[[i]][, 1], ")", sep = "")
    col2 <- paste(raw_table[[i]][, 2], " (", prop_svy[[i]][, 2], ")", sep = "")
    table_ls[[i]] <- data.frame(
      var = rep(var_vec[i], length(rownames(raw_table[[i]]))),
      category = rownames(raw_table[[i]]),
      col1 = col1,
      col2 = col2,
      p_value = rep(p_value[[i]]$p.value, length(rownames(raw_table[[i]]))),
      stringsAsFactors = FALSE
    )
  }
  table_descr <- do.call(rbind, table_ls)
  names(table_descr)[c(3, 4)] <- sort(na.omit(unique(as.character(data[, response]))))
  table_descr
}

extract_pr_glm <- function(x, dec_or = 2, dec = 2) {
  coeff <- x[, 1]
  lab <- attr(x, "dimnames")[[1]]
  li <- coeff - 1.96 * x[, 2]
  ls <- coeff + 1.96 * x[, 2]
  conf <- paste0("(", round(exp(li), dec), "-", round(exp(ls), dec), ")")
  data.frame(
    variable_cat = lab,
    PR_CI = paste(round(exp(coeff), dec_or), conf),
    p_value = round(x[, 4], 3),
    stringsAsFactors = FALSE
  )
}

pool_sites_reml <- function(est, se) {
  ok <- is.finite(est) & is.finite(se) & se > 0
  if (sum(ok) < 2) {
    return(list(pct = NA_real_, lwr = NA_real_, upr = NA_real_, I2 = NA_real_))
  }
  m <- metafor::rma(yi = est[ok], sei = se[ok], method = "REML")
  list(
    pct = as.numeric(m$b) * 100,
    lwr = as.numeric(m$ci.lb) * 100,
    upr = as.numeric(m$ci.ub) * 100,
    I2 = as.numeric(m$I2)
  )
}

# LRT of a variance on the boundary (H0: sigma^2 = 0).
# Null distribution is 0.5 chi^2_0 + 0.5 chi^2_1 (Self and Liang 1987; Stram and Lee 1994).
# p_mixture = 0.5 * P(chi^2_1 > LRT). p_chisq1 is the conservative chi^2_1 p-value.
lrt_boundary_mixture <- function(m0, m1) {
  ll0 <- as.numeric(stats::logLik(m0))
  ll1 <- as.numeric(stats::logLik(m1))
  lrt <- as.numeric(2 * (ll1 - ll0))
  if (!is.finite(lrt) || lrt <= 0) {
    return(data.frame(
      logLik_null = ll0,
      logLik_re = ll1,
      LRT = max(0, lrt, na.rm = TRUE),
      p_mixture = 1,
      p_chisq1 = 1,
      stringsAsFactors = FALSE
    ))
  }
  p1 <- stats::pchisq(lrt, df = 1, lower.tail = FALSE)
  data.frame(
    logLik_null = ll0,
    logLik_re = ll1,
    LRT = lrt,
    p_mixture = 0.5 * p1,
    p_chisq1 = p1,
    stringsAsFactors = FALSE
  )
}
