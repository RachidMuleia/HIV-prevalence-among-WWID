# Table 1. Crude and pooled RDS-weighted prevalence of HIV, HBsAg, HCV, and co-infections.

source("analysis/01_prepare_data.R")

outcomes <- c(
  HIV = "RESUL_HIV_PREV",
  HBsAg = "RESUL_HBV_CAT",
  HCV = "RESUL_HCV_CAT",
  HIV_HBsAg = "COINF_HIV_HBV",
  HIV_HCV = "COINF_HIV_HBC",
  HIV_HBsAg_HCV = "COINF_HIV_HBCV"
)

site_names <- vapply(data_list_rds, function(d) {
  as.character(unique(d$SITE_CITY))[1]
}, character(1))

table1_site <- do.call(rbind, lapply(names(outcomes), function(lab) {
  ov <- outcomes[[lab]]
  do.call(rbind, lapply(seq_along(data_list_rds), function(i) {
    rds_dat <- data_list_rds[[i]]
    woman <- is_woman(rds_dat)
    message(lab, " / ", site_names[i])
    fit <- gile_domain(rds_dat, ov, pop_size[i], woman, n_boot)
    pos <- grepl("POSITIVO|1_SIM", as.character(rds_dat[[ov]]))
    n_pos <- sum(pos[woman], na.rm = TRUE)
    n_tot <- sum(!is.na(rds_dat[[ov]][woman]))
    data.frame(
      outcome = lab,
      site = site_names[i],
      n_pos = n_pos,
      n = n_tot,
      crude_pct = round(100 * n_pos / n_tot, 1),
      est = fit$est,
      se = fit$se,
      stringsAsFactors = FALSE
    )
  }))
}))

table1 <- do.call(rbind, lapply(names(outcomes), function(lab) {
  d <- table1_site[table1_site$outcome == lab, ]
  p <- pool_sites_reml(d$est, d$se)
  n_pos <- sum(d$n_pos)
  n <- sum(d$n)
  data.frame(
    outcome = lab,
    n_N = paste0(n_pos, "/", n),
    crude_pct = round(100 * n_pos / n, 1),
    weighted_pct = round(p$pct, 1),
    ci = paste0(round(p$lwr, 1), "-", round(p$upr, 1)),
    I2 = round(p$I2, 1),
    stringsAsFactors = FALSE
  )
}))

print(table1)
rio::export(table1, file.path(output_dir, "table1_overall_prevalence.xlsx"))
rio::export(table1_site, file.path(output_dir, "table1_by_city.xlsx"))
