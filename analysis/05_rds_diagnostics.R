# RDS diagnostics referred to in the Methods (seeds, waves, sex homophily).

source("analysis/01_prepare_data.R")

diag_one <- function(rds_dat, N_site, site) {
  sex <- factor(
    ifelse(is_woman(rds_dat), "Woman", "Man"),
    levels = c("Woman", "Man")
  )
  is_seed <- as.character(rds_dat$recruiter.id) == "seed"
  wave <- RDS::get.wave(rds_dat)
  deg <- RDS::get.net.size(rds_dat)
  data.frame(
    site = site,
    n = nrow(rds_dat),
    n_women = sum(sex == "Woman"),
    n_seeds = sum(is_seed),
    n_seeds_women = sum(is_seed & sex == "Woman"),
    max_wave = max(wave, na.rm = TRUE),
    mean_degree_women = round(mean(deg[sex == "Woman"], na.rm = TRUE), 2),
    mean_degree_men = round(mean(deg[sex == "Man"], na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
}

site_names <- vapply(data_list_rds, function(d) {
  as.character(unique(d$SITE_CITY))[1]
}, character(1))

diagnostics <- do.call(rbind, lapply(seq_along(data_list_rds), function(i) {
  diag_one(data_list_rds[[i]], pop_size[i], site_names[i])
}))

print(diagnostics)
rio::export(diagnostics, file.path(output_dir, "rds_diagnostics.xlsx"))
