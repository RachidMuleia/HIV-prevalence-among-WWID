# Reproduce the published WWID tables from a local copy of DATA_PID_CLUSTER.csv.
# Run from the repository root.

pkgs <- c(
  "RDS", "metafor", "survey", "sandwich", "lmtest", "geepack",
  "dplyr", "rio"
)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install packages: ", paste(missing, collapse = ", "))
}

invisible(lapply(pkgs, library, character.only = TRUE))

source("analysis/02_table1_prevalence.R")
source("analysis/03_table2_stratified_regression.R")
source("analysis/04_cascade.R")
source("analysis/05_rds_diagnostics.R")

message("Finished. Tables written to output/.")
