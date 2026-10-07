# Paths and survey constants for the WWID paper analyses.
#
# Individual-level BBS data are not stored in this repository.
# Set DATA_FILE to the local analytic file (same structure as DATA_PID_CLUSTER.csv).

data_file <- Sys.getenv(
  "WWID_DATA_FILE",
  unset = file.path("..", "DATA_PID_CLUSTER.csv")
)

output_dir <- file.path("output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# PWID population sizes by city, in the same order as
# sort(unique(SITE_CITY)): Maputo, Beira, Tete, Quelimane, Nampula.
pop_size <- c(990, 1029, 3219, 681, 1792)

# Bootstrap replicates for Gile's SS. The published tables used 1000.
# Increase to 15000 to match the longer sensitivity runs.
n_boot <- as.integer(Sys.getenv("WWID_N_BOOT", unset = "1000"))

options(survey.lonely.psu = "average")
