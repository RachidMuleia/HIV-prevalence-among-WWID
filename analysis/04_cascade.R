# UNAIDS 95-95-95 cascade among WWID living with HIV (Results / Figure).

source("analysis/01_prepare_data.R")

hiv_pos <- wwid[as.character(wwid$RESUL_HIV_PREV) == "1_POSITIVO", ]
known <- hiv_pos[!is.na(hiv_pos$CONHECIMENTO_SERO_ESTADO), ]
n_known <- sum(known$CONHECIMENTO_SERO_ESTADO == "1_CONHECE")
n_status <- nrow(known)
first95 <- n_known / n_status

on_art <- known[known$CONHECIMENTO_SERO_ESTADO == "1_CONHECE", ]
n_art <- sum(on_art$TOMA_TARV2 == "1_SIM", na.rm = TRUE)
n_aware <- nrow(on_art)
second95 <- n_art / n_aware

suppressed <- on_art[on_art$TOMA_TARV2 == "1_SIM" & !is.na(on_art$CARGA_VIRAL), ]
n_vs <- sum(suppressed$CARGA_VIRAL == "1_SUPREMIDO")
n_on_art <- nrow(suppressed)
third95 <- n_vs / n_on_art

cascade <- data.frame(
  step = c(
    "Aware of HIV-positive status (1st 95)",
    "On ART among those aware (2nd 95)",
    "Virally suppressed among those on ART (3rd 95)",
    "Cumulatively suppressed among all WWID living with HIV"
  ),
  n_N = c(
    paste0(n_known, "/", n_status),
    paste0(n_art, "/", n_aware),
    paste0(n_vs, "/", n_on_art),
    paste0(n_vs, "/", n_status)
  ),
  percent = round(100 * c(first95, second95, third95, first95 * second95 * third95), 1),
  unaids_target = c(95, 95, 95, round(100 * 0.95^3, 1)),
  stringsAsFactors = FALSE
)

print(cascade)
rio::export(cascade, file.path(output_dir, "figure_hiv_cascade_95_95_95.xlsx"))
