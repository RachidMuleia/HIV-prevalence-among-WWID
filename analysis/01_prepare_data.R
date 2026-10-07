# Load analytic file, recode factors used in Tables 1–2, build city RDS objects.
setwd("/Users/rachidmuleia/Dropbox/INS/PID/PAPER WWID_AURIA/HIV-prevalence-among-WWID")

pkgs <- c("RDS", "metafor", "survey", "sandwich", "lmtest", "geepack", "dplyr", "rio")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Please install: ", p)
  library(p, character.only = TRUE)
}))

source("R/00_config.R")
source("R/01_functions.R")

if (!file.exists(data_file)) {
  stop(
    "Analytic file not found: ", data_file,
    "\nSet WWID_DATA_FILE to the local DATA_PID_CLUSTER.csv path."
  )
}

pid <- read.csv(data_file, header = TRUE, na.strings = "")
pid$ELINJMU[is.na(pid$ELINJMU) | pid$ELINJMU < 1] <- 1

pid$response <- ifelse(
  as.character(pid$RESUL_HIV_PREV) == "1_POSITIVO", 1,
  ifelse(as.character(pid$RESUL_HIV_PREV) == "2_NEGATIVO", 0, NA)
)

pid$EDU_CAT1 <- dplyr::recode(
  as.character(pid$EDU_CAT1),
  "1_SEM_ESCOLA" = "2_PRIMARIO",
  "2_PRIMARIO" = "2_PRIMARIO",
  "3_SECUNDARIO_SUP" = "3_SECUNDARIO_SUP"
)
if ("EDU_CAT1" %in% names(pid)) {
  pid$EDU_CAT1 <- factor(pid$EDU_CAT1)
  pid$EDU_CAT1 <- stats::relevel(pid$EDU_CAT1, ref = "2_PRIMARIO")
}

pid$DEACT <- factor(pid$DEACT, levels = c("Nao", "Sim"))


pid$DRUG6_CAT <- dplyr::case_when(
  as.character(pid$DRUG6_CAT) == "1_HEROINA" ~ "1_HEROINA",
  as.character(pid$DRUG6_CAT) %in% c("2_COCAINA", "3_OUTRA") ~ "2_COCAINA_OUTRA",
  TRUE ~ NA_character_
)

# Table 2 and GEE: Daily combines 1–3 and ≥5 times/day; reference is No_inject.
pid$FREQ_INJ <- dplyr::case_when(
  as.character(pid$ID6_FRQ_CAT) == "2_4_vezes_por_mes" ~ "Month",
  as.character(pid$ID6_FRQ_CAT) == "3_2_7_vezes_por_semana" ~ "Weekly",
  as.character(pid$ID6_FRQ_CAT) %in% c("4_1_3_vezes_por_dia", "5+_vezes_por_dia") ~ "Daily",
  as.character(pid$ID6_FRQ_CAT) == "1_NAO_INJECTOU" ~ "No_inject",
  TRUE ~ NA_character_
)
pid$FREQ_INJ <- factor(pid$FREQ_INJ, levels = c("No_inject", "Month", "Weekly", "Daily"))

# AUDIT-C: women, hazardous if total >= 3.
pid <- pid |>
  mutate(
    ALCOHOL1_PT = case_when(
      ALFRQ == "Nenhuma" ~ 0,
      ALFRQ == "Mensal_ou_menos" ~ 1,
      ALFRQ == "2_4_vezes_por_mes" ~ 2,
      ALFRQ %in% c("2_3_vezes_por_semana", "Semanal") ~ 3,
      ALFRQ %in% c("4_ou_mais_vezes_por_semana", "Diariamente") ~ 4,
      .default = NA
    ),
    ALCOHOL2_PT = case_when(
      ALDAY == "1_ou_2" ~ 0,
      ALDAY == "3_ou_4" ~ 1,
      ALDAY == "5_ou_6" ~ 2,
      ALDAY == "6_a_9" ~ 3,
      ALDAY == "10+" ~ 4,
      .default = NA
    ),
    ALCOHOL3_PT = case_when(
      ALBNGE == "Nunca" ~ 0,
      ALBNGE == "Menos_de_uma_vez_por_mes" ~ 1,
      ALBNGE == "Pelo_menos_uma_vez_por_mes" ~ 2,
      ALBNGE == "Pelo_menos_uma_vez_po_semana" ~ 3,
      ALBNGE == "Diari_ou_quase_todos_os_dias" ~ 4,
      .default = NA
    )
  ) |>
  rowwise() |>
  mutate(ALCOHOL_PTSOMA = sum(c_across(ALCOHOL1_PT:ALCOHOL3_PT), na.rm = TRUE)) |>
  mutate(
    AUDIT_SCORE = case_when(
      ALCOHOL_PTSOMA < 3 ~ "2_NAO_ABUSIVO",
      ALCOHOL_PTSOMA >= 3 ~ "1_ABUSIVO",
      ALFRQ1 == "Nao" ~ "2_NAO_ABUSIVO",
      .default = NA
    )
  ) |>
  ungroup()

covariates <- c(
  "AGE_CAT", "SITE_CITY", "EDU_CAT1", "DEMARSTA_CAT", "DEACT",
  "IDADE_SEX_CAT", "MSEX6_CAT", "LASTREL6_CAT", "ORIENTA_SEX1",
  "MSEX7_CAT", "AUDIT_SCORE", "PARTILHA_SERINGA", "IDNSTE_CAT1",
  "IDADE_DROGA_INJ", "DRUG6_CAT", "ANOS_INJECT_CAT", "FREQ_INJ",
  "IDNSTE1_CAT", "ODNAX_CAT", "TRATRDN_CAT", "PEEREDU1_CAT",
  "GRAVIDEZ", "SINTOMAS_ITS", "STGPFSX_CAT", "STGXCLP_CAT"
)

for (v in intersect(covariates, names(pid))) {
  pid[[v]] <- factor(pid[[v]])
}

# Table 2 multivariable specification and reference levels.
pid$DEACT <- factor(pid$DEACT, levels = c("Nao", "Sim"))
pid$ANOS_INJECT_CAT <- factor(
  pid$ANOS_INJECT_CAT,
  levels = c("<1_ANO", "2_1-3_ANOS", "3_4-6_ANOS", "4_7+")
)
pid$FREQ_INJ <- factor(pid$FREQ_INJ, levels = c("No_inject", "Month", "Weekly", "Daily"))
pid$DRUG6_BIN <- dplyr::case_when(
  as.character(pid$DRUG6_CAT) == "1_HEROINA" ~ "Heroin",
  as.character(pid$DRUG6_CAT) == "2_COCAINA_OUTRA" ~ "Cocaine_other",
  TRUE ~ NA_character_
)
pid$DRUG6_BIN <- factor(pid$DRUG6_BIN, levels = c("Heroin", "Cocaine_other"))
pid$MSEX6_CAT <- factor(pid$MSEX6_CAT, levels = c("1_SIM", "2_NAO"))
pid$AGE_CAT <- factor(pid$AGE_CAT, levels = c("1_16-17", "2_18-24", "3_25-34", "4_35+"))
pid$SITE_CITY <- factor(pid$SITE_CITY, levels = sort(unique(as.character(pid$SITE_CITY))))
pid$DEMARSTA_CAT <- factor(
  pid$DEMARSTA_CAT,
  levels = c("1_SOLTEIRO", "2_CASADOS", "3_SEPARADO_VIUVO")
)
pid$PEEREDU1_CAT <- factor(pid$PEEREDU1_CAT, levels = c("1_SIM", "2_NAO"))
pid$GRAVIDEZ <- factor(pid$GRAVIDEZ, levels = c("1_SIM", "2_NAO"))
pid$SINTOMAS_ITS <- factor(pid$SINTOMAS_ITS, levels = c("1_SIM", "2_NAO"))
pid$STGXCLP_CAT <- factor(pid$STGXCLP_CAT, levels = c("1_SIM", "2_NAO"))



sites <- sort(unique(as.character(pid$SITE_CITY)))
if (length(sites) != length(pop_size)) {
  stop("pop_size must have one entry per city in sort(SITE_CITY).")
}

data_list_rds <- rdsat_format(
  data = pid,
  site_list = sites,
  network_size = "ELINJMU",
  recruiter_ide = "recruiter.id",
  var_site = "SITE_CITY",
  max.coupon = 5,
  id = "CODPARTICIPANTE_dem"
)

wwid <- pid[is_woman(pid), ]

cov_levels <- lapply(stats::setNames(covariates, covariates), function(v) {
  x <- pid[[v]]
  if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])))
})

message("Loaded ", nrow(pid), " PWID; ", nrow(wwid), " WWID; ", length(sites), " cities.")

