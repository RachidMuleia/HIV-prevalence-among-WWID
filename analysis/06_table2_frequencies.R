# Frequency distribution 

source("analysis/01_prepare_data.R")

freq_one <- function(data, v) {
  x <- data[[v]]
  ok <- !is.na(x) & !(as.character(x) %in% c("", "NA"))
  xc <- as.character(x[ok])
  levs <- if (is.factor(x)) levels(x) else unique(xc)
  tab <- table(factor(xc, levels = levs), useNA = "no")
  tab <- tab[as.integer(tab) > 0]
  n_known <- sum(tab)
  data.frame(
    var = v,
    cat = names(tab),
    n = as.integer(tab),
    n_known = n_known,
    pct = round(100 * as.integer(tab) / n_known, 1),
    n_pct = paste0(as.integer(tab), " (", round(100 * as.integer(tab) / n_known, 1), ")"),
    stringsAsFactors = FALSE
  )
}

table2_freq <- do.call(rbind, lapply(covariates, function(v) freq_one(wwid, v)))
rownames(table2_freq) <- NULL

print(table2_freq[, c("var", "cat", "n_pct")])
rio::export(table2_freq, file.path(output_dir, "table2_frequencies.xlsx"))
message("WWID N = ", nrow(wwid), "; percentages exclude missing.")
