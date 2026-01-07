extract_lab <- function(data, variable) {
  if (variable %in% names(data)) {
    # 1. Convert to character
    all_values <- as.character(data[[variable]])
    
    # 2. Filter out NA values (represented as "NA" character string)
    #    and also filter out literal NA values if the conversion didn't handle them.
    categories <- unique(all_values[!is.na(all_values) & all_values != "NA"])
    
    # 3. Handle cases where the variable might be factor/labeled and contain a level called 'NA' 
    #    that is not a missing value, but for demographic analysis, it's safer to exclude.
    #    However, given the RDS context, we assume "NA" means missing data.
    
    return(categories)
  } else {
    stop(paste("Variable", variable, "not found in data."))
  }
}

rds_prev <- function(data, x_var, var_outcome, population_size, dec, dec.i, ...){
  
  categories <- extract_lab(data = data, variable = x_var)
  
  # --- Safe Extraction Helper ---
  safe_extract <- function(result, type = c("estimate", "se")) {
    is_failed <- (length(result) == 1 && is.na(result)) | is.null(result)
    
    if (is_failed) {
      return(NA) # Returns a numeric NA
    }
    
    type <- match.arg(type)
    
    if (type == "estimate") {
      return(as.numeric(result$estimate[1]))
    } else if (type == "se") {
      return(as.numeric(attr(result, 'bsresult')$se_estimate[1]))
    }
  }
  
  # 1. Calculate prevalence for each category using withCallingHandlers
  prevalence <- lapply(categories, function(x){
    
    # Use withCallingHandlers to selectively ignore non-fatal warnings
    prevalence_result <- withCallingHandlers(
      {
        # Attempt the RDS calculation
        rds_res <- RDS.bootstrap.intervals(data,
                                           outcome.variable = var_outcome,
                                           weight.type = "Gile's SS",
                                           uncertainty = 'Gile', confidence.level = 0.95,
                                           number.of.bootstrap.samples = 1000,
                                           subset = data[, x_var] == x,
                                           N = population_size,...)
        
        # Check if the result is invalid after the call
        if(is.null(rds_res) || (length(rds_res) == 1 && is.na(rds_res))) {
          return(NA) 
        }
        return(rds_res)
      },
      # --- Warning Handler ---
      warning = function(w) {
        # Fatal Warning: Outcome has only one level (prevalence is 0% or 100%)
        if(grepl("only one level", w$message)) {
          message(paste("Skipping", x, "due to one-level outcome:", w$message))
          invokeRestart("muffleWarning") # Suppress printing
          return(NA) # Force return NA
        }
        
        # Non-Fatal Warnings: Missing data handled by RDS package
        if (grepl("trait variable values were missing|network sizes were missing|presume these are", w$message)) {
          # Accept warning, allow execution to continue and return the result
          message(paste("RDS Warning accepted for", x, ":", w$message))
          invokeRestart("muffleWarning") # Suppress printing
          return(NULL) # Return NULL to indicate the handler didn't override the result
        } else {
          # Other critical warnings (treated as failure)
          message(paste("Critical RDS Warning:", w$message))
          invokeRestart("muffleWarning")
          return(NA)
        }
      },
      # --- Error Handler ---
      error = function(e) {
        message(paste("Error for", x, ":", e$message))
        return(NA)
      }
    )
    # If the handler returned NA, use NA. Otherwise, use the calculated result.
    return(prevalence_result)
  })
  
  # 2. Extract estimates and SEs using the safe_extract helper
  prev_raw <- sapply(prevalence, safe_extract, type = "estimate")
  se_raw <- sapply(prevalence, safe_extract, type = "se")
  
  # --- MODIFICATION: Create separate columns for Prevalence and SE ---
  prev_estimate <- round(prev_raw * 100, dec)
  se_estimate_pc <- round(se_raw * 100, dec.i) # SE expressed as percentage
  
  # 3. Create final data frame with new structure (no CI string, no p_value)
  prev_df <- data.frame(
    var = rep(x_var, length(categories)),
    cat = categories,
    prevalence = prev_estimate, # New column name
    se = se_estimate_pc# New column name
  )
  
  return(prev_df)
}

# rds_prev_all remains the same to stack all variables' categories for ONE province
rds_prev_all <- function(data, vec_var, var_outcome, ...){
  prev_all <- lapply(vec_var, \(x){
    rds_prev(data = data, x_var = x, var_outcome = var_outcome,...)
  })
  all_prev <- do.call(rbind, prev_all)
  return(all_prev)
}







library(dplyr)
library(purrr) # Needed for the reduce function

rds_prov_all <- function(data_list, var_vector, var_outcome, pop_size, dec, dec.i, site){
  
  # ... (Error check remains the same)
  
  # 1. Calculate results for ALL provinces (same as original)
  all_results_list <- lapply(seq_along(data_list), function(i){
    # ... (Inside calculation remains the same, calculating df and renaming columns)
    
    # Determine the site identifier (same as original)
    data <- data_list[[i]]
    current_pop_size <- pop_size[i]
    
    if (site %in% names(data)) {
      current_site_name <- unique(as.character(data[[site]]))[1]
    } else {
      current_site_name <- paste("Site", i)
    }
    
    # Calculate RDS prevalence for the current province/site
    df <- rds_prev_all(
      data = data, 
      vec_var = var_vector, 
      var_outcome = var_outcome,
      population_size = current_pop_size, 
      dec = dec, 
      dec.i = dec.i
    )
    
    # Rename the results columns (same as original)
    names(df)[names(df) %in% c("prevalence", "se")] <- 
      paste(c("prevalence", "se"), current_site_name, sep = "_")
    
    return(df)
  })
  
  # --- 2. MODIFICATION: Use full_join for safe merging ---
  # Use reduce and full_join to combine all results by matching 'var' and 'cat'
  # This automatically handles differing row counts by inserting NA where a category is missing.
  combined_df <- all_results_list %>%
    purrr::reduce(dplyr::full_join, by = c("var", "cat"))
  
  combined_df <- combined_df %>%
    
    # --- Sorting Logic ---
    # 3a. Extract the numeric prefix for reliable sorting (e.g., '1' from '1_16-17')
    # This uses a regular expression to grab the digit right after the colon in 'cat'
    mutate(SORT_ORDER = as.numeric(gsub("^.*:(\\d+)_.*", "\\1", cat))) %>%
    
    # 3b. Sort the dataframe based on the extracted number
    # If SORT_ORDER is NA (for non-numbered categories), the original order is maintained.
    arrange(var, SORT_ORDER) %>% 
    
    # 3c. Remove the temporary sort column
    select(!all_of("SORT_ORDER")) %>% 
    
    # 3d. Relocate 'var' and 'cat' to the front
    relocate(var, .before = 1) %>% 
    relocate(cat, .after = var)
  
  return(combined_df)
}

rds_prov_all(data_list = data_list_rds, var_vector = var_wid[20:27], var_outcome = 'RESUL_HIV_PREV',
             pop_size = pop_size,dec=1, dec.i=1,site = 'SITE_CITY')


var_wid
library(metafor)
library(dplyr)

library(metafor)
library(dplyr)

meta_analyze_rds <- function(rds_prov_all_output) {
  
  # 1. Identify Prevalence and SE Columns (unchanged)
  prev_cols <- names(rds_prov_all_output)[grepl("^prevalence_", names(rds_prov_all_output))]
  se_cols <- names(rds_prov_all_output)[grepl("^se_", names(rds_prov_all_output))]
  
  if (length(prev_cols) == 0 || length(prev_cols) != length(se_cols)) {
    stop("Input data must contain matching columns named 'prevalence_SiteX' and 'se_SiteX'.")
  }
  
  meta_results <- rds_prov_all_output %>%
    rowwise() %>% 
    mutate(
      study_data_raw = list(
        data.frame(
          yi_prop = c_across(all_of(prev_cols)) / 100,
          vi_prop = (c_across(all_of(se_cols)) / 100)^2
        )
      ),
      
      # Filter out NA rows (sites that failed)
      study_data_clean = list(
        study_data_raw %>% 
          filter(!is.na(yi_prop) & !is.na(vi_prop))
      ),
      
      n_valid_studies = nrow(study_data_clean),
      
      # Run the Fixed-Effects Meta-Analysis 
      meta_model = list(
        if (n_valid_studies >= 2) {
          tryCatch({
            rma(yi = yi_prop, vi = vi_prop, data = study_data_clean, method = "FE")
          }, error = function(e) {
            warning(paste("Meta-analysis failed for", SITE_CATEGORY, ":", e$message))
            return(NULL) # <--- **FIXED** to return NULL from the tryCatch handler
          })
        } else {
          NULL # <--- **FIXED**: Using NULL instead of return(NULL)
        }
      )
    ) %>%
    
    # 3. Extract and Format Results (unchanged logic)
    mutate(
      model_exists = !is.null(meta_model),
      
      pooled_prevalence = ifelse(model_exists, round(as.numeric(meta_model$b) * 100, 2), NA),
      ci_lower = ifelse(model_exists, round(as.numeric(meta_model$ci.lb) * 100, 2), NA),
      ci_upper = ifelse(model_exists, round(as.numeric(meta_model$ci.ub) * 100, 2), NA),
      q_pvalue = ifelse(model_exists && !is.null(meta_model$Qp), round(as.numeric(meta_model$Qp), 3), NA),
      
      pooled_prevalence_ci = paste0(
        pooled_prevalence, " (", ci_lower, " - ", ci_upper, ")"
      )
    ) %>%
    
    # 4. Select Final Columns (unchanged)
    select(
      SITE_CATEGORY, 
      pooled_prevalence_ci, 
      pooled_prevalence, 
      ci_lower, 
      ci_upper, 
      q_pvalue,
      n_valid_studies
    ) %>%
    as.data.frame()
  
  return(meta_results)
}


descriptive_data <- rds_prov_all(data_list = data_list_rds, var_vector = var_wid[-19], var_outcome = 'RESUL_HIV_PREV',
             pop_size = pop_size,dec=1, dec.i=1,site = 'SITE_CITY')



RDS.bootstrap.intervals(data_list_rds[[1]],
                        outcome.variable = "RESUL_HIV_PREV",
                        weight.type = "Gile's SS",
                        uncertainty = 'Gile', confidence.level =0.95,
                        number.of.bootstrap.samples=1000,
                        subset =data_list_rds[[1]][, var_wid[6]] == "2+",
                        N = pop_size[1])


descriptive_data <- descriptive_data|>
  mutate(SITE_CATEGORY = paste(var, cat, sep = ":"))|>
  relocate(SITE_CATEGORY, before = var)


meta_analyze_rds(descriptive_data)

export(meta_analyze_rds(descriptive_data), file = "weighted_estimates_HIV_prevalence.xlsx")


descriptive_data <- rds_prov_all(data_list = data_list_rds, var_vector = var_wid[-19], var_outcome = 'RESUL_HIV_PREV',
                                 pop_size = pop_size,dec=1, dec.i=1,site = 'SITE_CITY')


descriptive_data_HCV <- rds_prov_all(data_list = data_list_rds, var_vector = var_wid[-19], var_outcome = 'RESUL_HCV_CAT',
                                 pop_size = pop_size,dec=1, dec.i=1,site = 'SITE_CITY')

descriptive_data_HCV <- descriptive_data_HCV|>
  mutate(SITE_CATEGORY = paste(var, cat, sep = ":"))|>
  relocate(SITE_CATEGORY, before = var)




meta_analyze_rds(descriptive_data_HCV)

export(meta_analyze_rds(descriptive_data_HCV), file = "weighted_estimates_HICV_prevalence.xlsx")


