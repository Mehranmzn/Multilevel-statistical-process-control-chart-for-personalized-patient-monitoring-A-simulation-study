###############################################################
# Libraries
# Note: dplyr is required for distinct(), sample_n(), pull(), bind_rows(), etc.
###############################################################
library(qcc)
library(lme4)
library(dplyr)

###############################################################
# mean_estimator_ind_patient
# Purpose:
#   Online mean estimation for a single patient using a mixed-effects model
#   calibrated with a sampled "bank" of patients; computes calibrated residuals
#   and an evolving sigma.
# Arguments:
#   df            : Data frame containing the target patient's rows (incl. ID)
#   df_normal     : Population data for sampling calibration patients
#   variable_name : "Flow" or "Motor_power"
#   patient_number: Patient ID to process
#   s1            : Online calibration window length
#   s2            : Number of rows to take per sampled bank patient
#   N1            : Number of bank patients to sample
# Returns:
#   Data frame augmented with calibrated_resid, mean_estimation, sigma and
#   calibration meta-data (bank_p_size, bank_size, calib_size).
###############################################################
mean_estimator_ind_patient <- function(df, df_normal, variable_name, patient_number, s1, s2, N1) {
  
  # -- Select a calibration "bank" of patients --------------------------------
  selected_patients <- df_normal %>%
    distinct(ID) %>%
    sample_n(N1) %>%
    pull(ID)
  
  bank_data <- df_normal %>%
    filter(ID %in% selected_patients) %>%
    group_by(ID) %>%
    slice(1:s2) %>%
    ungroup()
  
  # -- Filter target patient rows and initialize fields -----------------------
  df <- df %>% filter(ID == patient_number)
  df$calibrated_resid <- 0
  df$mean_estimation  <- 0
  df$sigma            <- 0
  
  kl   <- nrow(df)
  iter <- min(s1, kl)
  online_mean <- numeric(kl)
  
  # -- Helper: fit mixed-effects model (REML = FALSE, bobyqa) -----------------
  fit_model <- function(data, formula) {
    lmer(
      formula,
      data    = data,
      REML    = FALSE,
      control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
    )
  }
  
  # -- Model formula by variable_name (as provided) ---------------------------
  if (variable_name == "Flow") {
    formula      <- Flow ~ (1 | ID) + HCT
    response_var <- "Flow"
  } else if (variable_name == "Motor_power") {
    formula      <- Motor_power ~ (1 | ID)
    response_var <- "Motor_power"
  } else {
    stop("Invalid variable_name")
  }
  
  # -- Online estimation over first 'iter' points -----------------------------
  for (j in seq_len(iter)) {
    a <- df[1:j, ]
    new_df <- bind_rows(bank_data, a)
    new_model <- fit_model(new_df, formula)
    the_estimated_mean <- predict(new_model, a)
    
    online_mean[j]            <- the_estimated_mean[j]
    df[j, "calibrated_resid"] <- df[[response_var]][j] - the_estimated_mean[j]
    df[j, "mean_estimation"]  <- the_estimated_mean[j]
  }
  
  # -- Hold last online mean for remaining points -----------------------------
  if (iter < kl) {
    df[(iter + 1):kl, "calibrated_resid"] <- df[(iter + 1):kl, response_var] - online_mean[iter]
    df[(iter + 1):kl, "mean_estimation"]  <- online_mean[iter]
  }
  
  # -- Sigma: grow over online window, then carry forward ---------------------
  df$sigma[1] <- 2
  for (s in 2:iter) {
    df$sigma[s] <- sqrt(var(df[1:s, "calibrated_resid"]))
  }
  if (iter < kl) {
    df[(iter + 1):kl, "sigma"] <- df[iter, "sigma"]
  }
  
  # -- Calibration meta-data ---------------------------------------------------
  df$bank_p_size <- length(selected_patients)
  df$bank_size   <- s2
  df$calib_size  <- s1
  
  return(df)
}

###############################################################
# patient_specific_residual
# Purpose:
#   Apply mean_estimator_ind_patient to all patient IDs in df and row-bind.
# Arguments:
#   df, df_normal, variable_name, s1, s2, N1 : see above.
# Returns:
#   Data frame with all patients' residualized outputs.
###############################################################
patient_specific_residual <- function(df, df_normal, variable_name, s1, s2, N1) {
  u <- sort(unique(df$ID))
  
  residuals_list <- lapply(u, function(patient_id) {
    mean_estimator_ind_patient(
      df             = df[df$ID == patient_id, ],
      df_normal      = df_normal,
      patient_number = patient_id,
      variable_name  = variable_name,
      s1 = s1, s2 = s2, N1 = N1
    )
  })
  
  df_contain_resid_all_patient <- do.call(rbind, residuals_list)
  return(df_contain_resid_all_patient)
}

###############################################################
# ooc_optimized
# Purpose:
#   Build dynamic LCL/UCL envelopes over a window (length = s1*2 - 1),
#   then carry forward; compute OOC flags via EWMA residuals.
# Arguments:
#   df_contain_resid : Data with columns ID, ewma_resid, sigma, etc.
#   L1, L2           : Control limit parameters
#   lambda           : EWMA smoothing parameter
#   Variable_name    : Unused here (kept for signature compatibility)
#   s1               : Base window length (transformed to 2*s1 - 1)
# Returns:
#   Row-bound data frame with lcl/ucl/ooc and annotated L1/L2/lambda.
###############################################################
ooc_optimized <- function(df_contain_resid, L1, L2, lambda, Variable_name, s1) {
  
  # -- Setup -------------------------------------------------------------------
  patients_ID <- unique(df_contain_resid$ID)
  result_list <- vector("list", length(patients_ID))
  coef        <- lambda / (2 - lambda)
  s1_window   <- s1 * 2 - 1
  
  # -- Per patient -------------------------------------------------------------
  for (j in seq_along(patients_ID)) {
    
    df_patient <- df_contain_resid[df_contain_resid$ID == patients_ID[j], ]
    number_of_obs <- nrow(df_patient)
    
    # Local window length per patient
    s1_local <- min(s1_window, number_of_obs)
    
    # Initialize limits
    df_patient$lcl <- numeric(number_of_obs)
    df_patient$ucl <- numeric(number_of_obs)
    
    # Step size between L1 and L2 across the window
    diif        <- abs(L1 - L2)
    step_factor <- diif / s1_local
    iterator    <- 1
    
    # Build LCL/UCL for first s1_local points using EWMA mean and sigma_k
    for (k in 1:s1_local) {
      a          <- 1 - ((1 - lambda)^(2 * k))
      mean_resid <- mean(df_patient[1:k, ]$ewma_resid)
      sigma_k    <- df_patient$sigma[k]
      lcl_ucl_factor <- L1 - step_factor * iterator
      
      df_patient$lcl[k] <- mean_resid - lcl_ucl_factor * sigma_k * sqrt(coef * a)
      df_patient$ucl[k] <- mean_resid + lcl_ucl_factor * sigma_k * sqrt(coef * a)
      
      iterator <- iterator + 1
    }
    
    # Carry forward remaining limits (if any)
    if (s1_local < number_of_obs) {
      df_patient$lcl[s1_local:number_of_obs] <- df_patient$lcl[s1_local]
      df_patient$ucl[s1_local:number_of_obs] <- df_patient$ucl[s1_local]
    }
    
    # OOC flag
    df_patient$ooc <- as.integer(
      df_patient$ewma_resid >= df_patient$ucl |
        df_patient$ewma_resid <= df_patient$lcl
    )
    
    # Annotate parameters
    df_patient$L1     <- L1
    df_patient$L2     <- L2
    df_patient$lambda <- lambda
    
    result_list[[j]] <- df_patient
  }
  
  df_all_patients_with_ooc <- do.call(rbind, result_list)
  return(df_all_patients_with_ooc)
}
