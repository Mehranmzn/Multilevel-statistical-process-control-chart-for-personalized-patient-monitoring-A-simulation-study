###############################################################
# Libraries
# Note: 'dplyr' is required for data wrangling in functions below.
###############################################################
library(qcc)
library(lme4)
library(dplyr)

###############################################################
# Mean Estimation for an Individual Patient (Online, Mixed-Effects)
# Purpose:
#   Estimate per-step mean and residuals for a single patient using a mixed-
#   effects model calibrated with a sampled "bank" of patients.
# Arguments:
#   df            : Data frame of the target patient's observations (must include ID)
#   df_normal     : Population data frame for sampling calibration patients
#   variable_name : "Flow" or "Motor_power" (determines model formula)
#   patient_number: ID of the target patient
#   s1            : Calibration window length (online steps)
#   s2            : Number of rows per sampled bank patient
#   N1            : Number of patients sampled for the bank
# Returns:
#   Data frame with added columns: calibrated_resid, mean_estimation, sigma,
#   bank_p_size, bank_size, calib_size
###############################################################
mean_estimator_ind_patient <- function(df, df_normal, variable_name, patient_number, s1, s2, N1) {
  
  # -- Select calibration cohort ("bank") --------------------------------------
  selected_patients <- df_normal %>%
    dplyr::distinct(ID) %>%
    dplyr::sample_n(N1) %>%
    dplyr::pull(ID)
  
  bank_data <- df_normal %>%
    dplyr::filter(ID %in% selected_patients) %>%
    dplyr::group_by(ID) %>%
    dplyr::slice(1:s2) %>%
    dplyr::ungroup()
  
  # -- Slice target patient and initialize containers --------------------------
  df <- df %>% dplyr::filter(ID == patient_number)
  df$calibrated_resid <- 0
  df$mean_estimation  <- 0
  df$sigma            <- 0
  
  kl   <- nrow(df)
  iter <- min(s1, kl)
  online_mean <- numeric(kl)
  
  # -- Local helper to fit the mixed-effects model -----------------------------
  fit_model <- function(data, formula) {
    lmer(formula, data = data, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  }
  
  # -- Choose formula by variable_name -----------------------------------------
  if (variable_name == "Flow") {
    formula      <- Flow ~ (1 + norm_speed | ID) + HCT
    response_var <- "Flow"
  } else if (variable_name == "Motor_power") {
    formula      <- Motor_power ~ (1 + norm_speed | ID)
    response_var <- "Motor_power"
  } else {
    stop("Invalid variable_name")
  }
  
  # -- Online estimation over the first 'iter' points --------------------------
  for (j in seq_len(iter)) {
    a <- df[1:j, ]
    new_df <- dplyr::bind_rows(bank_data, a)
    new_model <- fit_model(new_df, formula)
    the_estimated_mean <- predict(new_model, a)
    
    online_mean[j] <- the_estimated_mean[j]
    df[j, "calibrated_resid"] <- df[[response_var]][j] - the_estimated_mean[j]
    df[j, "mean_estimation"]  <- the_estimated_mean[j]
  }
  
  # -- For remaining points, hold last online mean constant --------------------
  if (iter < kl) {
    df[(iter + 1):kl, "calibrated_resid"] <- df[(iter + 1):kl, response_var] - online_mean[iter]
    df[(iter + 1):kl, "mean_estimation"]  <- online_mean[iter]
  }
  
  # -- Rolling sigma over the online window, then carry forward ----------------
  df$sigma[1] <- 2
  for (s in 2:iter) df$sigma[s] <- sqrt(var(df[1:s, "calibrated_resid"]))
  if (iter < kl) df[(iter + 1):kl, "sigma"] <- df[iter, "sigma"]
  
  # -- Provenance of calibration -----------------------------------------------
  df$bank_p_size <- length(selected_patients)
  df$bank_size   <- s2
  df$calib_size  <- s1
  
  return(df)
}

###############################################################
# Residuals for All Patients
# Purpose:
#   Apply mean_estimator_ind_patient to each patient ID.
# Returns:
#   Row-bound data frame of per-patient outputs.
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
# OOC Detection (Original Variant)
# Purpose:
#   Compute EWMA residuals and dynamic LCL/UCL per patient, then flag OOC points.
# Notes:
#   Keeps original loop structure; uses df_patient$sigma consistently.
###############################################################
ooc_for_all_patients <- function(df_contain_resid, L1, L2, lambda, Variable_name, s1) {
  
  patients_ID <- unique(df_contain_resid$ID)
  df_all_patients_with_ooc <- data.frame()
  
  s1 <- s1 * 2 - 1
  
  for (j in seq_along(patients_ID)) {
    df_patient <- df_contain_resid[df_contain_resid$ID == patients_ID[j], ]
    number_of_obs <- nrow(df_patient)
    
    # EWMA residuals (uses 'qcc::ewmaSmooth')
    df_patient$ewma_resid <- ewmaSmooth(
      x = seq_len(number_of_obs),
      y = df_patient$calibrated_resid,
      lambda = lambda
    )$y
    
    df_patient$lcl <- 0
    df_patient$ucl <- 0
    
    # Iterate over observations to build windowed limits
    for (p in 1:number_of_obs) {
      
      start <- p
      end   <- p + s1
      
      iterator <- 1
      coef <- lambda / (2 - lambda)
      
      if (L1 >= L2) {
        diif <- L1 - L2
        for (k in 1:s1) {
          a <- 1 - ((1 - lambda)^(2 * k))
          df_patient[k, "lcl"] <- mean(df_patient[1:k, "calibrated_resid"][[1]]) -
            (L1 - (diif * iterator) / 60) * df_patient$sigma[k] * sqrt(coef * a)
          df_patient[k, "ucl"] <- mean(df_patient[1:k, "calibrated_resid"][[1]]) +
            (L1 - (diif * iterator) / 60) * df_patient$sigma[k] * sqrt(coef * a)
          iterator <- iterator + 1
        }
      } else {
        diif <- L2 - L1
        for (k in 1:s1) {
          a <- 1 - ((1 - lambda)^(2 * k))
          df_patient[k, "lcl"] <- mean(df_patient[1:k, "calibrated_resid"][[1]]) -
            (L1 + (diif * iterator) / s1) * df_patient$sigma[k] * sqrt(coef * a)
          df_patient[k, "ucl"] <- mean(df_patient[1:k, "calibrated_resid"][[1]]) +
            (L1 + (diif * iterator) / s1) * df_patient$sigma[k] * sqrt(coef * a)
          iterator <- iterator + 1
        }
      }
    }
    
    # First-point guard (carry from second)
    df_patient$ucl[1] <- df_patient$ucl[2]
    df_patient$lcl[1] <- df_patient$lcl[2]
    
    # Propagate missing/zero limits and compute OOC
    df_patient$ooc <- 0L
    for (p in 1:number_of_obs) {
      if (df_patient$ucl[p] == 0 || is.na(df_patient$ucl[p] == TRUE)) {
        df_patient$ucl[p] <- df_patient$ucl[p - 1]
        df_patient$lcl[p] <- df_patient$lcl[p - 1]
      }
      if (df_patient$ewma_resid[p] >= df_patient$ucl[p] ||
          df_patient$ewma_resid[p] <= df_patient$lcl[p]) {
        df_patient$ooc[p] <- 1L
      }
    }
    
    # Annotate parameters
    df_patient$L1     <- L1
    df_patient$L2     <- L2
    df_patient$lambda <- lambda
    
    df_all_patients_with_ooc <- rbind(df_all_patients_with_ooc, df_patient)
  }
  
  return(df_all_patients_with_ooc)
}

###############################################################
# OOC Detection (Optimized Variant)
# Purpose:
#   Vectorized per-patient LCL/UCL build for first 's1' points,
#   then carry forward; compute OOC flags.
# Notes:
#   Preserves behavior; uses precomputed 'coef' and absolute step.
###############################################################
ooc_optimized <- function(df_contain_resid, L1, L2, lambda, Variable_name, s1) {
  
  patients_ID <- unique(df_contain_resid$ID)
  result_list <- vector("list", length(patients_ID))
  
  coef <- lambda / (2 - lambda)
  s1   <- s1 * 2 - 1
  
  for (j in seq_along(patients_ID)) {
    df_patient <- df_contain_resid[df_contain_resid$ID == patients_ID[j], ]
    
    number_of_obs <- nrow(df_patient)
    s1_local <- min(s1, number_of_obs)
    
    # Initialize limits
    df_patient$lcl <- numeric(number_of_obs)
    df_patient$ucl <- numeric(number_of_obs)
    
    diif        <- abs(L1 - L2)
    step_factor <- diif / s1_local
    iterator    <- 1
    
    # Build limits over first 's1_local' points
    for (k in 1:s1_local) {
      a          <- 1 - ((1 - lambda)^(2 * k))
      mean_resid <- mean(df_patient[1:k, ]$ewma_resid)
      sigma_k    <- df_patient$sigma[k]
      lcl_ucl_factor <- L1 - step_factor * iterator
      
      df_patient$lcl[k] <- mean_resid - lcl_ucl_factor * sigma_k * sqrt(coef * a)
      df_patient$ucl[k] <- mean_resid + lcl_ucl_factor * sigma_k * sqrt(coef * a)
      
      iterator <- iterator + 1
    }
    
    # Carry forward the last computed limits
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
