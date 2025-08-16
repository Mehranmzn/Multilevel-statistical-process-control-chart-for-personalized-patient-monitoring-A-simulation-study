###############################################################
# evaluator
# Purpose:
#   Summarize OOC performance and related statistics per patient.
# Arguments:
#   df_given : Data frame with per-observation fields, incl.:
#              ID, ooc, calibrated_resid, mean_estimation, ewma_resid,
#              lcl, ucl, L1, L2, lambda, bank_p_size, bank_size, calib_size,
#              and the measurement column named by 'var'.
#   s1       : Calibration window length.
#   position : Index where evaluation window begins.
#   var      : Name of measured variable (e.g., "Flow", "Motor_power").
# Returns:
#   A data frame of per-patient metrics and day-bucket OOC counts.
###############################################################
evaluator <- function(df_given, s1, position, var) {
  
  num_patients <- unique(df_given$ID)
  
  # Convenience indices around 'position'
  posit2  <- position - 1
  posit15 <- position - 2
  posit16 <- position - 4
  
  # Evaluation end index (kept as in original code)
  end_data <- 56
  
  # Output container (column names preserved)
  results <- data.frame(
    ID = character(),
    FalsePositives = integer(),
    TruePositives = integer(),
    NumOnesInOOC = integer(),
    EarliestOOC = integer(),
    AvgCalibratedResidual = numeric(),
    MeanEstimated = numeric(),
    AvgEwmasCritical = numeric(),
    AvgEwmasNormal = numeric(),
    BankPSize = numeric(),
    BankSize = numeric(),
    CalibSize = numeric(),
    EwmaResid = numeric(),
    Lcl = numeric(),
    Ucl = numeric(),
    Ooc = integer(),
    L1 = numeric(),
    L2 = numeric(),
    Lambda = numeric(),
    stringsAsFactors = FALSE,
    true_mean_calib = numeric(),
    true_mean_full = numeric(),
    true_mean_after_cal = numeric(),
    estimated_mean_full = numeric(),
    estimated_mean_after_cal = numeric(),
    estimated_mean_cal = numeric()
  )
  
  # Per-patient aggregation
  for (k in num_patients) {
    df_temp <- df_given[df_given$ID == k, ]
    
    # Core metrics (guard NAs as in original)
    false_positives        <- sum(df_temp$ooc[1:posit2] == 1, na.rm = TRUE)
    true_positives         <- ifelse(any(df_temp$ooc[position:end_data] == 1, na.rm = TRUE), 1, 0)
    num_ones_in_ooc        <- sum(df_temp$ooc[position:end_data] == 1, na.rm = TRUE)
    avg_calibrated_resid   <- mean(df_temp$calibrated_resid[1:s1], na.rm = TRUE)
    mean_estimated         <- df_temp$mean_estimation[posit2]
    avg_ewmas_critical     <- mean(df_temp$ewma_resid[position:end_data], na.rm = TRUE)
    avg_ewmas_normal       <- mean(df_temp$ewma_resid[posit2], na.rm = TRUE)  # single index by design
    
    # Specific day flags before position (names preserved)
    sum_day_15             <- sum(df_temp$ooc[posit15] == 1, na.rm = TRUE)
    sum_day_16             <- sum(df_temp$ooc[posit16] == 1, na.rm = TRUE)
    
    # True vs. estimated means over windows
    true_mean_calib        <- mean(df_temp[[var]][1:s1], na.rm = TRUE)
    true_mean_full         <- mean(df_temp[[var]], na.rm = TRUE)
    true_mean_after_cal    <- mean(df_temp[[var]][(s1 + 1):posit2], na.rm = TRUE)
    estimated_mean_full    <- mean(df_temp$mean_estimation[1:posit2], na.rm = TRUE)
    estimated_mean_after_cal <- mean(df_temp$mean_estimation[(s1 + 1):posit2], na.rm = TRUE)
    estimated_mean_cal     <- mean(df_temp$mean_estimation[1:s1], na.rm = TRUE)
    
    # Day-bucket OOC counts from 'end_data' down to 'position' in steps of 2
    day_values <- list()
    for (i in seq(end_data, position, by = -2)) {
      day_name <- paste0("d_", (end_data - i) / 2 + 1)  # dynamic day column
      day_values[[day_name]] <- sum(df_temp$ooc[(i - 2):i] >= 1, na.rm = TRUE)
    }
    day_df <- as.data.frame(t(unlist(day_values)))
    
    # Earliest OOC index in the evaluation window
    ooc_window   <- df_temp$ooc[position:end_data]
    earliest_ooc <- ifelse(any(ooc_window == 1, na.rm = TRUE),
                           position + which(ooc_window == 1)[1] - 1,
                           NA)
    
    # Carry meta-parameters (first value per patient)
    bank_p_size <- df_temp$bank_p_size[1]
    bank_size   <- df_temp$bank_size[1]
    calib_size  <- df_temp$calib_size[1]
    ewma_resid  <- mean(df_temp$ewma_resid, na.rm = TRUE)
    lcl         <- mean(df_temp$lcl, na.rm = TRUE)
    ucl         <- mean(df_temp$ucl, na.rm = TRUE)
    ooc_total   <- sum(df_temp$ooc, na.rm = TRUE)
    l1          <- df_temp$L1[1]
    l2          <- df_temp$L2[1]
    lambda      <- df_temp$lambda[1]
    
    # Bind a single-row result with dynamic day columns
    results <- rbind(
      results,
      data.frame(
        ID = k,
        FalsePositives = false_positives,
        TruePositives = true_positives,
        NumOnesInOOC = num_ones_in_ooc,
        EarliestOOC = earliest_ooc,
        AvgCalibratedResidual = avg_calibrated_resid,
        MeanEstimated = mean_estimated,
        AvgEwmasCritical = avg_ewmas_critical,
        AvgEwmasNormal = avg_ewmas_normal,
        BankPSize = bank_p_size,
        BankSize = bank_size,
        CalibSize = calib_size,
        EwmaResid = ewma_resid,
        Lcl = lcl,
        Ucl = ucl,
        Ooc = ooc_total,
        L1 = l1,
        L2 = l2,
        Lambda = lambda,
        location_15 = sum_day_15,   # kept original name
        locaiton_16 = sum_day_16,   # typo preserved intentionally
        day_df,
        true_mean_calib = true_mean_calib,
        true_mean_full = true_mean_full,
        true_mean_after_cal = true_mean_after_cal,
        estimated_mean_full = estimated_mean_full,
        estimated_mean_after_cal = estimated_mean_after_cal,
        estimated_mean_cal = estimated_mean_cal,
        check.names = FALSE
      )
    )
  }
  
  return(results)
}
