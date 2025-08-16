###############################################################
# shift_imposer.R
# Purpose:
#   Impose a specified shift pattern on a time series from a given
#   position onward per patient, then update calibrated residuals
#   and EWMA residuals.
# Arguments:
#   df_bf_shift : Data frame containing time-series with columns:
#                 ID, calibrated_resid, mean_estimation, and the target var
#   position    : Integer index at which the shift starts (1-based)
#   shift_mag   : Numeric magnitude of the shift (sign allowed)
#   var         : Column name (string) to be shifted (e.g., "Flow")
#   type        : One of c("mean_change", "sudden_jump", "random",
#                          "slow_trend_nonl", "slow_trend_l")
#   lambda      : EWMA smoothing parameter in (0, 1)
# Returns:
#   Data frame with updated columns for each patient:
#   - var (shifted as specified)
#   - calibrated_resid
#   - ewma_resid
#   - type_shift
#   - shift_mag
###############################################################

library(parallel)  # retained for compatibility
library(dplyr)     # bind_rows, dplyr verbs
library(qcc)       # ewmaSmooth

shift_imposer <- function(df_bf_shift,
                          position,
                          shift_mag,
                          var,
                          type,
                          lambda) {
  
  # -- Unique patient IDs ------------------------------------------------------
  u <- unique(df_bf_shift$ID)
  
  # -- Collector for processed patients ---------------------------------------
  DF <- data.frame()
  
  # -- Process each patient independently -------------------------------------
  for (k in u) {
    
    df <- df_bf_shift[df_bf_shift$ID == k, ]
    obs_n <- nrow(df)
    
    # -- Compute baseline EWMA of residuals (pre-shift) -----------------------
    df[, "ewma_resid"] <- ewmaSmooth(
      x = seq(1, obs_n, 1),
      y = df$calibrated_resid,
      lambda = lambda
    )$y
    
    # -- Require at least 360 observations (original logic preserved) ---------
    if (obs_n < 360) next
    
    # -- Shift patterns --------------------------------------------------------
    if (type == "mean_change") {
      # Apply constant shift from 'position' onward
      df[position:obs_n, var] <- df[position:obs_n, var] + shift_mag
      
      # Recompute residuals and EWMA from 'position' onward
      df[position:obs_n, "calibrated_resid"] <-
        df[position:obs_n, var] - df[position:obs_n, "mean_estimation"]
      
      df[position:obs_n, "ewma_resid"] <- ewmaSmooth(
        seq(1, obs_n - position + 1, 1),
        df[position:obs_n, ]$calibrated_resid,
        lambda = lambda,
        start  = mean(df[1:position, ]$calibrated_resid)
      )$y
      
    } else if (type == "sudden_jump") {
      # Apply a single jump at 'position'
      df[position, var] <- df[position, var] + shift_mag
      
      # Recompute residual at the jump and EWMA onward
      df[position, "calibrated_resid"] <- df[position, var] - df[position, "mean_estimation"]
      
      df[position:obs_n, "ewma_resid"] <- ewmaSmooth(
        seq(1, obs_n - position + 1, 1),
        df[position:obs_n, ]$calibrated_resid,
        lambda = lambda,
        start  = mean(df[1:position, ]$calibrated_resid)
      )$y
      
    } else if (type == "random") {
      # Apply random Gaussian shifts from 'position' onward (sd = |shift_mag|)
      random_shifts <- rnorm(obs_n - position + 1, mean = 0, sd = abs(shift_mag))
      df[position:obs_n, var] <- df[position:obs_n, var] + random_shifts
      
      # Recompute residuals and EWMA onward
      df[position:obs_n, "calibrated_resid"] <-
        df[position:obs_n, var] - df[position:obs_n, "mean_estimation"]
      
      df[position:obs_n, "ewma_resid"] <- ewmaSmooth(
        seq(1, obs_n - position + 1, 1),
        df[position:obs_n, ]$calibrated_resid,
        lambda = lambda,
        start  = mean(df[1:position, ]$calibrated_resid)
      )$y
      
    } else if (type == "slow_trend_nonl") {
      # Apply a nonlinear (exponential) slow trend from 'position' onward
      trend_shift <- shift_mag * (1 - exp(-(seq(0, obs_n - position) / (obs_n - position) * 5)))
      df[position:obs_n, var] <- df[position:obs_n, var] + trend_shift
      
      # Recompute residuals and EWMA onward
      df[position:obs_n, "calibrated_resid"] <-
        df[position:obs_n, var] - df[position:obs_n, "mean_estimation"]
      
      df[position:obs_n, "ewma_resid"] <- ewmaSmooth(
        seq(1, obs_n - position + 1, 1),
        df[position:obs_n, ]$calibrated_resid,
        lambda = lambda,
        start  = mean(df[1:position, ]$calibrated_resid)
      )$y
      
    } else if (type == "slow_trend_l") {
      # Apply a linear slow trend from 'position' onward
      trend_shift <- seq(0, shift_mag, length.out = (obs_n - position + 1))
      df[position:obs_n, var] <- df[position:obs_n, var] + trend_shift
      
      # Recompute residuals and EWMA onward
      df[position:obs_n, "calibrated_resid"] <-
        df[position:obs_n, var] - df[position:obs_n, "mean_estimation"]
      
      df[position:obs_n, "ewma_resid"] <- ewmaSmooth(
        seq(1, obs_n - position + 1, 1),
        df[position:obs_n, ]$calibrated_resid,
        lambda = lambda,
        start  = mean(df[1:position, ]$calibrated_resid)
      )$y
      
    } else {
      stop("Invalid type specified")
    }
    
    # -- Annotate shift meta-data ---------------------------------------------
    df[, "type_shift"] <- type
    df[, "shift_mag"]  <- shift_mag
    
    # -- Accumulate results ----------------------------------------------------
    DF <- dplyr::bind_rows(DF, df)
  }
  
  return(DF)
}
