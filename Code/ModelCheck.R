# =============================================================================
# Model Reliability Diagnostics for Bat Trend Models (glmmTMB)
# TREND-FOCUSED VERSION v3
#
# Key updates from v2:
#   - Dispersion stat now assessed directionally:
#       stat << 1 => underdispersion (simulated variance > observed)
#       stat >> 1 => overdispersion  (observed variance > simulated)
#   - Underdispersion with high RE variance is treated as a DHARMa/RE
#     interaction artefact rather than a true model failure
#   - Revised disqualifier thresholds reflect this distinction
#
# Overall flag categories:
#   "Report"                   - clean model, trend estimate reliable
#   "Report with Caveats"      - minor issues, note in methods
#   "Interpret Direction Only" - fit issues; direction informative,
#                                magnitude unreliable
#   "Do Not Report"            - convergence failure, insufficient data,
#                                true severe overdispersion, or YearS
#                                completely unreliable
# =============================================================================

library(glmmTMB)
library(DHARMa)
library(dplyr)
library(tibble)

# -----------------------------------------------------------------------------
# SECTION 1: Helper functions
# -----------------------------------------------------------------------------

#' Check convergence of a glmmTMB model
check_convergence <- function(fit) {
  if (is.null(fit)) {
    return(list(pass = FALSE, message = "Model failed to fit (NULL)"))
  }
  conv_code <- tryCatch(fit$fit$convergence, error = function(e) NA)
  hess_warning <- tryCatch({
    pdHess <- fit$sdr$pdHess
    if (is.null(pdHess)) NA else pdHess
  }, error = function(e) NA)
  
  if (is.na(conv_code)) {
    return(list(pass = FALSE, message = "Could not extract convergence code"))
  }
  if (conv_code != 0) {
    return(list(pass = FALSE, message = paste0("Non-zero convergence code: ", conv_code)))
  }
  if (!is.na(hess_warning) && !hess_warning) {
    return(list(pass = FALSE, message = "Hessian not positive definite (singular/near-singular fit)"))
  }
  return(list(pass = TRUE, message = "Converged cleanly"))
}

#' Extract the YearS fixed effect estimate and SE
check_years_effect <- function(fit) {
  if (is.null(fit)) {
    return(list(estimate = NA, se = NA, ci_lo = NA, ci_hi = NA,
                se_ratio = NA, p_value = NA))
  }
  coef_table <- tryCatch(
    summary(fit)$coefficients$cond,
    error = function(e) NULL
  )
  if (is.null(coef_table) || !"YearS" %in% rownames(coef_table)) {
    return(list(estimate = NA, se = NA, ci_lo = NA, ci_hi = NA,
                se_ratio = NA, p_value = NA))
  }
  est  <- coef_table["YearS", "Estimate"]
  se   <- coef_table["YearS", "Std. Error"]
  pval <- coef_table["YearS", "Pr(>|z|)"]
  
  ci_lo    <- est - 1.96 * se
  ci_hi    <- est + 1.96 * se
  se_ratio <- if (abs(est) > 0) se / abs(est) else NA
  
  list(
    estimate = round(est, 4),
    se       = round(se, 4),
    ci_lo    = round(ci_lo, 4),
    ci_hi    = round(ci_hi, 4),
    p_value  = round(pval, 4),
    se_ratio = round(se_ratio, 2)
  )
}

#' Run DHARMa dispersion, KS uniformity, and zero-inflation tests
#'
#' Dispersion flag logic:
#'   stat << 1 (< 1):  underdispersion relative to simulation.
#'                     With high RE variance this commonly reflects the fact
#'                     that DHARMa simulates from marginal predictions while
#'                     random effects absorb substantial spatial variance.
#'                     KS test is used as the primary fit check in this case.
#'   stat >> 1 (> 1):  overdispersion. Model not capturing enough variance.
#'                     More likely a genuine fit problem worth acting on.
check_dharma <- function(fit, n_sim = 500, seed = 42) {
  if (is.null(fit)) {
    return(list(
      dispersion_stat  = NA, dispersion_p    = NA,
      dispersion_flag  = "Do Not Report", dispersion_dir  = NA,
      ks_stat          = NA, ks_p            = NA, ks_flag = "Do Not Report",
      zeroinfl_stat    = NA, zeroinfl_p      = NA, zeroinfl_flag = "Do Not Report"
    ))
  }
  
  sim_res <- tryCatch({
    set.seed(seed)
    simulateResiduals(fittedModel = fit, n = n_sim, plot = FALSE)
  }, error = function(e) {
    cat("  DHARMa simulation failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(sim_res)) {
    return(list(
      dispersion_stat = NA, dispersion_p = NA,
      dispersion_flag = "Could not simulate", dispersion_dir = NA,
      ks_stat = NA, ks_p = NA, ks_flag = "Could not simulate",
      zeroinfl_stat = NA, zeroinfl_p = NA, zeroinfl_flag = "Could not simulate"
    ))
  }
  
  # --- Dispersion test ---
  disp_test <- tryCatch(testDispersion(sim_res, plot = FALSE), error = function(e) NULL)
  disp_stat <- if (!is.null(disp_test)) round(disp_test$statistic, 4) else NA
  disp_p    <- if (!is.null(disp_test)) round(disp_test$p.value, 4)   else NA
  
  # Direction: underdispersed (stat < 1) vs overdispersed (stat > 1)
  disp_dir  <- if (is.na(disp_stat)) NA else if (disp_stat < 1) "under" else "over"
  
  disp_flag <- if (is.na(disp_p) || is.na(disp_stat)) {
    "Could not test"
    
  } else if (disp_p >= 0.05) {
    "OK"
    
    # --- Underdispersion branch ---
    # stat << 1: simulated variance >> observed variance
    # Commonly an artefact of high RE variance inflating DHARMa's simulated
    # spread. Severity binned for transparency; KS test is primary check.
  } else if (disp_dir == "under" && disp_stat < 0.1) {
    "Underdispersion — likely RE/DHARMa artefact (stat very low; check RE variance)"
  } else if (disp_dir == "under" && disp_stat < 0.5) {
    "Underdispersion — moderate; check RE variance"
  } else if (disp_dir == "under") {
    "Mild underdispersion (significant, stat 0.5-1)"
    
    # --- Overdispersion branch ---
    # stat >> 1: observed variance >> simulated variance
    # Genuine fit concern; model not capturing enough spread.
  } else if (disp_dir == "over" && disp_stat > 2.0) {
    "Severe overdispersion (stat > 2)"
  } else if (disp_dir == "over" && disp_stat > 1.5) {
    "Moderate overdispersion (stat 1.5-2)"
  } else {
    "Mild overdispersion (significant, stat 1-1.5)"
  }
  
  # --- KS uniformity test (corresponds to QQ plot) ---
  # Tests whether the overall marginal distribution matches the chosen family.
  # This is the primary distributional fit check when underdispersion is flagged.
  ks_res  <- tryCatch(testUniformity(sim_res, plot = FALSE), error = function(e) NULL)
  ks_stat <- if (!is.null(ks_res)) round(ks_res$statistic, 4) else NA
  ks_p    <- if (!is.null(ks_res)) round(ks_res$p.value, 4)   else NA
  ks_flag <- if (is.na(ks_p)) {
    "Could not test"
  } else if (ks_p < 0.001) {
    "Distributional fit poor (KS p < 0.001)"
  } else if (ks_p < 0.05) {
    "Distributional fit marginal (KS p 0.001-0.05)"
  } else {
    "OK"
  }
  
  # --- Zero-inflation test ---
  zi_test  <- tryCatch(testZeroInflation(sim_res, plot = FALSE), error = function(e) NULL)
  zi_stat  <- if (!is.null(zi_test)) round(zi_test$statistic, 3) else NA
  zi_p     <- if (!is.null(zi_test)) round(zi_test$p.value, 4)   else NA
  zi_flag  <- if (is.na(zi_p)) {
    "Could not test"
  } else if (zi_p < 0.05) {
    "Zero-inflation detected"
  } else {
    "OK"
  }
  
  list(
    dispersion_stat = disp_stat,
    dispersion_p    = disp_p,
    dispersion_flag = disp_flag,
    dispersion_dir  = disp_dir,
    ks_stat         = ks_stat,
    ks_p            = ks_p,
    ks_flag         = ks_flag,
    zeroinfl_stat   = zi_stat,
    zeroinfl_p      = zi_p,
    zeroinfl_flag   = zi_flag
  )
}

#' Count GRTS cells used in a fitted model
count_grts_cells <- function(fit) {
  if (is.null(fit)) return(NA_integer_)
  tryCatch({
    ngrps <- summary(fit)$ngrps$cond
    if ("GRTS.Cell.ID" %in% names(ngrps)) as.integer(ngrps[["GRTS.Cell.ID"]]) else NA_integer_
  }, error = function(e) NA_integer_)
}

#' Extract total number of observations
count_obs <- function(fit) {
  if (is.null(fit)) return(NA_integer_)
  tryCatch(nobs(fit), error = function(e) NA_integer_)
}

#' Extract random effect variance for GRTS.Cell.ID
check_re_variance <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  tryCatch({
    vc <- VarCorr(fit)$cond
    if ("GRTS.Cell.ID" %in% names(vc)) {
      round(as.numeric(vc[["GRTS.Cell.ID"]]), 4)
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
}

# -----------------------------------------------------------------------------
# SECTION 2: Overall reliability flag
#
# Underdispersion with high RE variance (>= 0.5) is treated as a likely
# DHARMa/RE artefact. The KS test becomes the primary distributional fit
# check in that scenario.
#
# Overdispersion (stat > 1) is treated as a genuine fit concern and
# escalates the flag accordingly.
# -----------------------------------------------------------------------------

assign_overall_flag <- function(conv_pass, n_grts, n_obs,
                                disp_flag, disp_dir, ks_flag, zi_flag,
                                se_ratio, re_var) {
  
  # === HARD DISQUALIFIERS ===
  if (!conv_pass)                                            return("Do Not Report")
  #if (!is.na(n_grts) && n_grts < 10)                        return("Do Not Report") #keep for BC full models but remove for regional models
  if (!is.na(n_obs)  && n_obs  < 150)                        return("Do Not Report")
  # Only severe *over*dispersion disqualifies — underdispersion handled below
  if (grepl("Severe overdispersion", disp_flag))             return("Do Not Report")
  # YearS estimate completely swamped by uncertainty
  if (!is.na(se_ratio) && se_ratio > 3)                      return("Do Not Report")
  
  # === INTERPRET DIRECTION ONLY ===
  direction_only <- FALSE
  
  # Genuine overdispersion at moderate level: magnitude unreliable
  if (grepl("Moderate overdispersion", disp_flag))           direction_only <- TRUE
  
  # Underdispersion + poor KS: likely true distributional misfit,
  # not just a DHARMa/RE artefact
  if (!is.na(disp_dir) && disp_dir == "under" &&
      grepl("poor", ks_flag))                                direction_only <- TRUE
  
  # Poor KS on its own (distributional family a poor fit)
  if (grepl("Distributional fit poor", ks_flag))             direction_only <- TRUE
  
  # Zero-inflation biases counts toward zero, affecting magnitude
  if (grepl("Zero-inflation", zi_flag))                      direction_only <- TRUE
  
  # Wide but non-extreme SE: direction visible, CI too wide for magnitude
  if (!is.na(se_ratio) && se_ratio > 1 && se_ratio <= 3)    direction_only <- TRUE
  
  if (direction_only)                                        return("Interpret Direction Only")
  
  # === REPORT WITH CAVEATS ===
  caution <- FALSE
  if (!is.na(n_grts) && n_grts < 25)                        caution <- TRUE
  if (!is.na(n_obs)  && n_obs  < 100)                       caution <- TRUE
  # Mild overdispersion: note in methods
  if (grepl("Mild overdispersion", disp_flag))               caution <- TRUE
  # Underdispersion with clean KS: likely artefact, worth noting
  if (!is.na(disp_dir) && disp_dir == "under" &&
      grepl("OK", ks_flag))                                  caution <- TRUE
  # Marginal KS: distributional fit slightly off but not disqualifying
  if (grepl("marginal", ks_flag))                            caution <- TRUE
  if (!is.na(re_var) && re_var < 0.001)                      caution <- TRUE
  
  if (caution)                                               return("Report with Caveats")
  
  return("Report")
}

# -----------------------------------------------------------------------------
# SECTION 3: Main diagnostic loop
# -----------------------------------------------------------------------------

run_model_diagnostics <- function(fit_list, n_sim = 500) {
  
  cat("Running model diagnostics for", length(fit_list), "species...\n")
  cat("DHARMa simulations:", n_sim, "per model\n")
  cat("Mode: TREND-FOCUSED v3 (underdispersion/overdispersion treated separately)\n\n")
  
  results <- lapply(names(fit_list), function(sp) {
    
    cat("---", sp, "---\n")
    fit <- fit_list[[sp]]
    
    # 1. Convergence
    conv   <- check_convergence(fit)
    cat("  Convergence:", conv$message, "\n")
    
    # 2. Sample sizes
    n_grts <- count_grts_cells(fit)
    n_obs  <- count_obs(fit)
    re_var <- check_re_variance(fit)
    cat("  GRTS cells:", n_grts, "| Observations:", n_obs,
        "| RE variance:", re_var, "\n")
    
    # 3. YearS effect
    ye <- check_years_effect(fit)
    cat("  YearS estimate:", ye$estimate, "| SE:", ye$se,
        "| SE ratio:", ye$se_ratio, "| p:", ye$p_value, "\n")
    
    # 4. DHARMa
    if (conv$pass) {
      dh <- check_dharma(fit, n_sim = n_sim)
    } else {
      dh <- list(
        dispersion_stat = NA, dispersion_p = NA,
        dispersion_flag = "Skipped (convergence failure)",
        dispersion_dir  = NA,
        ks_stat = NA, ks_p = NA,
        ks_flag = "Skipped (convergence failure)",
        zeroinfl_stat = NA, zeroinfl_p = NA,
        zeroinfl_flag = "Skipped (convergence failure)"
      )
    }
    cat("  Dispersion:", dh$dispersion_flag,
        "(stat =", dh$dispersion_stat, ", p =", dh$dispersion_p, ")\n")
    cat("  KS uniformity:", dh$ks_flag,
        "(stat =", dh$ks_stat, ", p =", dh$ks_p, ")\n")
    cat("  Zero-inflation:", dh$zeroinfl_flag,
        "(stat =", dh$zeroinfl_stat, ", p =", dh$zeroinfl_p, ")\n")
    
    # 5. Overall flag
    overall <- assign_overall_flag(
      conv_pass = conv$pass,
      n_grts    = n_grts,
      n_obs     = n_obs,
      disp_flag = dh$dispersion_flag,
      disp_dir  = dh$dispersion_dir,
      ks_flag   = dh$ks_flag,
      zi_flag   = dh$zeroinfl_flag,
      se_ratio  = ye$se_ratio,
      re_var    = re_var
    )
    cat("  >> Overall:", overall, "\n\n")
    
    tibble(
      Species              = sp,
      Convergence          = conv$message,
      N_GRTS_cells         = n_grts,
      N_observations       = n_obs,
      RE_variance_GRTS     = re_var,
      YearS_estimate       = ye$estimate,
      YearS_SE             = ye$se,
      YearS_CI_lo          = ye$ci_lo,
      YearS_CI_hi          = ye$ci_hi,
      YearS_p              = ye$p_value,
      SE_ratio             = ye$se_ratio,
      Dispersion_stat      = dh$dispersion_stat,
      Dispersion_p         = dh$dispersion_p,
      Dispersion_direction = dh$dispersion_dir,
      Dispersion_flag      = dh$dispersion_flag,
      KS_stat              = dh$ks_stat,
      KS_p                 = dh$ks_p,
      KS_flag              = dh$ks_flag,
      ZeroInflation_stat   = dh$zeroinfl_stat,
      ZeroInflation_p      = dh$zeroinfl_p,
      ZeroInflation_flag   = dh$zeroinfl_flag,
      Overall_flag         = overall
    )
  })
  
  bind_rows(results)
}

# -----------------------------------------------------------------------------
# SECTION 4: Run diagnostics and save outputs
# -----------------------------------------------------------------------------

diag_table <- run_model_diagnostics(best_fits, n_sim = 500)

# Print summary
cat("\n=== SUMMARY TABLE ===\n")
print(
  diag_table %>%
    select(Species, N_GRTS_cells, N_observations,
           YearS_estimate, YearS_SE, YearS_p, SE_ratio,
           Dispersion_flag, Dispersion_direction,
           KS_flag, ZeroInflation_flag,
           Overall_flag) %>%
    as.data.frame(),
  row.names = FALSE
)

# Flag counts
cat("\n=== FLAG COUNTS ===\n")
print(table(diag_table$Overall_flag))

# Legend
cat("\n=== FLAG LEGEND ===\n")
cat("Report                   : All diagnostics pass.\n")
cat("Report with Caveats      : Minor issues (sample size, mild overdispersion,\n")
cat("                           or underdispersion with clean KS). Note in methods.\n")
cat("Interpret Direction Only : Moderate fit problems or wide SE ratio.\n")
cat("                           Report trend direction and significance only;\n")
cat("                           do not emphasise exact magnitude of estimate.\n")
cat("Do Not Report            : Convergence failure, severe overdispersion,\n")
cat("                           insufficient data, or YearS completely unreliable.\n")
cat("\nNote on underdispersion: stat << 1 with high RE variance is commonly a\n")
cat("DHARMa/RE interaction artefact. KS uniformity test is used as the primary\n")
cat("distributional fit check in this scenario.\n")

# Save
write.csv(diag_table, "model_diagnostics_full.csv", row.names = FALSE)
cat("\nFull diagnostics saved to: model_diagnostics_full.csv\n")

# -----------------------------------------------------------------------------
# SECTION 5: DHARMa plots for flagged models
# -----------------------------------------------------------------------------

flagged_species <- diag_table %>%
  filter(Overall_flag != "Report") %>%
  pull(Species)

if (length(flagged_species) > 0) {
  cat("\nGenerating DHARMa residual plots for flagged species:",
      paste(flagged_species, collapse = ", "), "\n")
  
  pdf("model_diagnostics_flagged_plots.pdf", width = 10, height = 6)
  for (sp in flagged_species) {
    fit <- bc.fits.cov.singlet[[sp]]
    if (!is.null(fit)) {
      tryCatch({
        set.seed(42)
        sim_res <- simulateResiduals(fit, n = 500, plot = FALSE)
        par(mfrow = c(1, 2))
        plot(sim_res)
        flag <- diag_table %>% filter(Species == sp) %>% pull(Overall_flag)
        mtext(paste0(sp, "  [", flag, "]"),
              side = 3, line = -1.5, outer = TRUE, cex = 1.1, font = 2)
      }, error = function(e) {
        cat("  Could not plot", sp, ":", conditionMessage(e), "\n")
      })
    }
  }
  dev.off()
  cat("Plots saved to: model_diagnostics_flagged_plots.pdf\n")
}  
  
  #--------------------------------------------CHECK REGIONAL MODELS-----------------------------------------------------------#
# -----------------------------------------------------------------------------
# SECTION 3: Main diagnostic loop (updated for nested region → species list)
# -----------------------------------------------------------------------------

run_model_diagnostics <- function(fit_list_nested, n_sim = 500) {
  
  cat("Running model diagnostics across", length(fit_list_nested), "regions...\n")
  cat("DHARMa simulations:", n_sim, "per model\n")
  cat("Mode: TREND-FOCUSED v3 (underdispersion/overdispersion treated separately)\n\n")
  
  results <- lapply(names(fit_list_nested), function(region) {
    
    cat("\n========== REGION:", region, "==========\n")
    fit_list <- fit_list_nested[[region]]
    
    region_results <- lapply(names(fit_list), function(sp) {
      
      cat("---", sp, "---\n")
      fit <- fit_list[[sp]]
      
      # 1. Convergence
      conv   <- check_convergence(fit)
      cat("  Convergence:", conv$message, "\n")
      
      # 2. Sample sizes
      n_grts <- count_grts_cells(fit)
      n_obs  <- count_obs(fit)
      re_var <- check_re_variance(fit)
      cat("  GRTS cells:", n_grts, "| Observations:", n_obs,
          "| RE variance:", re_var, "\n")
      
      # 3. YearS effect
      ye <- check_years_effect(fit)
      cat("  YearS estimate:", ye$estimate, "| SE:", ye$se,
          "| SE ratio:", ye$se_ratio, "| p:", ye$p_value, "\n")
      
      # 4. DHARMa
      if (conv$pass) {
        dh <- check_dharma(fit, n_sim = n_sim)
      } else {
        dh <- list(
          dispersion_stat = NA, dispersion_p = NA,
          dispersion_flag = "Skipped (convergence failure)",
          dispersion_dir  = NA,
          ks_stat = NA, ks_p = NA,
          ks_flag = "Skipped (convergence failure)",
          zeroinfl_stat = NA, zeroinfl_p = NA,
          zeroinfl_flag = "Skipped (convergence failure)"
        )
      }
      cat("  Dispersion:", dh$dispersion_flag,
          "(stat =", dh$dispersion_stat, ", p =", dh$dispersion_p, ")\n")
      cat("  KS uniformity:", dh$ks_flag,
          "(stat =", dh$ks_stat, ", p =", dh$ks_p, ")\n")
      cat("  Zero-inflation:", dh$zeroinfl_flag,
          "(stat =", dh$zeroinfl_stat, ", p =", dh$zeroinfl_p, ")\n")
      
      # 5. Overall flag
      overall <- assign_overall_flag(
        conv_pass = conv$pass,
        n_grts    = n_grts,
        n_obs     = n_obs,
        disp_flag = dh$dispersion_flag,
        disp_dir  = dh$dispersion_dir,
        ks_flag   = dh$ks_flag,
        zi_flag   = dh$zeroinfl_flag,
        se_ratio  = ye$se_ratio,
        re_var    = re_var
      )
      cat("  >> Overall:", overall, "\n\n")
      
      tibble(
        Region               = region,
        Species              = sp,
        Convergence          = conv$message,
        N_GRTS_cells         = n_grts,
        N_observations       = n_obs,
        RE_variance_GRTS     = re_var,
        YearS_estimate       = ye$estimate,
        YearS_SE             = ye$se,
        YearS_CI_lo          = ye$ci_lo,
        YearS_CI_hi          = ye$ci_hi,
        YearS_p              = ye$p_value,
        SE_ratio             = ye$se_ratio,
        Dispersion_stat      = dh$dispersion_stat,
        Dispersion_p         = dh$dispersion_p,
        Dispersion_direction = dh$dispersion_dir,
        Dispersion_flag      = dh$dispersion_flag,
        KS_stat              = dh$ks_stat,
        KS_p                 = dh$ks_p,
        KS_flag              = dh$ks_flag,
        ZeroInflation_stat   = dh$zeroinfl_stat,
        ZeroInflation_p      = dh$zeroinfl_p,
        ZeroInflation_flag   = dh$zeroinfl_flag,
        Overall_flag         = overall
      )
    })
    
    bind_rows(region_results)
  })
  
  bind_rows(results)
}

# -----------------------------------------------------------------------------
# SECTION 4: Run diagnostics and save outputs
# -----------------------------------------------------------------------------

diag_table <- run_model_diagnostics(all.region.fits, n_sim = 500)

# Print summary
cat("\n=== SUMMARY TABLE ===\n")
print(
  diag_table %>%
    select(Region, Species, N_GRTS_cells, N_observations,
           YearS_estimate, YearS_SE, YearS_p, SE_ratio,
           Dispersion_flag, Dispersion_direction,
           KS_flag, ZeroInflation_flag,
           Overall_flag) %>%
    as.data.frame(),
  row.names = FALSE
)

# Flag counts (overall and by region)
cat("\n=== FLAG COUNTS (OVERALL) ===\n")
print(table(diag_table$Overall_flag))

cat("\n=== FLAG COUNTS BY REGION ===\n")
print(table(diag_table$Region, diag_table$Overall_flag))

# Legend
cat("\n=== FLAG LEGEND ===\n")
cat("Report                   : All diagnostics pass.\n")
cat("Report with Caveats      : Minor issues (sample size, mild overdispersion,\n")
cat("                           or underdispersion with clean KS). Note in methods.\n")
cat("Interpret Direction Only : Moderate fit problems or wide SE ratio.\n")
cat("                           Report trend direction and significance only;\n")
cat("                           do not emphasise exact magnitude of estimate.\n")
cat("Do Not Report            : Convergence failure, severe overdispersion,\n")
cat("                           insufficient data, or YearS completely unreliable.\n")

# Save
write.csv(diag_table, "Data/Analyzed/model_diagnostics_regional_full.csv", row.names = FALSE)
cat("\nFull diagnostics saved to: model_diagnostics_regional_full.csv\n")

# -----------------------------------------------------------------------------
# SECTION 5: DHARMa plots for flagged models
# -----------------------------------------------------------------------------

flagged_models <- diag_table %>%
  filter(Overall_flag != "Report") %>%
  select(Region, Species, Overall_flag)

if (nrow(flagged_models) > 0) {
  cat("\nGenerating DHARMa residual plots for", nrow(flagged_models), "flagged models...\n")
  
  pdf("model_diagnostics_regional_flagged_plots.pdf", width = 10, height = 6)
  for (i in seq_len(nrow(flagged_models))) {
    region <- flagged_models$Region[i]
    sp     <- flagged_models$Species[i]
    flag   <- flagged_models$Overall_flag[i]
    
    fit <- all.region.fits[[region]][[sp]]
    if (!is.null(fit)) {
      tryCatch({
        set.seed(42)
        sim_res <- simulateResiduals(fit, n = 500, plot = FALSE)
        par(mfrow = c(1, 2))
        plot(sim_res)
        mtext(paste0(sp, " | ", region, "  [", flag, "]"),
              side = 3, line = -1.5, outer = TRUE, cex = 1.1, font = 2)
      }, error = function(e) {
        cat("  Could not plot", sp, "in", region, ":", conditionMessage(e), "\n")
      })
    }
  }
  dev.off()
  cat("Plots saved to: model_diagnostics_regional_flagged_plots.pdf\n")
}  

#Get the cleaned version (removing any model results that are sus) of the model outputs for report
do_not_report <- diag_table %>%
  filter(Overall_flag == "Do Not Report") %>%
  select(Region, Species)

# Filter the summary table
all.region.slope.summary.singlet.filtered <- all.region.slope.summary.singlet %>%
  anti_join(do_not_report, by = c("Region" = "Region", "SpeciesGroup" = "Species"))

write.csv(all.region.slope.summary.singlet.filtered,"Data/Analyzed/region.slope.summary.singlet.csv")
