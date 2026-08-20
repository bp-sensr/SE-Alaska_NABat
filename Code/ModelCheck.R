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
  } else if (ks_p < 0.01) {
    "Distributional fit poor (KS p < 0.01)"
  } else if (ks_p < 0.05) {
    "Distributional fit marginal (KS p 0.01-0.05)"
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
#AK.fits.cov.singlet <- readRDS("Data/Analyzed/AK_fits_cov_singlet.rds") #leaving this here. Should work if you have saved the output so you dont need to re-run the model selction script
#best_fits <- lapply(AK.fits.cov.singlet, function(x) x$best_fit)

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
write.csv(diag_table, "Data/Analyzed/model_diagnostics_full.csv", row.names = FALSE)
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
  
  pdf("Data/Analyzed/model_diagnostics_flagged_plots.pdf", width = 10, height = 6)
  for (sp in flagged_species) {
    fit <- best_fits[[sp]]
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
  
  #--------------------------------------------CHECK Transect Models-----------------------------------------------------------#
# =============================================================================
# Signal & robustness diagnostics — TRANSECT abundance models
# For sparse transect abundance data the useful question about model fit is different:
# "is the YearS trend real and identifiable, or is it fitting noise / one
# lucky cell / an unidentified NB?" Those are robustness questions, and they
# can be checked honestly at small n because they do not rely on asymptotic
# distributional theory. This script replaces the DHARMa battery with:
#
#   1. Parametric bootstrap LRT for YearS  — the honest replacement for the
#      Wald p-value. Simulates from the null (no-YearS) fit and builds the
#      LRT null distribution empirically. This is the primary "is there
#      signal" check.
#   2. Leave-one-cell-out (LOCO)           — refit dropping each GRTS cell in
#      turn. If one cell flips the sign or halves the slope, the "trend" is
#      that cell, not a regional pattern.
#   3. Profile-likelihood CI for YearS     — sparse-count likelihoods are
#      asymmetric, so the Wald SE (and the SE_ratio guard) can mislead. The
#      profile CI validates or undermines that guard.
#   4. Poisson vs NB comparison            — with mostly 0/1/2 counts the NB
#      dispersion is often unidentified; a stably-fit Poisson then beats a
#      wobbly NB. Compared by AIC.
#
#
# SELF-CONTAINED: the per-species data (cell.summary) is not retained by the
# selection loop, so this script rebuilds a clean modelling frame from each
# fit's own $frame. No need to re-run or edit the selection script.
# =============================================================================

library(glmmTMB)
library(dplyr)
library(tibble)

# -----------------------------------------------------------------------------
# SECTION 1: Basic extractors (carried over unchanged from the DHARMa script)
# -----------------------------------------------------------------------------

check_convergence <- function(fit) {
  if (is.null(fit)) {
    return(list(pass = FALSE, message = "Model failed to fit (NULL)"))
  }
  conv_code <- tryCatch(fit$fit$convergence, error = function(e) NA)
  hess_ok <- tryCatch({
    pdHess <- fit$sdr$pdHess
    if (is.null(pdHess)) NA else pdHess
  }, error = function(e) NA)
  
  if (is.na(conv_code)) {
    return(list(pass = FALSE, message = "Could not extract convergence code"))
  }
  if (conv_code != 0) {
    return(list(pass = FALSE, message = paste0("Non-zero convergence code: ", conv_code)))
  }
  if (!is.na(hess_ok) && !hess_ok) {
    return(list(pass = FALSE, message = "Hessian not positive definite (singular/near-singular fit)"))
  }
  list(pass = TRUE, message = "Converged cleanly")
}

check_years_effect <- function(fit) {
  if (is.null(fit)) {
    return(list(estimate = NA, se = NA, p_value = NA, se_ratio = NA))
  }
  coef_table <- tryCatch(summary(fit)$coefficients$cond, error = function(e) NULL)
  if (is.null(coef_table) || !"YearS" %in% rownames(coef_table)) {
    return(list(estimate = NA, se = NA, p_value = NA, se_ratio = NA))
  }
  est  <- coef_table["YearS", "Estimate"]
  se   <- coef_table["YearS", "Std. Error"]
  pval <- coef_table["YearS", "Pr(>|z|)"]
  se_ratio <- if (abs(est) > 0) se / abs(est) else NA
  list(
    estimate = est,   # unrounded for downstream math; rounded only in output
    se       = se,
    p_value  = round(pval, 4),
    se_ratio = round(se_ratio, 2)
  )
}

count_grts_cells <- function(fit) {
  if (is.null(fit)) return(NA_integer_)
  tryCatch({
    ngrps <- summary(fit)$ngrps$cond
    if ("GRTS.Cell.ID" %in% names(ngrps)) as.integer(ngrps[["GRTS.Cell.ID"]]) else NA_integer_
  }, error = function(e) NA_integer_)
}

count_obs <- function(fit) {
  if (is.null(fit)) return(NA_integer_)
  tryCatch(nobs(fit), error = function(e) NA_integer_)
}

check_re_variance <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  tryCatch({
    vc <- VarCorr(fit)$cond
    if ("GRTS.Cell.ID" %in% names(vc)) round(as.numeric(vc[["GRTS.Cell.ID"]]), 4) else NA_real_
  }, error = function(e) NA_real_)
}

# -----------------------------------------------------------------------------
# SECTION 2: Rebuild a clean, refittable frame from a fit object
#
# The fit references its data by name, which is gone after the ddply loop.
# glmmTMB stores the model frame in fit$frame, but transformed terms like
# log(transect.length) live there as a column literally named
# "log(transect.length)", which cannot be re-evaluated by the formula.
# This helper sanitises names and rebuilds a matching formula + family so we
# can refit full, null, subset, and Poisson variants freely.
# -----------------------------------------------------------------------------

rebuild_from_frame <- function(fit) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required (nobars/findbars). Please install it.")
  }
  fr        <- fit$frame
  full_form <- formula(fit)
  fixed_form <- lme4::nobars(full_form)
  bars       <- lme4::findbars(full_form)
  
  resp         <- all.vars(fixed_form)[1]
  fixed_labels <- attr(terms(fixed_form), "term.labels")
  
  clean_name <- function(x) {
    x2 <- gsub("^log\\((.*)\\)$", "log_\\1", x)   # log(a.b) -> log_a.b
    x2 <- gsub("[^A-Za-z0-9._]", "_", x2)
    x2
  }
  
  names(fr) <- vapply(names(fr), clean_name, character(1))
  
  resp_clean  <- clean_name(resp)
  fixed_clean <- vapply(fixed_labels, clean_name, character(1))
  bar_terms   <- vapply(bars, function(b) paste0("(", deparse(b), ")"), character(1))
  
  rhs      <- paste(c(fixed_clean, bar_terms), collapse = " + ")
  new_form <- stats::as.formula(paste(resp_clean, "~", rhs))
  
  fam_name <- fit$modelInfo$family$family
  fam_link <- fit$modelInfo$family$link
  fam <- switch(fam_name,
                nbinom1 = glmmTMB::nbinom1(link = fam_link),
                nbinom2 = glmmTMB::nbinom2(link = fam_link),
                poisson = stats::poisson(link = fam_link),
                genpois = glmmTMB::genpois(link = fam_link),
                stop(paste("Unhandled family:", fam_name))
  )
  
  list(data = fr, formula = new_form, family = fam,
       family_name = fam_name, family_link = fam_link,
       response = resp_clean)
}

# -----------------------------------------------------------------------------
# SECTION 3: The four signal / robustness checks
# -----------------------------------------------------------------------------

# 3.1 Parametric bootstrap LRT for YearS -------------------------------------
# Re-fits full and null from scratch on each simulated response. (An earlier
# version used glmmTMB::refit(newresp=), which does not work on glmmTMB
# objects the way lme4::refit does — every replicate errored, giving 0 valid
# reps. Re-fitting directly is slower but robust.)
bootstrap_lrt <- function(rb, full2, null2, B = 499, seed = 42) {
  obs_lrt <- tryCatch(as.numeric(2 * (logLik(full2) - logLik(null2))),
                      error = function(e) NA_real_)
  
  sims <- tryCatch(as.data.frame(simulate(null2, nsim = B, seed = seed)),
                   error = function(e) NULL)
  if (is.null(sims)) return(list(obs = round(obs_lrt, 3), p = NA, nvalid = 0))
  
  full_form <- formula(full2)
  null_form <- formula(null2)
  resp      <- rb$response
  dat       <- rb$data
  
  sim_lrt <- vapply(seq_len(ncol(sims)), function(b) {
    dat[[resp]] <- sims[[b]]
    fb <- tryCatch(suppressWarnings(
      glmmTMB::glmmTMB(full_form, family = rb$family, data = dat)),
      error = function(e) NULL)
    if (is.null(fb)) return(NA_real_)
    nb <- tryCatch(suppressWarnings(
      glmmTMB::glmmTMB(null_form, family = rb$family, data = dat)),
      error = function(e) NULL)
    if (is.null(nb)) return(NA_real_)
    ll_f <- tryCatch(as.numeric(logLik(fb)), error = function(e) NA_real_)
    ll_n <- tryCatch(as.numeric(logLik(nb)), error = function(e) NA_real_)
    if (!is.finite(ll_f) || !is.finite(ll_n)) return(NA_real_)
    2 * (ll_f - ll_n)
  }, numeric(1))
  
  valid  <- sim_lrt[is.finite(sim_lrt)]
  nvalid <- length(valid)
  p <- if (nvalid == 0) NA else (1 + sum(valid >= obs_lrt)) / (1 + nvalid)
  list(obs = round(obs_lrt, 3), p = round(p, 4), nvalid = nvalid)
}

# 3.2 Leave-one-cell-out ------------------------------------------------------
loco_years <- function(rb, full_est) {
  cells <- levels(droplevels(as.factor(rb$data[["GRTS.Cell.ID"]])))
  if (length(cells) < 2) {
    return(list(min = NA, max = NA, sign_stable = NA, max_rel = NA))
  }
  ests <- vapply(cells, function(cl) {
    d <- rb$data[as.character(rb$data[["GRTS.Cell.ID"]]) != cl, , drop = FALSE]
    d[["GRTS.Cell.ID"]] <- droplevels(as.factor(d[["GRTS.Cell.ID"]]))
    f <- tryCatch(suppressWarnings(
      glmmTMB::glmmTMB(rb$formula, family = rb$family, data = d)),
      error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    ct <- tryCatch(summary(f)$coefficients$cond, error = function(e) NULL)
    if (is.null(ct) || !"YearS" %in% rownames(ct)) return(NA_real_)
    ct["YearS", "Estimate"]
  }, numeric(1))
  
  valid <- ests[is.finite(ests)]
  if (length(valid) == 0) {
    return(list(min = NA, max = NA, sign_stable = NA, max_rel = NA))
  }
  sign_stable <- all(sign(valid) == sign(full_est))
  max_rel     <- if (abs(full_est) > 0) max(abs(valid - full_est) / abs(full_est)) else NA
  list(min = round(min(valid), 4), max = round(max(valid), 4),
       sign_stable = sign_stable, max_rel = round(max_rel, 3))
}

# 3.3 Profile-likelihood CI for YearS ----------------------------------------
profile_ci <- function(full2, full_est) {
  try_ci <- function(method) {
    tryCatch(suppressWarnings(
      confint(full2, parm = "YearS", method = method, level = 0.95)),
      error = function(e) NULL)
  }
  ci <- try_ci("profile"); method <- "profile"
  if (is.null(ci)) { ci <- try_ci("uniroot"); method <- "uniroot" }
  if (is.null(ci)) { ci <- try_ci("wald");    method <- "wald"    }
  if (is.null(ci)) {
    return(list(lo = NA, hi = NA, method = "none",
                excludes_zero = NA, asym = NA))
  }
  lo <- ci[1, 1]; hi <- ci[1, 2]
  excludes_zero <- (lo > 0 & hi > 0) | (lo < 0 & hi < 0)
  left  <- full_est - lo
  right <- hi - full_est
  asym  <- if ((hi - lo) > 0) abs(right - left) / (hi - lo) else NA
  list(lo = round(lo, 4), hi = round(hi, 4), method = method,
       excludes_zero = excludes_zero, asym = round(asym, 3))
}

# 3.4 Poisson vs NB -----------------------------------------------------------
poisson_compare <- function(rb, full2) {
  pois <- tryCatch(suppressWarnings(
    glmmTMB::glmmTMB(rb$formula,
                     family = stats::poisson(link = rb$family_link),
                     data   = rb$data)),
    error = function(e) NULL)
  aic_nb <- tryCatch(AIC(full2), error = function(e) NA)
  aic_po <- if (is.null(pois)) NA else tryCatch(AIC(pois), error = function(e) NA)
  nb_disp <- tryCatch(as.numeric(sigma(full2)), error = function(e) NA)
  
  # NB collapses toward Poisson when the extra dispersion parameter buys
  # nothing: Poisson AIC within ~2 of NB means NB is effectively unidentified.
  preferred <- if (is.na(aic_po) || is.na(aic_nb)) FALSE else (aic_po <= aic_nb + 2)
  
  list(aic_nb = round(aic_nb, 2), aic_pois = round(aic_po, 2),
       nb_dispersion = round(nb_disp, 3), poisson_preferred = preferred)
}

# -----------------------------------------------------------------------------
# SECTION 4: Overall flag (signal-driven)
#
# "No Detectable Trend" is a first-class outcome here — it is the honest
# answer when the bootstrap says the data are not distinguishable from no
# trend, and it is more truthful than forcing a null result into a caveat.
# -----------------------------------------------------------------------------

assign_transect_flag <- function(conv_pass, n_grts, se_ratio,
                                 boot_p, loco_sign_stable, loco_max_rel,
                                 prof_excludes_zero, prof_asym,
                                 poisson_preferred) {
  
  # === HARD DISQUALIFIERS ===
  if (!conv_pass)                        return("Do Not Report")
  if (is.na(se_ratio) || se_ratio > 3)   return("Do Not Report")   # YearS swamped
  # Trend that flips sign when a single cell is removed is a one-cell artefact
  if (!is.na(loco_sign_stable) && !loco_sign_stable) return("Do Not Report")
  
  # === NO SIGNAL (honest null) ===
  # Bootstrap could not reject the no-trend null: the data are not saying
  # anything detectable about a YearS effect.
  if (!is.na(boot_p) && boot_p >= 0.05)  return("No Detectable Trend")
  
  # From here the bootstrap is significant (or NA and treated cautiously).
  
  # === INTERPRET DIRECTION ONLY ===
  direction_only <- FALSE
  if (is.na(boot_p))                                     direction_only <- TRUE  # couldn't bootstrap
  if (isTRUE(poisson_preferred))                        direction_only <- TRUE  # NB unidentified
  if (!is.na(se_ratio) && se_ratio > 1)                 direction_only <- TRUE  # wide CI
  if (!is.na(prof_asym) && prof_asym > 0.3)             direction_only <- TRUE  # skewed likelihood
  if (!is.na(prof_excludes_zero) && !prof_excludes_zero) direction_only <- TRUE # profile CI spans 0
  if (!is.na(loco_max_rel) && loco_max_rel > 0.5)       direction_only <- TRUE  # slope halves/doubles
  if (direction_only)                                   return("Interpret Direction Only")
  
  # === REPORT WITH CAVEATS ===
  caution <- FALSE
  if (!is.na(n_grts) && n_grts < 5)                     caution <- TRUE
  if (!is.na(loco_max_rel) && loco_max_rel > 0.25)      caution <- TRUE
  if (caution)                                          return("Report with Caveats")
  
  return("Report")
}

# -----------------------------------------------------------------------------
# SECTION 5: Per-species driver
# -----------------------------------------------------------------------------

check_transect_signal <- function(fit, sp, B = 499,
                                  do_loco = TRUE, do_profile = TRUE,
                                  do_poisson = TRUE, seed = 42) {
  
  cat("---", sp, "---\n")
  
  conv <- check_convergence(fit)
  cat("  Convergence:", conv$message, "\n")
  
  n_grts <- count_grts_cells(fit)
  n_obs  <- count_obs(fit)
  re_var <- check_re_variance(fit)
  ye     <- check_years_effect(fit)
  cat("  GRTS cells:", n_grts, "| Obs:", n_obs, "| RE var:", re_var, "\n")
  cat("  YearS est:", round(ye$estimate, 4), "| SE:", round(ye$se, 4),
      "| SE ratio:", ye$se_ratio, "| Wald p:", ye$p_value, "\n")
  
  # Defaults for when the fit is unusable
  boot <- list(obs = NA, p = NA, nvalid = 0)
  loco <- list(min = NA, max = NA, sign_stable = NA, max_rel = NA)
  prof <- list(lo = NA, hi = NA, method = "none", excludes_zero = NA, asym = NA)
  pois <- list(aic_nb = NA, aic_pois = NA, nb_dispersion = NA, poisson_preferred = NA)
  
  if (conv$pass) {
    rb    <- rebuild_from_frame(fit)
    full2 <- glmmTMB::glmmTMB(rb$formula, family = rb$family, data = rb$data)
    null2 <- glmmTMB::glmmTMB(update(rb$formula, . ~ . - YearS),
                              family = rb$family, data = rb$data)
    
    cat("  Bootstrap LRT (B =", B, ")...\n")
    boot <- bootstrap_lrt(rb, full2, null2, B = B, seed = seed)
    cat("    obs LRT:", boot$obs, "| boot p:", boot$p,
        "| valid reps:", boot$nvalid, "\n")
    
    if (do_loco) {
      cat("  Leave-one-cell-out...\n")
      loco <- loco_years(rb, full_est = ye$estimate)
      cat("    YearS range: [", loco$min, ",", loco$max, "] | sign stable:",
          loco$sign_stable, "| max rel change:", loco$max_rel, "\n")
    }
    if (do_profile) {
      prof <- profile_ci(full2, full_est = ye$estimate)
      cat("  Profile CI (", prof$method, "): [", prof$lo, ",", prof$hi,
          "] | excludes 0:", prof$excludes_zero, "| asymmetry:", prof$asym, "\n")
    }
    if (do_poisson) {
      pois <- poisson_compare(rb, full2)
      cat("  AIC  NB:", pois$aic_nb, "| Poisson:", pois$aic_pois,
          "| NB dispersion:", pois$nb_dispersion,
          "| Poisson preferred:", pois$poisson_preferred, "\n")
    }
  }
  
  overall <- assign_transect_flag(
    conv_pass          = conv$pass,
    n_grts             = n_grts,
    se_ratio           = ye$se_ratio,
    boot_p             = boot$p,
    loco_sign_stable   = loco$sign_stable,
    loco_max_rel       = loco$max_rel,
    prof_excludes_zero = prof$excludes_zero,
    prof_asym          = prof$asym,
    poisson_preferred  = pois$poisson_preferred
  )
  cat("  >> Overall:", overall, "\n\n")
  
  tibble(
    Species               = sp,
    Convergence           = conv$message,
    N_GRTS_cells          = n_grts,
    N_observations        = n_obs,
    RE_variance_GRTS      = re_var,
    YearS_estimate        = round(ye$estimate, 4),
    YearS_SE              = round(ye$se, 4),
    YearS_wald_p          = ye$p_value,
    SE_ratio              = ye$se_ratio,
    Boot_LRT_obs          = boot$obs,
    Boot_p                = boot$p,
    Boot_valid_reps       = boot$nvalid,
    LOCO_YearS_min        = loco$min,
    LOCO_YearS_max        = loco$max,
    LOCO_sign_stable      = loco$sign_stable,
    LOCO_max_rel_change   = loco$max_rel,
    Profile_CI_lo         = prof$lo,
    Profile_CI_hi         = prof$hi,
    Profile_method        = prof$method,
    Profile_excludes_zero = prof$excludes_zero,
    Profile_asymmetry     = prof$asym,
    NB_dispersion         = pois$nb_dispersion,
    AIC_NB                = pois$aic_nb,
    AIC_Poisson           = pois$aic_pois,
    Poisson_preferred     = pois$poisson_preferred,
    Overall_flag          = overall
  )
}

# -----------------------------------------------------------------------------
# SECTION 6: Main loop
# -----------------------------------------------------------------------------

run_transect_signal_diagnostics <- function(fit_list, B = 499,
                                            do_loco = TRUE, do_profile = TRUE,
                                            do_poisson = TRUE, seed = 42) {
  cat("Transect signal & robustness diagnostics for", length(fit_list), "species\n")
  cat("Bootstrap reps:", B,
      "| LOCO:", do_loco, "| Profile:", do_profile, "| Poisson cmp:", do_poisson, "\n")
  cat("Note: each species refits the model many times — expect this to run\n")
  cat("      slower than the DHARMa version, especially the bootstrap.\n\n")
  
  results <- lapply(names(fit_list), function(sp) {
    check_transect_signal(fit_list[[sp]], sp, B = B,
                          do_loco = do_loco, do_profile = do_profile,
                          do_poisson = do_poisson, seed = seed)
  })
  bind_rows(results)
}

# -----------------------------------------------------------------------------
# SECTION 7: Run and save
# -----------------------------------------------------------------------------
# If you saved the transect selection output you can rebuild best_fits without
# re-running selection:
#AK.auto.transect.fits.cov.singlet <- readRDS("Data/Analyzed/AK_transect_fits_cov_singlet.rds")
#best_fits_transect <- lapply(AK.auto.transect.fits.cov.singlet, function(x) x$best_fit)

# Bootstrap resolution: B = 499 gives p-values to ~0.002. Drop to 199 for a
# quick pass, raise to 999 for anything sitting near p = 0.05.
diag_signal <- run_transect_signal_diagnostics(
  best_fits_transect, B = 499,
  do_loco = TRUE, do_profile = TRUE, do_poisson = TRUE
)

cat("\n=== SIGNAL SUMMARY TABLE ===\n")
print(
  diag_signal %>%
    select(Species, N_GRTS_cells, N_observations,
           YearS_estimate, Boot_p,
           LOCO_sign_stable, LOCO_max_rel_change,
           Profile_excludes_zero, Profile_asymmetry,
           Poisson_preferred, Overall_flag) %>%
    as.data.frame(),
  row.names = FALSE
)

cat("\n=== FLAG COUNTS ===\n")
print(table(diag_signal$Overall_flag))

cat("\n=== FLAG LEGEND ===\n")
cat("Report                   : Bootstrap-significant, slope stable to cell removal,\n")
cat("                           profile CI excludes 0 and is roughly symmetric.\n")
cat("Report with Caveats      : Significant and directionally stable, but few cells\n")
cat("                           or moderate slope sensitivity. Note in methods.\n")
cat("Interpret Direction Only : Signal present but magnitude unreliable — wide/skewed\n")
cat("                           CI, profile spans 0, big LOCO swing, NB unidentified,\n")
cat("                           or bootstrap could not run. Report sign + significance.\n")
cat("No Detectable Trend      : Bootstrap could not reject the no-YearS null. The\n")
cat("                           honest 'data are not saying anything detectable' result.\n")
cat("Do Not Report            : Non-convergence, YearS swamped by SE, or the trend\n")
cat("                           flips sign when a single GRTS cell is dropped.\n")
cat("\nPrimary check is the bootstrap LRT. DHARMa QQ/KS is intentionally NOT used:\n")
cat("at this sample size its distributional verdict is not trustworthy.\n")

write.csv(diag_signal, "Data/Analyzed/transect_signal_diagnostics.csv", row.names = FALSE)
cat("\nSaved to: transect_signal_diagnostics.csv\n")
