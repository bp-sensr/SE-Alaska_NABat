##Model Selection for stationary models in report
#load raw data after running through the NABatTrends2025 code to clean up
bat.data.5yr <- read.csv("Data/bat.data.5yr.csv")
bat.transect.data.5yr <- read.csv("Data/bat.transect.data.5yr.Aug20.csv")

#first get the means for the covariates across the sites (to know)
AK.cov.site.mean <- plyr::ddply(bat.data.5yr, c("GRTS.Cell.ID","Year"), plyr::summarize,
                             #Nightly.Mean.Temp  = mean(as.numeric(Nightly.Mean.Temp),  na.rm=TRUE),
                             Nightly.Min.Temp  = mean(Nightly.Min.Temp, na.rm=TRUE),
                             Nightly.Max.Temp  = mean(Nightly.Max.Temp, na.rm=TRUE),
                             Nightly.Mean.RH = mean(as.numeric(Nightly.Mean.RH), na.rm=TRUE),
                             Nightly.Mean.Windsp  = mean(Nightly.Mean.Windsp,  na.rm=TRUE),
                             Nightly.Mean.Precipitaton =mean(Nightly.Precipitation,  na.rm=TRUE))
AK.cov.year.mean <- plyr::ddply(AK.cov.site.mean, c("Year"), plyr::summarize,
                                #Nightly.Mean.Temp  = mean(as.numeric(Nightly.Mean.Temp),  na.rm=TRUE),
                                Nightly.Min.Temp  = mean(Nightly.Min.Temp, na.rm=TRUE),
                                Nightly.Max.Temp  = mean(Nightly.Max.Temp, na.rm=TRUE),
                                Nightly.Mean.RH = mean(as.numeric(Nightly.Mean.RH), na.rm=TRUE),
                                Nightly.Mean.Windsp  = mean(Nightly.Mean.Windsp,  na.rm=TRUE),
                                Nightly.Mean.Precipitaton =mean(Nightly.Mean.Precipitaton,  na.rm=TRUE))

# =====================================================================
# Programmatic candidate-model generator for AIC selection
#   * sweeps all combinations of the available covariates
#   * never puts two collinear variables in the same model
#   * no interaction terms
# =====================================================================

# ---- 1.1. Covariate pool ----------------------------------------------
# Variables inside the same group are treated as collinear: at most ONE
# member of a group can appear in any single model.
collinear_groups <- list(
  clutter = c("Distance.to.Clutter..m.", "Percent.Clutter"),
  water   = c("Water.Nearby", "WaterDist"),
  temp    = c("Nightly.Min.Temp", "Nightly.Max.Temp"),
  rh      = c("Nightly.Mean.RH", "Nightly.Min.RH", "Nightly.Max.RH", "Nightly.Precipitation")
)

# Variables with no collinearity constraint (each free to be in or out).
free_vars <- c("Nightly.Mean.Windsp", "jNight")

# ---- 1.2. Generator ---------------------------------------------------
make_candidates <- function(collinear_groups, free_vars,
                            min_terms = 1, max_terms = 4) {
  
  # Each "slot" offers either nothing (NA) or exactly one of its members,
  # which is what guarantees no two collinear terms co-occur.
  choice_list <- c(
    lapply(collinear_groups, function(g) c(NA_character_, g)),
    lapply(free_vars,        function(v) c(NA_character_, v))
  )
  
  grid <- do.call(expand.grid, c(choice_list, stringsAsFactors = FALSE))
  
  # Collapse each row to its non-NA terms.
  term_sets <- lapply(seq_len(nrow(grid)), function(i) {
    x <- unlist(grid[i, ], use.names = FALSE)
    x[!is.na(x)]
  })
  
  # Filter by model size.
  k <- lengths(term_sets)
  keep <- k >= min_terms & k <= max_terms
  term_sets <- term_sets[keep]
  
  # Comma-delimited strings, ordered by complexity then name.
  models <- vapply(term_sets, paste, collapse = ",", FUN.VALUE = character(1))
  models <- models[order(lengths(term_sets), models)]
  
  setNames(as.list(models), paste0("m", seq_along(models)))
}

# ---- 1.3. Build the shared candidate set ------------------------------
model_set <- make_candidates(collinear_groups, free_vars,
                             min_terms = 1, max_terms = 4)

length(model_set)     
#head(model_set)     # peek at the first few

# ---- 1.4. Assign the same set to every species ------------------------
species <- c("LACI", "LANO", "MYCA", "MYEV", "MYLU", "MYVO", "MYYU")
candidate_vars <- setNames(
  lapply(species, function(sp) model_set),
  species
)

# ── 2.  Helper: fit one glmmTMB model given a var string ────────────────────

fit_one_model <- function(cell.summary, var_string) {
  
  # Parse variables and drop rows with NAs in any covariate
  vars_vec <- if (nchar(trimws(var_string)) == 0) {
    character(0)
  } else {
    strsplit(var_string, ",")[[1]]
  }
  
  check <- unlist(sapply(strsplit(vars_vec, split = "*", fixed = TRUE), `[`, simplify = "TRUE"))
  check <- check[!is.na(check) & nchar(trimws(check)) > 0]
  
  # Strip any I(...) wrappers before trying to use as column names
  check_cols <- gsub("I\\((.*)\\)", "\\1", check)
  check_cols <- gsub("\\^.*", "", check_cols)          # remove ^2 etc.
  check_cols <- trimws(check_cols)
  check_cols <- check_cols[check_cols %in% names(cell.summary)]
  
  data_clean <- tidyr::drop_na(cell.summary, any_of(check_cols))
  
  if (nrow(data_clean) < 5) return(NULL)   # not enough data to fit
  
  # Build formula
  rhs_covs <- if (length(vars_vec) > 0) paste("+", paste(vars_vec, collapse = " + ")) else ""
  form <- paste0(
    "total.detect ~ YearS + log(quad.nights) + (1|GRTS.Cell.ID) + (1|Quadrant)",
    rhs_covs
  )
  
  tryCatch(
    glmmTMB::glmmTMB(
      formula = formula(form),
      family  = nbinom2(link = "log"),
      data    = data_clean,
      se      = TRUE,
      verbose = FALSE
    ),
    error   = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) { cat("  WARNING:", conditionMessage(w), "\n"); NULL }
  )
}


# ── 3.  Main loop: fit all candidates per species, rank by AIC ───────────────
sploc <- read.csv("Data/sploc.csv",check.names = FALSE, row.names = 1)
AK.fits.cov.singlet <- plyr::dlply(bat.data.5yr, "SpeciesGroup", function(x) {
  
  sp <- x$SpeciesGroup[1]
  cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Species:", sp, "\n")
  
  # ── Use the shared programmatic candidate set for EVERY species
  sp_candidates <- model_set
  if (is.null(sp_candidates) || length(sp_candidates) == 0) {
    cat("  model_set is empty — skipping\n")
    return(NULL)
  }
  
  # ── Filter to surveyed GRTS cells
  spgr_row  <- sploc[sp, ]
  grts_list <- names(spgr_row)[spgr_row == TRUE]
  x         <- x[x$GRTS.Cell.ID %in% grts_list, ]
  cat("  GRTS cells:", n_distinct(x$GRTS.Cell.ID), "\n")
  
  # ── Summarise to cell × year × quadrant
  #    Covariate column names are chosen to MATCH model_set exactly.
  cell.summary <- plyr::ddply(
    x, c("Year", "GRTS.Cell.ID", "Quadrant"), plyr::summarize,
    YearS                   = Year[1] - 2019,
    total.detect            = sum(SpeciesSingleton),
    quad.nights             = length(SpeciesSingleton),
    Distance.to.Clutter..m. = median(Distance.to.Clutter..m., na.rm = TRUE),  # was Clutterdist
    Percent.Clutter         = mean(Percent.Clutter,           na.rm = TRUE),
    Water.Nearby            = median(Water.Nearby,            na.rm = TRUE),
    WaterDist               = median(WaterDist,               na.rm = TRUE),   # added
    Nightly.Min.Temp        = mean(Nightly.Min.Temp,          na.rm = TRUE),
    Nightly.Max.Temp        = mean(Nightly.Max.Temp,          na.rm = TRUE),
    Nightly.Mean.RH         = mean(Nightly.Mean.RH,           na.rm = TRUE),
    Nightly.Min.RH          = mean(Nightly.Min.RH,            na.rm = TRUE),   # added
    Nightly.Max.RH          = mean(Nightly.Max.RH,            na.rm = TRUE),   # added
    Nightly.Precipitation   = mean(Nightly.Precipitation,     na.rm = TRUE),   # added
    Nightly.Mean.Windsp     = mean(Nightly.Mean.Windsp,       na.rm = TRUE),
    jNight                  = mean(jNight),
    Feature.Sampled         = unique(Feature.Sampled, na.rm = TRUE)[1],
    Habitat.Type            = unique(Habitat.Type,    na.rm = TRUE)[1]
  )
  cell.summary$countperquad <- cell.summary$total.detect / cell.summary$quad.nights
  
  # ── Restrict to complete cases across ALL candidate covariates, so every
  #    model is fit on the SAME rows (required for AIC to be comparable).
  all_cov <- unique(unlist(strsplit(unlist(sp_candidates), ",")))
  all_cov <- all_cov[all_cov %in% names(cell.summary)]
  n_before <- nrow(cell.summary)
  cell.summary <- tidyr::drop_na(cell.summary, dplyr::any_of(all_cov))
  cat("  Rows:", n_before, "->", nrow(cell.summary),
      "after complete-case filter on", length(all_cov), "covariates\n")
  
  # ── Fit every candidate model
  cat("  Fitting", length(sp_candidates), "candidate models...\n")
  
  fits <- lapply(names(sp_candidates), function(model_name) {
    cat("  →", model_name, ":", sp_candidates[[model_name]], "\n")
    fit_one_model(cell.summary, sp_candidates[[model_name]])
  })
  names(fits) <- names(sp_candidates)
  
  # ── Build AIC comparison table
  aic_table <- do.call(rbind, lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    if (is.null(fit)) {
      data.frame(model = model_name, vars = sp_candidates[[model_name]],
                 AIC = NA, converged = FALSE, stringsAsFactors = FALSE)
    } else {
      data.frame(model = model_name, vars = sp_candidates[[model_name]],
                 AIC = AIC(fit), converged = TRUE, stringsAsFactors = FALSE)
    }
  }))
  
  aic_table <- aic_table[order(aic_table$AIC), ]
  aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC, na.rm = TRUE)
  
  cat("\n  AIC comparison for", sp, ":\n")
  print(aic_table, row.names = FALSE)
  
  best_model_name <- aic_table$model[1]
  cat("\n  ✔ Best model:", best_model_name, "(AIC =", round(aic_table$AIC[1], 2), ")\n")
  
  list(
    fits      = fits,
    aic_table = aic_table,
    best_fit  = fits[[best_model_name]],
    best_vars = sp_candidates[[best_model_name]]
  )
})


# ── 4.  Summary AIC table across all species ─────────────────────────────────

all_aic <- do.call(rbind, lapply(names(AK.fits.cov.singlet), function(sp) {
  res <- AK.fits.cov.singlet[[sp]]
  if (is.null(res)) return(NULL)
  cbind(species = sp, res$aic_table)
}))

cat("\n\n═══════════ FULL AIC SUMMARY ═══════════\n")
print(all_aic, row.names = FALSE)


# ── 5.  Extract best fits only (same structure as your original output) ───────

best_fits <- lapply(AK.fits.cov.singlet, `[[`, "best_fit")

AK.cov.slope.summary.singlet <- plyr::ldply(best_fits, function(x){
  if(is.null(x$fit))return(NULL)
  summary.table <- summary(x)
  #browser()
  slope   <- as.data.frame(summary.table$coefficients$cond[2,, drop=FALSE])
  slope$intercept <- summary.table$coefficients$cond[1,,drop=FALSE]
  slope$n <- nobs(x)
  slope
})


AK.cov.slope.summary.singlet <- plyr::rename(AK.cov.slope.summary.singlet, c("Pr(>|z|)"="p.value","Std. Error"="SE",".id"="SpeciesGroup"))
# AK.cov.slope.summary.singlet$model <- paste(AK.cov.slope.summary.singlet$Var1, AK.cov.slope.summary.singlet$Var2, AK.cov.slope.summary.singlet$Var3, AK.cov.slope.summary.singlet$Var4, AK.cov.slope.summary.singlet$Var5, AK.cov.slope.summary.singlet$Var6, AK.cov.slope.summary.singlet$Var7, sep = " + ")
temptable <- AK.cov.slope.summary.singlet[,c("SpeciesGroup","Estimate","SE","p.value")]

kable(temptable, row.names=FALSE, 
      caption="Estimated trends on logarithmic scale in bat detections for semipool, fullpool and singlet models",
      col.names=c("Species","Estimate","SE","P-value"),
      digits=c(0,  2,2,4))  %>% 
  add_header_above(c(" "=1, "AutoID AK W Covariates"=3)) %>%
  column_spec(column=c(1),       width="12cm") %>%
  column_spec(column=c(2:4),       width="2cm") %>%
  kable_styling("bordered",position = "center", full_width=FALSE, latex_options = "HOLD_position")  ####%$%$


write.csv(AK.cov.slope.summary.singlet, file="C:/Users/cami/Documents/SE-Alaska_NABat/Data/Analyzed/AK.cov.estimates.singlet.csv", row.names=FALSE)
saveRDS(AK.fits.cov.singlet,"Data/Analyzed/AK.all.fits.rds")

##-------------------------Make the vars lists best on the best fit models---------------------
varlists <- do.call(rbind, lapply(names(AK.fits.cov.singlet), function(sp) {
  res <- AK.fits.cov.singlet[[sp]]
  if (is.null(res)) return(NULL)
  data.frame(
    vars = res$best_vars,
    sp   = sp,
    stringsAsFactors = FALSE
  )
}))

write.csv(varlists,"C:/Users/cami/Documents/SE-Alaska_NABat/Data/Analyzed/ModelVariables.csv")



##-----------------------------Model Selection for transect models in report------------------------------------------------------------------------------------
# ── 1.1. Covariate pool ────────────────────────────────────────────────────
# Transect surveys cover a lot of ground, so site-specific covariates
# (clutter, water distance, etc.) are excluded. Only weather / temporal
# covariates are eligible here.
#
tcollinear_groups <- list(
  temp    = c("Nightly.Min.Temp", "Nightly.Max.Temp"),
  rh      = c("Nightly.Mean.RH", "Nightly.Min.RH", "Nightly.Max.RH", "Nightly.Precipitation")
)

# Weather / temporal covariates with no collinearity constraint.
tfree_vars <- c("Nightly.Mean.Windsp", "jNight")

# ── 1.2. Generator (identical to the stationary script) ─────────────────────
# If make_candidates() is already sourced in this session, delete this copy
# and keep a single source of truth.
make_candidates <- function(tcollinear_groups, tfree_vars,
                            min_terms = 1, max_terms = 2) {
  
  # Each "slot" offers either nothing (NA) or exactly one of its members,
  # which is what guarantees no two collinear terms co-occur.
  choice_list <- c(
    lapply(tcollinear_groups, function(g) c(NA_character_, g)),
    lapply(tfree_vars,        function(v) c(NA_character_, v))
  )
  
  grid <- do.call(expand.grid, c(choice_list, stringsAsFactors = FALSE))
  
  # Collapse each row to its non-NA terms.
  term_sets <- lapply(seq_len(nrow(grid)), function(i) {
    x <- unlist(grid[i, ], use.names = FALSE)
    x[!is.na(x)]
  })
  
  # Filter by model size.
  k <- lengths(term_sets)
  keep <- k >= min_terms & k <= max_terms
  term_sets <- term_sets[keep]
  
  # Comma-delimited strings, ordered by complexity then name.
  models <- vapply(term_sets, paste, collapse = ",", FUN.VALUE = character(1))
  models <- models[order(lengths(term_sets), models)]
  
  setNames(as.list(models), paste0("m", seq_along(models)))
}

# ── 1.3. Build the shared candidate set ─────────────────────────────────────
model_set_transect <- make_candidates(tcollinear_groups,
                                      tfree_vars,
                                      min_terms = 1, max_terms = 4)

length(model_set_transect)  
#head(model_set_transect)    # peek at the first few


# ── 2. Helper: fit one transect glmmTMB model given a var string ────────────
fit_one_transect_model <- function(cell.summary, var_string) {
  
  # Parse variables and drop rows with NAs in any covariate
  vars_vec <- if (nchar(trimws(var_string)) == 0) {
    character(0)
  } else {
    strsplit(var_string, ",")[[1]]
  }
  
  check <- unlist(sapply(strsplit(vars_vec, split = "*", fixed = TRUE), `[`, simplify = "TRUE"))
  check <- check[!is.na(check) & nchar(trimws(check)) > 0]
  
  # Strip any I(...) wrappers before trying to use as column names
  check_cols <- gsub("I\\((.*)\\)", "\\1", check)
  check_cols <- gsub("\\^.*", "", check_cols)
  check_cols <- trimws(check_cols)
  check_cols <- check_cols[check_cols %in% names(cell.summary)]
  
  data_clean <- tidyr::drop_na(cell.summary, any_of(check_cols))
  
  if (nrow(data_clean) < 5) return(NULL)
  
  # Build formula — keeps transect.length and single random effect
  rhs_covs <- if (length(vars_vec) > 0) paste("+", paste(vars_vec, collapse = " + ")) else ""
  form <- paste0(
    "total.detect ~ YearS + log(transect.length) + (1|GRTS.Cell.ID)",
    rhs_covs
  )
  
  # ── Fit: record warnings without discarding the fit, catch only true errors.
  #    A benign "NA/NaN function evaluation" warning is emitted when the
  #    optimizer probes a bad point mid-search but still converges. We keep
  #    the fit and disqualify it later only if it genuinely failed.
  fit <- withCallingHandlers(
    tryCatch(
      glmmTMB::glmmTMB(
        formula = formula(form),
        family  = nbinom1(link = "log"),
        data    = data_clean,
        se      = TRUE,
        verbose = FALSE
      ),
      error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
    ),
    warning = function(w) {
      cat("  (warning, fit kept):", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }
  )
  
  # ── Keep the fit only if it actually converged with finite SEs.
  #    This is the real disqualifier for a count model: a non-positive-definite
  #    Hessian, non-finite SEs, or a nonzero convergence code.
  if (!is.null(fit)) {
    ok <- tryCatch({
      pdHess    <- isTRUE(fit$sdr$pdHess)
      ses       <- sqrt(diag(vcov(fit)$cond))
      finite_se <- all(is.finite(ses))
      conv_ok   <- fit$fit$convergence == 0
      pdHess && finite_se && conv_ok
    }, error = function(e) FALSE)
    if (!ok) {
      cat("  DROPPED: non-convergent or non-finite SEs\n")
      fit <- NULL
    }
  }
  
  fit
}


# ── 3. Main loop: fit all candidates per species, rank by AIC ───────────────

AK.auto.transect.fits.cov.singlet <- plyr::dlply(bat.transect.data.5yr, "SpeciesGroup", function(x) {
  
  sp <- x$SpeciesGroup[1]
  cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Species:", sp, "\n")
  
  # ── Use the shared programmatic candidate set for EVERY species
  sp_candidates <- model_set_transect
  if (is.null(sp_candidates) || length(sp_candidates) == 0) {
    cat("  model_set_transect is empty — skipping\n")
    return(NULL)
  }
  
  # ── Summarise to cell × year
  cell.summary <- plyr::ddply(x, c("Year", "GRTS.Cell.ID"), plyr::summarize,
                              YearS                 = Year[1] - 2015,
                              total.detect          = sum(SpeciesSingleton),
                              transect.nights       = length(SpeciesSingleton),
                              transect.length       = mean(TLength),
                              Nightly.Min.Temp      = mean(Nightly.Min.Temp,      na.rm = TRUE),
                              Nightly.Max.Temp      = mean(Nightly.Max.Temp,      na.rm = TRUE),
                              Nightly.Mean.RH       = mean(Nightly.Mean.RH,       na.rm = TRUE),
                              Nightly.Min.RH        = mean(Nightly.Min.RH,        na.rm = TRUE),
                              Nightly.Max.RH        = mean(Nightly.Max.RH,        na.rm = TRUE),
                              Nightly.Precipitation = mean(Nightly.Precipitation, na.rm = TRUE),
                              Nightly.Mean.Windsp   = mean(Nightly.Mean.Windsp,   na.rm = TRUE),
                              jNight                = mean(jNight,                na.rm = TRUE))
  
  cell.summary$GRTS.Cell.ID <- as.factor(cell.summary$GRTS.Cell.ID)
  
  # ── Restrict to complete cases across ALL candidate covariates, so every
  #    model is fit on the SAME rows (required for AIC to be comparable).
  all_cov <- unique(unlist(strsplit(unlist(sp_candidates), ",")))
  all_cov <- all_cov[all_cov %in% names(cell.summary)]
  n_before <- nrow(cell.summary)
  cell.summary <- tidyr::drop_na(cell.summary, dplyr::any_of(all_cov))
  cat("  Rows:", n_before, "->", nrow(cell.summary),
      "after complete-case filter on", length(all_cov), "covariates\n")
  
  # ── Fit every candidate model
  cat("  Fitting", length(sp_candidates), "candidate models...\n")
  
  fits <- lapply(names(sp_candidates), function(model_name) {
    cat("  →", model_name, ":", sp_candidates[[model_name]], "\n")
    fit_one_transect_model(cell.summary, sp_candidates[[model_name]])
  })
  names(fits) <- names(sp_candidates)
  
  # ── Build AIC comparison table
  aic_table <- do.call(rbind, lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    if (is.null(fit)) {
      data.frame(model = model_name, vars = sp_candidates[[model_name]],
                 AIC = NA, converged = FALSE, stringsAsFactors = FALSE)
    } else {
      data.frame(model = model_name, vars = sp_candidates[[model_name]],
                 AIC = AIC(fit), converged = TRUE, stringsAsFactors = FALSE)
    }
  }))
  
  aic_table <- aic_table[order(aic_table$AIC), ]
  aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC, na.rm = TRUE)
  
  cat("\n  AIC comparison for", sp, ":\n")
  print(aic_table, row.names = FALSE)
  
  best_model_name <- aic_table$model[1]
  cat("\n  ✔ Best model:", best_model_name, "(AIC =", round(aic_table$AIC[1], 2), ")\n")
  
  # ── Return everything
  list(
    fits      = fits,
    aic_table = aic_table,
    best_fit  = fits[[best_model_name]],
    best_vars = sp_candidates[[best_model_name]]
  )
})


# ── 4. Summary AIC table across all species ──────────────────────────────────

all_aic_transect <- do.call(rbind, lapply(names(AK.auto.transect.fits.cov.singlet), function(sp) {
  res <-AK.auto.transect.fits.cov.singlet[[sp]]
  if (is.null(res)) return(NULL)
  cbind(species = sp, res$aic_table)
}))

cat("\n\n═══════════ FULL AIC SUMMARY ═══════════\n")
print(all_aic_transect, row.names = FALSE)


# ── 5. Extract best fits only ─────────────────────────────────────────────────

best_fits_transect <- lapply(AK.auto.transect.fits.cov.singlet, `[[`, "best_fit")

AK.transect.cov.slope.summary.singlet <- plyr::ldply(best_fits_transect, function(x){
  if (is.null(x) || !inherits(x, "glmmTMB")) return(NULL)
  summary.table <- summary(x)
  slope   <- as.data.frame(summary.table$coefficients$cond[2, , drop = FALSE])
  slope$intercept <- summary.table$coefficients$cond[1, , drop = FALSE]
  slope
})

AK.transect.cov.slope.summary.singlet <- plyr::rename(
  AK.transect.cov.slope.summary.singlet,
  c("Pr(>|z|)" = "p.value", "Std. Error" = "SE", ".id" = "SpeciesGroup")
)
#AK.slope.summary
AK.transect.cov.slope.summary.singlet <- plyr::rename(AK.transect.cov.slope.summary.singlet,
                                                      c("Pr(>|z|)"="p.value","Std. Error"="SE",".id"="SpeciesGroup"))


AK.transect.cov.slope.summary.table.singlet <- AK.transect.cov.slope.summary.singlet[,c("SpeciesGroup","Estimate","SE","p.value")]
AK.transect.cov.slope.summary.table.singlet[,"p.value"] <- insight::format_p(AK.transect.cov.slope.summary.table.singlet[,"p.value"])

write.csv(AK.transect.cov.slope.summary.table.singlet, file="Data/Analyzed/AK.transect.cov.estimatesMV.csv", row.names=FALSE)


# ── 6. Save variable list for best models ────────────────────────────────────

varlists_transect <- do.call(rbind, lapply(names(AK.auto.transect.fits.cov.singlet), function(sp) {
  res <- AK.auto.transect.fits.cov.singlet[[sp]]
  if (is.null(res)) return(NULL)
  data.frame(
    vars = res$best_vars,
    sp   = sp,
    stringsAsFactors = FALSE
  )
}))

write.csv(varlists_transect, "Data/Analyzed/ModelVariables_TransectMV.csv")