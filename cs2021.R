# =============================================================
# Extension: Callaway & Sant'Anna (2021) Estimator Applied to
# Dube, Lester & Reich (2010)
#
# CS (2021) addresses two problems with TWFE in staggered DiD:
#   1. Already-treated units used as controls (diagnosed by GB)
#   2. Negative and non-interpretable weights on ATT(g,t) cells
#
# We estimate group-time ATTs using not-yet-treated counties as
# the control group, with doubly-robust estimation conditioning
# on log population and log total private employment.
#
# Treatment definition: absorbing binary indicator equal to 1
# from the first year a county's state MW ever exceeded the
# federal floor. This matches the GB script exactly so results
# are directly comparable.
#
# NOTE on treatment absorbing-ness: the actual MW differential
# is episodic — 21 counties (in DE and NJ) see above_federal
# drop back to 0 when the federal floor was raised in 1996-98
# to match the state level. The treat_year variable is still
# correctly absorbing in the CS sense (first-year of adoption),
# but CS is estimating "effect of ever having had a state MW
# above federal starting in year g", not "effect of currently
# having a MW differential". This differs from Dube's continuous
# log(MW) approach and should be noted in any write-up.
#
# Control group choice: not-yet-treated (rather than never-
# treated) to avoid reintroducing Dube's spatial heterogeneity
# problem — never-treated counties are disproportionately
# Southern states that only ever faced the federal floor, which
# are on systematically different employment trajectories.
#
# References:
#   Callaway, B. & Sant'Anna, P.H.C. (2021). "Difference-in-
#   Differences with Multiple Time Periods." Journal of
#   Econometrics, 225(2), 200-230.
#
# Data file required:
#   QCEWindustry_minwage_all.dta
# =============================================================

library(haven)
library(dplyr)
library(did)      # install.packages("did") if needed
library(fastglm)
library(ggplot2)

# ---- Compatibility patches ----------------------------------
# did 2.3.0 calls fastglm via the formula API and calls predict()
# without newdata. fastglm 0.0.4 has no formula S3 method, and its
# predict() crashes on NULL newdata. Both are fixed below.

# Patch 1: fastglm.formula — convert formula to matrix and dispatch
# to the existing matrix API (fastglm.default).
.fastglm_formula_patch <- function(formula, data = NULL,
                                   family = gaussian(), ...) {
  env        <- environment(formula)
  y          <- eval(formula[[2]], envir = env)
  pred_names <- all.vars(formula[-2])
  pred_list  <- lapply(pred_names,
                       function(nm) eval(as.name(nm), envir = env))
  x          <- do.call(cbind, pred_list)
  fastglm::fastglm(x = x, y = y, family = family, intercept = FALSE, ...)
}
registerS3method("fastglm", "formula", .fastglm_formula_patch,
                 envir = asNamespace("fastglm"))

# Patch 2: predict.fastglm — fall back to model$linear.predictors
# when newdata is NULL (did 2.3.0 calls predict() without newdata to
# recover fitted propensity scores; fastglm stores them already).
.predict_fastglm_patch <- function(object, newdata = NULL,
    type = c("link", "response"), se.fit = FALSE,
    dispersion = NULL, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    if (se.fit) stop("se.fit not supported by fastglm")
    dims <- dim(newdata)
    if (is.null(dims)) { newdata <- as.matrix(newdata); dims <- dim(newdata) }
    beta <- object$coefficients
    if (dims[2] != length(beta)) stop("newdata does not match model")
    eta <- drop(newdata %*% beta)
  }
  if (type == "response") eta <- family(object)$linkinv(eta)
  eta
}
registerS3method("predict", "fastglm", .predict_fastglm_patch,
                 envir = asNamespace("fastglm"))

# ---- Load data ----------------------------------------------
dat_all <- read_dta("QCEWindustry_minwage_all.dta") |> zap_labels()

# ---- Prepare annual panel -----------------------------------
# FIX 1: The raw QCEW panel runs 1990Q1-2006Q2 (66 quarters).
# The year 2006 contains only Q1-Q2, so its annual average uses
# 2 quarters vs. 4 for every other year. Including 2006 creates
# a systematic measurement shift in the last period that would
# contaminate ATT(g, 2006) estimates. We drop 2006 and restrict
# the analysis to 1990-2005 (16 full calendar years).

d_annual <- dat_all |>
  filter(nonmissing_rest_both == 66) |>
  group_by(countyreal, year) |>
  summarise(
    lnemp      = mean(lnemp_rest_both, na.rm = TRUE),
    lnpop      = mean(lnpop,           na.rm = TRUE),
    lnemp_TOT  = mean(lnemp_TOT,       na.rm = TRUE),
    federalmin = mean(federalmin,      na.rm = TRUE),
    minwage    = mean(minwage,         na.rm = TRUE),
    state_fips = first(state_fips),
    .groups    = "drop"
  ) |>
  filter(year <= 2005)   # drop incomplete 2006

# FIX 3: CS (2021) requires covariates to be pre-determined
# (pre-treatment). lnpop and lnemp_TOT are time-varying, so we
# pin them to their 1990 baseline values, making them fixed per
# county and unambiguously pre-treatment for every cohort.
# (The did package does internally extract pre-treatment covariate
# values when panel=TRUE, but using a fixed baseline is cleaner
# and eliminates any ambiguity about which period is used.)
baseline_covs <- d_annual |>
  filter(year == 1990) |>
  select(countyreal, lnpop_base = lnpop, lnemp_TOT_base = lnemp_TOT)

d_cs <- d_annual |>
  left_join(baseline_covs, by = "countyreal") |>
  mutate(above_federal = as.integer(minwage > federalmin + 0.001)) |>
  arrange(countyreal, year) |>
  group_by(countyreal) |>
  mutate(
    first_treated_year = if (any(above_federal == 1))
      min(year[above_federal == 1]) else NA_real_
  ) |>
  ungroup() |>
  mutate(
    # treat_year = 0 for:
    #   (a) never-treated counties (first_treated_year is NA)
    #   (b) always-treated counties (first treated in 1990, no pre-period)
    #   (c) counties first treated in 2006: panel ends in 2005 after the
    #       fix above, so there is no post-treatment period to estimate
    #       ATT(g, t) for this cohort.  [FIX 5]
    #   (d) cohorts g=1993 (1 county) and g=1996 (2 counties): the DR
    #       estimator fits a propensity-score logit within each (g,t) cell.
    #       With only 1-2 treated observations the logit degenerates and
    #       DRDID throws "argument 'y' is missing". The did package itself
    #       warns about these groups. We drop them from identified cohorts;
    #       they remain in the data as part of the not-yet-treated control
    #       pool for other cohorts.
    treat_year = case_when(
      is.na(first_treated_year)              ~ 0,
      first_treated_year == 1990             ~ 0,
      first_treated_year >= 2006             ~ 0,   # no post-treatment period
      first_treated_year %in% c(1993, 1996)  ~ 0,   # too few treated units for DR
      TRUE                                   ~ first_treated_year
    )
  ) |>
  mutate(county_id = as.integer(factor(countyreal)))

# Report what was recoded
n_never   <- d_cs |> filter(treat_year == 0,  is.na(first_treated_year))  |>
  distinct(county_id) |> nrow()
n_always  <- d_cs |> filter(treat_year == 0, !is.na(first_treated_year),
                              first_treated_year == 1990) |>
  distinct(county_id) |> nrow()
n_no_post <- d_cs |> filter(treat_year == 0, !is.na(first_treated_year),
                              first_treated_year >= 2006) |>
  distinct(county_id) |> nrow()
n_treated <- d_cs |> filter(treat_year > 0) |> distinct(county_id) |> nrow()

cat("\n")
cat("=============================================================\n")
cat("CALLAWAY & SANT'ANNA (2021)\n")
cat("Dube, Lester & Reich (2010) — Extension\n")
cat("=============================================================\n")

cat("\nSample summary:\n")
cat(sprintf("  Counties (total):                    %d\n", n_distinct(d_cs$county_id)))
cat(sprintf("  Years:                               %d (%d to %d)\n",
            n_distinct(d_cs$year), min(d_cs$year), max(d_cs$year)))
cat(sprintf("  Never-treated (treat_year=0):        %d counties\n", n_never))
cat(sprintf("  Always-treated (treat_year=0):       %d counties (excluded: no pre-period)\n", n_always))
cat(sprintf("  No post-period (treat_year=0):       %d counties (excluded: first treated 2006)\n", n_no_post))
cat(sprintf("  Identified treated cohorts:          %d counties\n", n_treated))
cat(sprintf("  Treatment cohorts:                   %s\n",
            paste(sort(unique(d_cs$treat_year[d_cs$treat_year > 0])), collapse = ", ")))

# Warn about thin cohorts  [FIX 6]
cat("\nNOTE — thin cohorts:\n")
thin <- d_cs |> filter(treat_year > 0) |>
  distinct(county_id, treat_year) |>
  count(treat_year) |>
  filter(n < 10)
if (nrow(thin) > 0) {
  for (i in seq_len(nrow(thin))) {
    cat(sprintf("  Cohort g=%d has only %d county/counties — estimates will be noisy.\n",
                thin$treat_year[i], thin$n[i]))
  }
}

# ---- Estimate ATT(g,t) [doubly robust] ----------------------
# FIX 2: use est_method = "dr" (doubly robust) throughout.
# "reg" (outcome regression only) was a diagnostic pass only.
#
# FIX 7: base_period = "universal" normalises all pre-treatment
# relative periods to a single common base period rather than
# the period immediately before treatment (base_period =
# "varying"). "universal" is more conservative for pre-trend
# testing but shifts the event-study normalisation slightly.
# We keep "universal" here for clean pre-trend assessment;
# switch to "varying" if the goal is a canonical event-study
# centered exactly at g-1.

cat("\nEstimating ATT(g,t) [doubly robust, not-yet-treated controls]...\n")

cs_out <- att_gt(
  yname                  = "lnemp",
  tname                  = "year",
  idname                 = "county_id",
  gname                  = "treat_year",   # FIX 8: column holding first treatment year (0 = never/always treated)
  data                   = d_cs,
  xformla                = ~ lnpop_base + lnemp_TOT_base,   # FIX 3: fixed 1990 baseline covariates
  control_group          = "notyettreated",
  est_method             = "dr",           # FIX 2: doubly robust (preferred)
  clustervars            = "state_fips",
  panel                  = TRUE,
  allow_unbalanced_panel = FALSE,
  base_period            = "universal",
  faster_mode            = FALSE
)

cat("\nATT(g,t) estimation complete.\n")
print(summary(cs_out))

# ---- Aggregations -------------------------------------------

# 1. Simple aggregation: overall ATT
cat("\n--- Overall ATT (simple aggregation) ---\n")
agg_simple <- aggte(cs_out, type = "simple")
print(summary(agg_simple))

# 2. Dynamic aggregation: event-study plot
cat("\n--- Dynamic (event-study) aggregation ---\n")
agg_dyn <- aggte(cs_out, type = "dynamic")
print(summary(agg_dyn))

# 3. Calendar time aggregation
cat("\n--- Calendar time aggregation ---\n")
agg_cal <- aggte(cs_out, type = "calendar")
print(summary(agg_cal))

# 4. Group aggregation (ATT by treatment cohort)
cat("\n--- Group (cohort) aggregation ---\n")
agg_grp <- aggte(cs_out, type = "group")
print(summary(agg_grp))

# ---- Reference point for comparison -------------------------
cat("\n--- Summary comparison ---\n")
cat(sprintf("  CS (2021) overall ATT:               %6.4f  (SE: %6.4f)\n",
            agg_simple$overall.att, agg_simple$overall.se))
cat(sprintf("  GB clean comparisons only:          %6.4f\n", -0.0543))
cat(sprintf("  GB TWFE (all comparisons):          %6.4f\n", -0.0252))
cat(sprintf("  Dube et al. Spec 1 (TWFE w/ ctrls): %6.4f\n", -0.1766))
cat(sprintf("  Dube et al. Spec 6 (preferred):     %6.4f\n", -0.0169))
cat("\nNote: CS and GB use binary absorbing treatment on an annual\n")
cat("panel (1990-2005). Dube Spec 1/6 use continuous log(MW) on\n")
cat("quarterly panel with controls. Magnitudes are not directly\n")
cat("comparable but directional inference is.\n")

# ---- Plots --------------------------------------------------

# Event study plot
# FIX 2: subtitle now correctly reflects est_method = "dr"
p_es <- ggdid(agg_dyn) +
  geom_hline(yintercept = -0.0169, linetype = "dotted",
             colour = "darkred", linewidth = 0.8) +
  annotate("text", x = max(agg_dyn$egt) * 0.6,
           y = -0.0169 - 0.02,
           label = "Dube et al. Spec 6 (-0.017)",
           size = 3, colour = "darkred") +
  labs(
    title    = "CS (2021): Event Study — Effect of State MW Adoption",
    subtitle = "Log restaurant employment | Not-yet-treated control group | DR estimator",
    x        = "Years relative to first treatment",
    y        = "ATT",
    caption  = paste0(
      "Treatment: first year state MW > federal floor (absorbing). Panel: 1990-2005.\n",
      "Covariates: 1990 baseline log population and log total private employment.\n",
      "SEs clustered at state level. Dotted line = Dube et al. Spec 6."
    )
  ) +
  theme_bw(base_size = 11)

ggsave("dube2010_cs2021_eventstudy.png", p_es,
       width = 8, height = 5.5, dpi = 150)
cat("\nEvent study plot saved to: dube2010_cs2021_eventstudy.png\n")

# Group ATT plot
p_grp <- ggdid(agg_grp) +
  labs(
    title    = "CS (2021): ATT by Treatment Cohort",
    subtitle = "Each cohort = first year state MW exceeded federal floor",
    x        = "Treatment cohort (first year above federal floor)",
    y        = "ATT",
    caption  = "SEs clustered at state level. Panel: 1990-2005."
  ) +
  theme_bw(base_size = 11)

ggsave("dube2010_cs2021_cohort.png", p_grp,
       width = 8, height = 5.5, dpi = 150)
cat("Cohort ATT plot saved to: dube2010_cs2021_cohort.png\n")

cat("\n=============================================================\n")
cat("Done.\n")
cat("=============================================================\n")
