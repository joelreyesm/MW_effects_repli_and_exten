# =============================================================
# Extension: Callaway, Goodman-Bacon & Sant'Anna (2024/2025)
# Continuous Treatment DiD Applied to Dube, Lester & Reich (2010)
#
# CGBS addresses a limitation of both TWFE and CS (2021): when
# treatment is continuous, a single coefficient conflates the
# extensive margin (treated vs. untreated) with the intensive
# margin (how much treatment). CGBS estimates a dose-response
# function ATT(d) — the average effect of receiving dose d
# relative to no treatment — nonparametrically via B-splines.
#
# This script estimates two target parameters:
#   ATT(d)  [target_parameter = "level"]  — effect of dose d
#            vs. untreated, aggregated across cohorts/periods
#   ACRT(d) [target_parameter = "slope"]  — marginal causal
#            response: derivative of ATT(d) w.r.t. dose
#
# And two aggregations:
#   "dose"       — ATT/ACRT as a function of the MW dose
#   "eventstudy" — dynamic ATT/ACRT averaged across doses
#
# TREATMENT VARIABLE (dname):
#   Log MW in the first year the county's state exceeded the
#   federal floor, held constant across all periods for that
#   unit. This is required by contdid: dname must be time-
#   invariant and set to its actual value (not 0) in pre-
#   treatment periods. For never/always-treated units the
#   variable is defined but does not drive identification.
#
# TIMING GROUP (gname):
#   First year state MW > federal floor (0 = never/always
#   treated, same as CS script).
#
# ============================================================
# SPATIAL LIMITATION
# ============================================================
# This script uses the all-county sample. contdid does not
# support the pair×period fixed effects that are the core of
# Dube et al.'s identification strategy. As a consequence:
#
#   - The parallel trends assumption is conditional only on
#     timing (not-yet-treated control group), not geography.
#   - Counties that only ever faced the federal MW (793 never-
#     treated counties, disproportionately Southern) enter the
#     control group and may be on systematically different
#     employment trends — exactly Dube's spatial heterogeneity
#     problem.
#   - The pre-trend violations seen in the CS (2021) results
#     for early cohorts will likely persist here for the same
#     reason.
#
# The spatial dimension is addressed in the next script
# (dCDH), which applies did_multiplegt_dyn to the pre-spatially-
# demeaned contiguous county sample, combining the geographic
# fix with heterogeneity-robust continuous treatment estimation.
# ============================================================
#
# References:
#   Callaway, B., Goodman-Bacon, A. & Sant'Anna, P.H.C. (2025).
#   "Difference-in-Differences with a Continuous Treatment."
#   arXiv: 2107.02637.
#
#   contdid R package (alpha): github.com/bcallaway11/contdid
#   Install: devtools::install_github("bcallaway11/contdid")
#
# Data file required:
#   QCEWindustry_minwage_all.dta
# =============================================================

library(haven)
library(dplyr)
library(contdid)   # devtools::install_github("bcallaway11/contdid")
library(ggplot2)

# ---- Load data ----------------------------------------------
dat_all <- read_dta("QCEWindustry_minwage_all.dta") |> zap_labels()

# ---- Prepare annual panel -----------------------------------
# Key data requirements for contdid:
#   dname : time-invariant dose, set to actual value in ALL
#           periods (not 0 in pre-treatment). For treated units
#           this is log(MW) at first treatment year.
#   gname : first year of treatment (0 = never/always treated).
#   idname: integer unit id.
#   tname : integer year.
#
# We use the same gname construction as CS (2021) script for
# direct comparability across the three estimators.

d_cgbs <- dat_all |>
  filter(nonmissing_rest_both == 66) |>
  group_by(countyreal, year) |>
  summarise(
    lnemp      = mean(lnemp_rest_both, na.rm = TRUE),
    lnMW       = mean(lnMW,            na.rm = TRUE),
    federalmin = mean(federalmin,      na.rm = TRUE),
    minwage    = mean(minwage,         na.rm = TRUE),
    state_fips = first(state_fips),
    .groups    = "drop"
  ) |>
  mutate(above_federal = as.integer(minwage > federalmin + 0.001)) |>
  arrange(countyreal, year) |>
  group_by(countyreal) |>
  mutate(
    first_treated_year = if (any(above_federal == 1))
      min(year[above_federal == 1]) else NA_real_,
    # lnMW at first treatment year (time-invariant dose)
    lnMW_at_treat = if (any(above_federal == 1))
      lnMW[year == min(year[above_federal == 1])][1] else NA_real_
  ) |>
  ungroup() |>
  mutate(
    # treat_year: 0 for never-treated and always-treated (1990).
    # Do NOT name "gname" — conflicts with contdid/did internals.
    treat_year = case_when(
      is.na(first_treated_year)  ~ 0,
      first_treated_year == 1990 ~ 0,
      TRUE ~ first_treated_year
    ),
    # dose: log MW at treatment onset, time-invariant.
    # For never/always-treated units, set to their average lnMW
    # (required to be non-missing and time-invariant by contdid;
    # the value for g=0 units does not affect identification).
    dose = case_when(
      treat_year > 0 ~ lnMW_at_treat,
      TRUE ~ mean(lnMW, na.rm = TRUE)   # placeholder for g=0 units
    )
  ) |>
  # dose must be exactly time-invariant per unit — fill within county
  group_by(countyreal) |>
  mutate(dose = first(dose[!is.na(dose)])) |>
  ungroup() |>
  mutate(county_id = as.integer(factor(countyreal)))

cat("\n")
cat("=============================================================\n")
cat("CALLAWAY, GOODMAN-BACON & SANT'ANNA (2024/2025)\n")
cat("Continuous Treatment DiD — Dube, Lester & Reich (2010)\n")
cat("=============================================================\n")

cat("\nSample summary:\n")
cat(sprintf("  Counties (total):            %d\n", n_distinct(d_cgbs$county_id)))
cat(sprintf("  Years:                       %d (%d to %d)\n",
            n_distinct(d_cgbs$year), min(d_cgbs$year), max(d_cgbs$year)))
cat(sprintf("  Never-treated (treat_year=0):%d counties\n",
            d_cgbs |> filter(treat_year == 0, is.na(first_treated_year)) |>
              distinct(county_id) |> nrow()))
cat(sprintf("  Always-treated (treat_year=0):%d counties (excluded)\n",
            d_cgbs |> filter(treat_year == 0, !is.na(first_treated_year)) |>
              distinct(county_id) |> nrow()))
cat(sprintf("  Identified treated cohorts:  %d counties\n",
            d_cgbs |> filter(treat_year > 0) |> distinct(county_id) |> nrow()))
cat(sprintf("  Treatment cohorts:           %s\n",
            paste(sort(unique(d_cgbs$treat_year[d_cgbs$treat_year > 0])),
                  collapse = ", ")))

dose_range <- d_cgbs |> filter(treat_year > 0) |> pull(dose)
cat(sprintf("  Dose (lnMW at treat onset):  min=%.3f  median=%.3f  max=%.3f\n",
            min(dose_range, na.rm=TRUE),
            median(dose_range, na.rm=TRUE),
            max(dose_range, na.rm=TRUE)))

cat("\n[SPATIAL LIMITATION] This uses the all-county sample.\n")
cat("contdid does not support pair x period FEs. The spatial\n")
cat("heterogeneity problem identified by Dube et al. and\n")
cat("visible in the CS (2021) pre-trend violations is NOT\n")
cat("addressed here. See the dCDH script for the spatial fix.\n")

# ---- Why CGBS is not feasible with this dataset -------------
#
# CGBS requires within-cohort dose variation: units in the same
# timing group must have different doses so the package can fit
# ATT(d) as a function of d within each 2x2 subset.
#
# In this dataset, the minimum wage is set at the STATE level.
# Every county in a given cohort has exactly the same dose:
#
#   Cohort 1992: 19 counties, 1 unique dose (1.6194)
#   Cohort 1993:  1 county,   1 unique dose (1.4997)
#   Cohort 1996:  2 counties, 1 unique dose (1.5197)
#   Cohort 2004: 63 counties, 1 unique dose (1.7047)
#   Cohort 2005: 101 counties,2 unique doses (1.7721, 1.7918)
#   Cohort 2006: 22 counties, 1 unique dose (1.8165)
#
# Only 7 distinct dose values across all 208 treated counties.
# When contdid tries to fit B-splines within each (g,t) cell,
# the treated units all have the same x-value — boundary knots
# cannot be established from a degenerate single-point
# distribution, causing: "Cannot set boundary knots from x."
#
# This is not a code bug or a fixable parameter choice. It is
# a structural incompatibility: CGBS is designed for settings
# where treatment intensity varies WITHIN timing groups (e.g.,
# counties receiving different pollution doses from a plant
# based on distance, or firms receiving different subsidy
# amounts). State-level minimum wages, by definition, assign
# the same dose to all counties in the same state-cohort.
#
# WHAT CGBS WOULD REQUIRE TO BE FEASIBLE HERE:
#   - County-level variation in effective MW (e.g. local MWs,
#     tip credit variation, or a county-level exposure measure)
#   - Or: define dose as the MW gap relative to the federal
#     floor at treatment onset, which still has limited cross-
#     county variation within states
#   - Or: use a different treatment variable with genuine
#     within-cohort variation (e.g. the Aladangady et al.
#     CBSA-level lock-in exposure instrument)
#
# CONCLUSION:
#   CGBS provides the right framework for this question —
#   does the employment effect vary with how much the MW
#   was raised? — but cannot be implemented with state-level
#   MW data using the onset-dose approach. The conceptual
#   contribution of CGBS to this replication is therefore
#   diagnostic: it establishes that the dose variation in
#   this dataset is insufficient for nonparametric dose-
#   response estimation, and that a richer treatment variable
#   would be needed to fully implement it.
#
#   The dCDH estimator (next script) sidesteps this by using
#   the time-varying log(MW) as the treatment rather than a
#   time-invariant onset dose, which gives far richer variation
#   and is compatible with the contiguous county spatial fix.

cat("\n")
cat("=============================================================\n")
cat("CGBS FEASIBILITY ASSESSMENT\n")
cat("=============================================================\n")
cat("\nCGBS requires within-cohort dose variation.\n")
cat("In this dataset, MW is set at the state level:\n\n")

dose_by_cohort <- d_cgbs |>
  filter(treat_year > 0) |>
  group_by(treat_year) |>
  summarise(
    n_counties   = n_distinct(county_id),
    n_doses      = n_distinct(round(dose, 4)),
    min_dose     = min(dose, na.rm = TRUE),
    max_dose     = max(dose, na.rm = TRUE),
    .groups = "drop"
  )

for (i in seq_len(nrow(dose_by_cohort))) {
  cat(sprintf("  Cohort %d: %3d counties, %d unique dose(s), range [%.4f, %.4f]\n",
              dose_by_cohort$treat_year[i],
              dose_by_cohort$n_counties[i],
              dose_by_cohort$n_doses[i],
              dose_by_cohort$min_dose[i],
              dose_by_cohort$max_dose[i]))
}

cat(sprintf("\n  Total distinct doses: %d across %d treated counties\n",
            n_distinct(round(dose_range, 4)),
            length(unique(d_cgbs$county_id[d_cgbs$treat_year > 0]))))

cat("\nConclusion: B-splines cannot be fit within each (g,t) cell\n")
cat("when all treated units share the same dose. CGBS is\n")
cat("infeasible with state-level MW data and onset-dose\n")
cat("construction. See comments above for what would be needed.\n")
cat("\nProceeding to dCDH (next script), which uses time-varying\n")
cat("log(MW) and the contiguous county spatial fix.\n")

cat("\n=============================================================\n")
cat("Done (CGBS infeasibility confirmed — see dCDH script).\n")
cat("=============================================================\n")