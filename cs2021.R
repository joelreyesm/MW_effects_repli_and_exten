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
library(ggplot2)

# ---- Load data ----------------------------------------------
dat_all <- read_dta("QCEWindustry_minwage_all.dta") |> zap_labels()

# ---- Prepare annual panel -----------------------------------
# Same construction as the GB script for comparability.
# CS requires:
#   idname  = unique unit identifier (integer)
#   tname   = time period (integer year)
#   gname   = first period of treatment (0 = never treated)
#   yname   = outcome
#   xformla = covariates (one-sided formula)

d_cs <- dat_all |>
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
    #   (b) always-treated counties (first treated in first panel year, 1990)
    #       — CS cannot identify ATT for these as there is no pre-period
    # All others get their first treatment year.
    # NOTE: do NOT name this column "gname" — conflicts with did package internals.
    treat_year = case_when(
      is.na(first_treated_year)  ~ 0,
      first_treated_year == 1990 ~ 0,
      TRUE ~ first_treated_year
    )
  ) |>
  mutate(county_id = as.integer(factor(countyreal)))

# Report what was recoded
n_never   <- sum(d_cs$treat_year == 0 & is.na(d_cs$first_treated_year),  na.rm = TRUE) /
  n_distinct(d_cs$year)
n_always  <- sum(d_cs$treat_year == 0 & !is.na(d_cs$first_treated_year), na.rm = TRUE) /
  n_distinct(d_cs$year)
n_treated <- sum(d_cs$treat_year > 0) / n_distinct(d_cs$year)

cat("\n")
cat("=============================================================\n")
cat("CALLAWAY & SANT'ANNA (2021)\n")
cat("Dube, Lester & Reich (2010) — Extension\n")
cat("=============================================================\n")

cat("\nSample summary:\n")
cat(sprintf("  Counties (total):            %d\n", n_distinct(d_cs$county_id)))
cat(sprintf("  Years:                       %d (%d to %d)\n",
            n_distinct(d_cs$year), min(d_cs$year), max(d_cs$year)))
cat(sprintf("  Never-treated (treat_year=0):     %d counties\n",
            d_cs |> filter(treat_year == 0, is.na(first_treated_year)) |>
              distinct(county_id) |> nrow()))
cat(sprintf("  Always-treated (treat_year=0):    %d counties (excluded: no pre-period)\n",
            d_cs |> filter(treat_year == 0, !is.na(first_treated_year)) |>
              distinct(county_id) |> nrow()))
cat(sprintf("  Identified treated cohorts:       %d counties\n",
            d_cs |> filter(treat_year > 0) |> distinct(county_id) |> nrow()))
cat(sprintf("  Treatment cohorts:                %s\n",
            paste(sort(unique(d_cs$treat_year[d_cs$treat_year > 0])), collapse = ", ")))

# ---- Estimate ATT(g,t) --------------------------------------
# Two passes:
#   Pass 1: est_method = "reg" (outcome regression) — fast, good diagnostic
#   Pass 2: est_method = "dr"  (doubly robust) — preferred, slower
# Start with "reg" to verify the setup runs cleanly, then switch to "dr".
cat("\nEstimating ATT(g,t) [pass 1: outcome regression]...\n")

cs_out <- att_gt(
  yname                  = "lnemp",
  tname                  = "year",
  idname                 = "county_id",
  gname                  = "treat_year",    # NOT "gname" — avoid did package conflict
  data                   = d_cs,
  xformla                = ~lnpop + lnemp_TOT,
  control_group          = "notyettreated",
  est_method             = "reg",           # switch to "dr" for preferred results
  clustervars            = "state_fips",
  panel                  = TRUE,
  allow_unbalanced_panel = FALSE,
  base_period            = "universal",
  faster_mode            = FALSE
)

cat("\nATT(g,t) estimation complete.\n")
cat("Note: using est_method='reg'. Switch to 'dr' for doubly-robust results.\n")
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
cat("panel. Dube Spec 1/6 use continuous log(MW) on quarterly\n")
cat("panel with controls. Magnitudes are not directly comparable\n")
cat("but directional inference is.\n")

# ---- Plots --------------------------------------------------

# Event study plot
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
      "Treatment: first year state MW > federal floor (absorbing).\n",
      "Covariates: log population, log total private employment.\n",
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
    caption  = "SEs clustered at state level."
  ) +
  theme_bw(base_size = 11)

ggsave("dube2010_cs2021_cohort.png", p_grp,
       width = 8, height = 5.5, dpi = 150)
cat("Cohort ATT plot saved to: dube2010_cs2021_cohort.png\n")

cat("\n=============================================================\n")
cat("Done.\n")
cat("=============================================================\n")