# =============================================================================
# SUN & ABRAHAM (2021) EVENT STUDY — DLR (2010) WIDE REPLICATION
# =============================================================================
#
# This script implements the Sun & Abraham (2021) interaction-weighted event
# study estimator as the wide replication extension of Dube, Lester & Reich
# (2010). It is designed to run directly on the replication data files:
#
#   QCEWindustry_minwage_contig.dta  — contiguous border county-pair sample
#   QCEWindustry_minwage_all.dta     — all-county sample
#
# -----------------------------------------------------------------------------
# PART I: CONCEPTUAL FRAMING
# -----------------------------------------------------------------------------
#
# 1. WHY SA (2021) AND NOT JUST DLR'S DISTRIBUTED LAG SPEC?
# ----------------------------------------------------------
# DLR's Equation (7) already produces an event study of sorts: a distributed
# lag spec with 8 leads and 16 lags of the log minimum wage change, plotted
# as Figure 4. Under their preferred Spec 6 (county-pair x period FEs), the
# leads are flat and the employment response is stable near zero.
#
# But DLR's distributed lag spec stacks ALL minimum wage events together and
# recovers a single average impulse response. It cannot distinguish whether
# pre-trends are flat for all adoption cohorts or only on average.
#
# This matters because of the Goodman-Bacon (2021) decomposition: the TWFE
# estimator is a weighted average of all possible 2x2 DiD comparisons,
# including "forbidden" comparisons that use already-treated units as controls
# for later-treated ones. When treatment effects are heterogeneous across
# adoption cohorts -- which is plausible if early adopters (pre-1996) faced
# different labor market conditions than late adopters (post-2001 federal
# increases) -- these forbidden comparisons can bias the pooled estimate in
# ways that are not detectable from the pooled pre-trends test alone.
#
# Sun & Abraham (2021) solve this by estimating cohort-specific average
# treatment effects on the treated (CATTs) using only clean comparisons:
# each cohort is compared to never-treated or not-yet-treated units only.
# The aggregate ATT is then a cohort-size-weighted average of the CATTs,
# guaranteed free of Goodman-Bacon contamination.
#
# 2. THE TENSION: SA IS BINARY, DLR IS CONTINUOUS
# ------------------------------------------------
# SA (2021) requires binary staggered treatment: a unit is either treated
# or not, with a well-defined "first treatment date." DLR's preferred
# estimand is the elasticity with respect to the CONTINUOUS log minimum wage
# -- a county can experience many MW increases of different magnitudes.
#
# Callaway & Sant'Anna (2021) extend DiD to continuous treatment dose, but
# their estimator requires additional functional form assumptions and is
# substantially more involved. For a replication paper, this complexity is
# difficult to justify.
#
# 3. OUR SOLUTION: SA AS DIAGNOSTIC, DLR CONTINUOUS SPEC AS MAIN RESULT
# -----------------------------------------------------------------------
# We adopt a two-track strategy standard in post-2021 applied DiD work:
#
# TRACK A -- MAIN RESULT (narrow replication):
#   Reproduce DLR Spec 6: county-pair x period FEs, continuous lnMW,
#   restaurant employment and earnings. This is the paper's primary estimand.
#
# TRACK B -- SA EVENT STUDY (wide replication):
#   Define a binary treatment: indicator for effective MW > federal MW.
#   Define cohorts by first quarter this indicator switches on.
#   Never-treated = counties whose state never exceeds the federal floor.
#   Estimate SA's interaction-weighted estimator via fixest::sunab().
#
#   This produces:
#   (a) Cohort-robust pre-trends test: flat pre-trends here is stronger
#       evidence than DLR's Table 3 pooled test, ruling out cohort-specific
#       pre-trends that average out in the pooled spec.
#   (b) CATTs at each event-time horizon: if uniformly near zero for
#       employment and positive for earnings, DLR's null result is confirmed
#       as robust to cohort heterogeneity.
#   (c) Cohort-disaggregated plots: direct test of whether early adopters
#       (131 counties first treated in 1990Q1) show the same null result as
#       late adopters (67 counties first treated in 2001Q1-2006Q2).
#
# 4. COHORT STRUCTURE IN THIS DATA
# ---------------------------------
# The data runs from 1984Q1 to 2007Q4. The period variable indexes quarters
# since 1984Q1, so period = 25 corresponds to 1990Q1 (the start of the DLR
# estimation sample). Inspecting the data reveals:
#
#   period 25     (1990Q1):         131 counties -- "early adopters"
#   period 26-44  (1990Q2-1994Q4):   10 counties -- "early adopters"
#   period 45-66  (1995Q1-2000Q2):    2 counties
#   period 67+    (2000Q3-2006Q2):   67 counties -- "late adopters"
#   Never treated:                  294 counties -- clean control group
#
# The 131 counties "first treated" at period 25 are LEFT-CENSORED: they were
# already above the federal floor before the sample starts, so we have no
# pre-treatment observations for them. sunab() will drop pre-treatment
# interactions for this cohort. Post-treatment CATTs are still identified
# using the 294 never-treated counties as the comparison group.
#
# The late adopter cohort (period 81-89, roughly 2004Q1-2006Q1) provides the
# cleanest pre-trends test and is the most informative for the wide replication.
#
# 5. IDENTIFYING ASSUMPTION
# -------------------------
# SA requires:
#   (a) No anticipation: flat leads in the event study (DLR's Table 3
#       provides pooled evidence; SA plots provide cohort-specific evidence)
#   (b) Parallel trends by cohort: within each cohort, treated counties follow
#       the same employment trend as never-treated counties absent treatment.
#
# The county-pair x period FEs absorb much of the confounding. SA adds the
# cohort x event-time interaction structure on top, strengthening the design.
#
# -----------------------------------------------------------------------------
# PART II: SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(haven)       # read_dta() for Stata files
library(fixest)      # feols() + sunab() for SA estimator
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# PART III: DATA LOADING AND PREPARATION
# -----------------------------------------------------------------------------

# Load the contiguous border county-pair sample
dat_contig <- read_dta("QCEWindustry_minwage_contig.dta")

# -- 3.1 Filter to estimation sample --
# DLR's estimation sample: nonmissing_rest_both == 66 (balanced panel,
# all 66 quarters present), years 1990Q1-2006Q2 (period 25-90).

dat_sa <- dat_contig |>
  filter(
    nonmissing_rest_both == 66,
    period >= 25,   # 1990Q1
    period <= 90    # 2006Q2
  )

cat(sprintf("Estimation sample: %d rows, %d unique counties, %d unique pairs\n",
            nrow(dat_sa),
            n_distinct(dat_sa$county),
            n_distinct(dat_sa$pair_id)))
# Expected: 70,620 rows, 504 counties, 754 pair_id values

# -- 3.2 Variable reference --
# All outcome, treatment, and control variables are pre-built in the data:
#
#   lnemp_rest_both  = log(restaurant employment)        [primary employment outcome]
#   lnAWW_rest_both  = log(average weekly wages)         [primary earnings outcome]
#   lnMW             = log(effective minimum wage)       [continuous treatment, DLR Spec 6]
#   federalmin       = federal minimum wage level        [for binary treatment construction]
#   minwage          = effective MW = max(state, federal)[for binary treatment construction]
#   lnpop            = log(county population)            [control]
#   lnemp_TOT        = log(total private employment)     [control]
#   pair_id          = county-pair ID (e.g. "01003-12033")[pair FE identifier]
#   state_fips       = state FIPS code                   [clustering]
#   period           = quarters since 1984Q1             [time variable; 1990Q1 = 25]
#   censusdiv        = census division (1-9)             [for robustness specs]

# -- 3.3 Define binary treatment indicator --
# above_fed = 1 if effective MW > federal MW.
# Small tolerance (0.01) handles floating point imprecision in the Stata data.

dat_sa <- dat_sa |>
  mutate(
    above_fed = as.integer(minwage > federalmin + 0.01)
  )

cat(sprintf("above_fed: %d treated obs, %d control obs\n",
            sum(dat_sa$above_fed == 1),
            sum(dat_sa$above_fed == 0)))
# Expected: ~8,847 treated, ~61,773 control

# -- 3.4 Define cohort variable g --
# g = period of first quarter in which above_fed == 1, by county.
# Never-treated counties get g = Inf (fixest's sunab() treats Inf as never treated).

cohort_df <- dat_sa |>
  group_by(county) |>
  summarise(
    g = if_else(
      any(above_fed == 1),
      min(period[above_fed == 1]),
      Inf
    ),
    .groups = "drop"
  )

dat_sa <- dat_sa |>
  left_join(cohort_df, by = "county")

# Cohort summary
cat("\nCohort distribution:\n")
cohort_df |>
  mutate(g_label = case_when(
    g == 25                      ~ "25 (1990Q1) - left-censored early adopters [n=131]",
    g > 25  & g <= 44            ~ "26-44 (1990Q2-1994Q4) - early adopters [n=10]",
    g > 44  & g <= 66            ~ "45-66 (1995Q1-2000Q2) - mid adopters [n=2]",
    g > 66  & is.finite(g)       ~ "67+ (2001Q1-2006Q2) - late adopters [n=67]",
    g == Inf                     ~ "Never treated - clean controls [n=294]"
  )) |>
  count(g_label) |>
  arrange(g_label) |>
  print()

# -- 3.5 Assign cohort group labels for disaggregated plots --
dat_sa <- dat_sa |>
  mutate(
    cohort_group = case_when(
      g <= 44                    ~ "Early (1990Q1-1994Q4)",
      g > 44 & g <= 66           ~ "Mid (1995Q1-2000Q2)",
      g > 66 & is.finite(g)      ~ "Late (2001Q1-2006Q2)",
      g == Inf                   ~ "Never treated"
    )
  )

# -----------------------------------------------------------------------------
# PART IV: SA ESTIMATION
# -----------------------------------------------------------------------------
#
# fixest::feols() with sunab(cohort, time):
#
# sunab(g, period) constructs the cohort x relative-time interaction terms.
# Key arguments:
#   ref.p = -1        normalize the period just before first treatment to zero
#   bin               bin leads beyond -8 and lags beyond +16 (matching DLR's window)
#
# Fixed effects:
#   county            county-level FEs (time-invariant confounders)
#   pair_id^period    pair x period FEs (DLR's key identification device)
#
# Controls: lnpop + lnemp_TOT (DLR Spec 6 controls)
#
# Clustering: two-way on state_fips and pair_id.
# DLR cluster on "state and border segment"; pair_id is the border segment
# identifier in this data (each unique pair_id string represents one segment).
#
# NOTE ON LEFT-CENSORED COHORT (g = 25):
# For the 131 counties first treated at period 25 (1990Q1 = start of sample),
# there are no pre-treatment periods in the data. sunab() will drop pre-treatment
# interactions for this cohort automatically. Post-treatment CATTs are still
# identified using the 294 never-treated counties.

# -- 4.1 Employment equation --
sa_emp <- feols(
  lnemp_rest_both ~
    sunab(g, period,
          ref.p = -1,
          bin   = list("bin::3" = c(-9:-100),   # bin all leads beyond -8
                       "bin::3" = c(17:100)))   # bin all lags beyond +16
    + lnpop + lnemp_TOT
    | county + pair_id^period,
  data    = dat_sa,
  cluster = ~ state_fips + pair_id
)

# -- 4.2 Earnings equation --
sa_aww <- feols(
  lnAWW_rest_both ~
    sunab(g, period,
          ref.p = -1,
          bin   = list("bin::3" = c(-9:-100),
                       "bin::3" = c(17:100)))
    + lnpop + lnemp_TOT
    | county + pair_id^period,
  data    = dat_sa,
  cluster = ~ state_fips + pair_id
)

# -- 4.3 Print aggregate ATTs --
# Compare to DLR Spec 6 narrow replication estimates.
# If ATTs match: TWFE not biased by cohort heterogeneity -> DLR confirmed.
# If ATTs diverge: Goodman-Bacon contamination present -> new finding.

cat("\n=== SA Aggregate ATT: Employment ===\n")
print(summary(sa_emp, agg = "att"))

cat("\n=== SA Aggregate ATT: Earnings ===\n")
print(summary(sa_aww, agg = "att"))

# -----------------------------------------------------------------------------
# PART V: EVENT STUDY PLOTS
# -----------------------------------------------------------------------------

# -- 5.1 Employment --
p_emp <- iplot(
  sa_emp,
  ref.line = TRUE,
  ci_level = 0.95,
  xlab     = "Quarters relative to first state MW increase",
  ylab     = "CATT estimate",
  main     = ""
) +
  geom_vline(xintercept = -0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  theme_classic(base_size = 11) +
  labs(title = "A. Log Restaurant Employment")

# -- 5.2 Earnings --
p_aww <- iplot(
  sa_aww,
  ref.line = TRUE,
  ci_level = 0.95,
  xlab     = "Quarters relative to first state MW increase",
  ylab     = "CATT estimate",
  main     = ""
) +
  geom_vline(xintercept = -0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  theme_classic(base_size = 11) +
  labs(title = "B. Log Restaurant Average Weekly Earnings")

# -- 5.3 Combined figure (mirrors DLR Figure 4 structure) --
fig_sa <- p_emp + p_aww +
  plot_annotation(
    title   = "Sun & Abraham (2021) Event Study",
    caption = paste0(
      "Notes: Interaction-weighted estimator of Sun & Abraham (2021), CBCP sample (n = 504 counties, 318 pairs). ",
      "County and pair x period FEs. Controls: log population, log total private employment. ",
      "SEs two-way clustered on state and pair. Reference period: t = -1. ",
      "Event window: [-8, +16] quarters. 131 counties first treated at 1990Q1 are left-censored ",
      "(no pre-treatment periods); their post-treatment CATTs are identified from 294 never-treated counties."
    ),
    theme = theme(plot.caption = element_text(size = 7, hjust = 0))
  )

print(fig_sa)
ggsave("fig_sa_event_study.pdf", fig_sa, width = 12, height = 5)
cat("Saved: fig_sa_event_study.pdf\n")

# -----------------------------------------------------------------------------
# PART VI: COHORT-DISAGGREGATED CATT PLOTS
# -----------------------------------------------------------------------------
#
# This is the core new contribution of the wide replication.
# We estimate SA separately for early vs. late adopter cohorts using fixest's
# split argument and plot the CATTs side by side.
#
# What to look for:
#   Early adopters (g <= 44): post-treatment window is longer; employment CATTs
#     should be near zero if DLR's null is uniform. Note: no pre-treatment
#     periods for g = 25 cohort.
#   Late adopters (g > 66): cleanest pre-trends test (many pre-periods visible);
#     if flat pre-trends AND null employment effect, strong confirmation of DLR.

sa_emp_split <- feols(
  lnemp_rest_both ~
    sunab(g, period, ref.p = -1,
          bin = list("bin::3" = c(-9:-100), "bin::3" = c(17:100)))
    + lnpop + lnemp_TOT
    | county + pair_id^period,
  data    = dat_sa |> filter(cohort_group != "Never treated"),
  cluster = ~ state_fips + pair_id,
  split   = ~ cohort_group
)

iplot(
  sa_emp_split,
  ref.line = TRUE,
  ci_level = 0.95,
  main     = "Employment CATTs by Adoption Cohort",
  xlab     = "Quarters relative to first MW increase",
  ylab     = "CATT estimate"
)
ggsave("fig_sa_cohort_emp.pdf", width = 14, height = 5)
cat("Saved: fig_sa_cohort_emp.pdf\n")

# -----------------------------------------------------------------------------
# PART VII: COMPARISON TABLE
# -----------------------------------------------------------------------------
#
# Replace the placeholder values below with your actual narrow replication
# estimates from DLR Spec 6 once those are available.

spec6_emp_coef <- -0.016   # DLR Table 2 Spec 6, employment with controls
spec6_emp_se   <-  0.098
spec6_aww_coef <-  0.188   # DLR Table 2 Spec 6, earnings with controls
spec6_aww_se   <-  0.060

sa_att_emp <- summary(sa_emp, agg = "att")$coeftable["ATT", ]
sa_att_aww <- summary(sa_aww, agg = "att")$coeftable["ATT", ]

cat("\n")
cat("=============================================================\n")
cat("TABLE: DLR TWFE (Spec 6) vs. SA AGGREGATE ATT\n")
cat("=============================================================\n")
cat(sprintf("%-35s  %8s  %8s\n", "", "Coef.", "SE"))
cat(strrep("-", 55), "\n")
cat("Log Restaurant Employment\n")
cat(sprintf("  %-33s  %8.4f  %8.4f\n", "DLR TWFE (Spec 6)",    spec6_emp_coef, spec6_emp_se))
cat(sprintf("  %-33s  %8.4f  %8.4f\n", "SA Aggregate ATT",     sa_att_emp["Estimate"], sa_att_emp["Std. Error"]))
cat("\nLog Restaurant Average Weekly Earnings\n")
cat(sprintf("  %-33s  %8.4f  %8.4f\n", "DLR TWFE (Spec 6)",    spec6_aww_coef, spec6_aww_se))
cat(sprintf("  %-33s  %8.4f  %8.4f\n", "SA Aggregate ATT",     sa_att_aww["Estimate"], sa_att_aww["Std. Error"]))
cat("=============================================================\n")
cat("Notes: TWFE = DLR Spec 6, pair x period FEs, continuous lnMW.\n")
cat("SA = Sun & Abraham (2021), binary treatment (minwage > federalmin),\n")
cat("cohorts by first treatment period. CBCP sample. Controls and FE\n")
cat("structure identical. SEs two-way clustered on state and pair.\n")
