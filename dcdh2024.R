# =============================================================
# Extension: de Chaisemartin & D'Haultfœuille (2024)
# did_multiplegt_dyn Applied to Dube, Lester & Reich (2010)
#
# dCDH addresses both limitations of TWFE and the gaps left
# by CS (2021) and CGBS (2025) in this application:
#
#   1. CONTINUOUS, TIME-VARYING TREATMENT: unlike CS/CGBS which
#      require binary absorbing or time-invariant dose, dCDH
#      uses the actual quarterly log(MW) which increases and
#      decreases multiple times as state and federal MWs change.
#      This gives far richer treatment variation than the 7
#      distinct onset-dose values available to CGBS.
#
#   2. NON-ABSORBING TREATMENT: the MW is not a switch — it
#      moves up and down. dCDH's switcher design compares
#      counties whose MW changes between t-1 and t against
#      counties with the same t-1 MW that do not change,
#      making no assumption about absorbing states.
#
#   3. SPATIAL FIX (Version 2 only): by pre-demeaning outcomes
#      and treatment within contiguous border county pairs,
#      we approximate Dube's pair×period fixed effects within
#      the dCDH framework. This combines:
#        - dCDH's heterogeneity-robust continuous treatment
#          estimation
#        - Dube's geographic control group selection
#      in a single analysis.
#
# TWO VERSIONS:
#   Version 1: All-county sample, no spatial fix
#              Baseline dCDH result, comparable to Dube Spec 1
#   Version 2: Spatially pre-demeaned contiguous sample
#              Closest available implementation combining all
#              three fixes (timing, continuous treatment, spatial)
#
# HOW SPATIAL PRE-DEMEANING WORKS:
#   For each contiguous county pair (i,j) at time t:
#     y_it_tilde = y_it - mean(y_it, y_jt)
#     D_it_tilde = D_it - mean(D_it, D_jt)
#   This removes the common pair×period component from both
#   outcome and treatment. Residual variation in D_it_tilde
#   reflects the within-pair MW gap — exactly Dube's identifying
#   variation. dCDH then compares counties whose within-pair
#   MW gap increases against counties whose gap is stable.
#
# LIMITATION OF VERSION 2:
#   Pre-demeaning is not identical to including pair×period FEs
#   in the dCDH estimator itself (which is not supported). The
#   approximation is valid if pair-specific trends are additive
#   and time-invariant, but if pairs have time-varying common
#   shocks not removed by simple demeaning, some bias remains.
#   This is stated explicitly in the write-up as an
#   approximation, not a full implementation of Dube's design.
#
# References:
#   de Chaisemartin, C. & D'Haultfœuille, X. (2024).
#   "Difference-in-Differences Estimators of Intertemporal
#   Treatment Effects." Review of Economics and Statistics.
#
#   R package: install.packages("DIDmultiplegtDYN")
#
# Data files required:
#   QCEWindustry_minwage_all.dta
#   QCEWindustry_minwage_contig.dta
# =============================================================

library(haven)
library(dplyr)
library(DIDmultiplegtDYN)   # install.packages("DIDmultiplegtDYN")
library(ggplot2)

# ---- Load data ----------------------------------------------
dat_all    <- read_dta("QCEWindustry_minwage_all.dta")    |> zap_labels()
dat_contig <- read_dta("QCEWindustry_minwage_contig.dta") |> zap_labels()

cat("\n")
cat("=============================================================\n")
cat("DE CHAISEMARTIN & D'HAULTFOEUILLE (2024)\n")
cat("did_multiplegt_dyn — Dube, Lester & Reich (2010)\n")
cat("=============================================================\n")

# =============================================================
# VERSION 1: ALL-COUNTY SAMPLE, NO SPATIAL FIX
# Quarterly panel, treatment = log(MW), group = county
# Directly comparable to Dube Spec 1 (TWFE) and GB result
# =============================================================

cat("\n--- VERSION 1: All-county, no spatial fix ---\n")
cat("Treatment: quarterly log(MW) | Group: county | Time: quarter\n")
cat("Comparable to Dube Spec 1\n\n")

d_all <- dat_all |>
  filter(nonmissing_rest_both == 66) |>
  select(countyreal, period, year, quarter,
         lnemp_rest_both, lnMW, lnpop, lnemp_TOT, state_fips) |>
  mutate(county_id = as.integer(factor(countyreal))) |>
  arrange(county_id, period)

cat(sprintf("All-county sample: %d counties, %d periods\n",
            n_distinct(d_all$county_id),
            n_distinct(d_all$period)))
cat(sprintf("Treatment (lnMW) range: %.3f to %.3f\n",
            min(d_all$lnMW), max(d_all$lnMW)))

# Number of effects and placebos
# With 66 quarters, use 8 post-treatment periods and 8 placebos
# to match Dube's distributed lag window
N_EFFECTS  <- 8
N_PLACEBO  <- 8

cat(sprintf("\nEstimating %d dynamic effects + %d placebo periods...\n",
            N_EFFECTS, N_PLACEBO))
cat("(Bootstrap SEs — this will take several minutes)\n\n")

mod_v1 <- did_multiplegt_dyn(
  df         = d_all,
  outcome    = "lnemp_rest_both",
  group      = "county_id",
  time       = "period",
  treatment  = "lnMW",
  effects    = N_EFFECTS,
  placebo    = N_PLACEBO,
  controls   = c("lnpop", "lnemp_TOT"),
  cluster    = "state_fips",
  breps      = 200,
  graph_off  = TRUE
)

cat("\n--- Version 1 results ---\n")
print(mod_v1)

# =============================================================
# VERSION 2: SPATIALLY PRE-DEMEANED CONTIGUOUS SAMPLE
# Approximate Dube's pair×period FEs via within-pair demeaning
# =============================================================

cat("\n\n--- VERSION 2: Spatially pre-demeaned contiguous sample ---\n")
cat("Outcome and treatment demeaned within each contiguous pair×period.\n")
cat("Residual variation = within-pair MW gap (Dube's identifying variation).\n\n")

# ---- Prepare contiguous data --------------------------------
d_ct_raw <- dat_contig |>
  filter(nonmissing_rest_both == 66) |>
  select(countyreal, pair_id, period, year, quarter,
         lnemp_rest_both, lnMW, lnpop, lnemp_TOT,
         state_fips, pair_id_period) |>
  arrange(pair_id, period, countyreal)

cat(sprintf("Contiguous sample: %d county-pair obs, %d unique counties\n",
            nrow(d_ct_raw),
            n_distinct(d_ct_raw$countyreal)))

# ---- Within-pair demeaning ----------------------------------
# For each (pair_id, period) cell, subtract the pair×period mean
# from both outcome and treatment. This removes the shared local
# economic trend that each county has in common with its neighbor.
d_ct_demeaned <- d_ct_raw |>
  group_by(pair_id, period) |>
  mutate(
    pair_mean_lnemp = mean(lnemp_rest_both, na.rm = TRUE),
    pair_mean_lnMW  = mean(lnMW,            na.rm = TRUE),
    pair_mean_lnpop = mean(lnpop,           na.rm = TRUE),
    pair_mean_lemp_TOT = mean(lnemp_TOT,    na.rm = TRUE),
    # Demeaned variables
    lnemp_dm  = lnemp_rest_both - pair_mean_lnemp,
    lnMW_dm   = lnMW            - pair_mean_lnMW,
    lnpop_dm  = lnpop           - pair_mean_lnpop,
    lempTOT_dm = lnemp_TOT      - pair_mean_lemp_TOT
  ) |>
  ungroup() |>
  # For counties in multiple pairs, average the demeaned values
  # (each county appears p times if it belongs to p pairs)
  group_by(countyreal, period) |>
  summarise(
    lnemp_dm   = mean(lnemp_dm,   na.rm = TRUE),
    lnMW_dm    = mean(lnMW_dm,    na.rm = TRUE),
    lnpop_dm   = mean(lnpop_dm,   na.rm = TRUE),
    lempTOT_dm = mean(lempTOT_dm, na.rm = TRUE),
    state_fips = first(state_fips),
    .groups    = "drop"
  ) |>
  mutate(county_id = as.integer(factor(countyreal))) |>
  arrange(county_id, period)

cat(sprintf("After demeaning: %d counties, %d periods\n",
            n_distinct(d_ct_demeaned$county_id),
            n_distinct(d_ct_demeaned$period)))

# Sanity check: within-pair MW gap variation
mw_gap_sd <- d_ct_demeaned |>
  group_by(period) |>
  summarise(sd_gap = sd(lnMW_dm, na.rm = TRUE)) |>
  pull(sd_gap) |> mean()
cat(sprintf("Average cross-period SD of within-pair MW gap: %.4f\n", mw_gap_sd))
cat("(Should be > 0; zero would mean no within-pair variation)\n\n")

cat(sprintf("Estimating %d dynamic effects + %d placebo periods...\n",
            N_EFFECTS, N_PLACEBO))
cat("(Bootstrap SEs — this will take several minutes)\n\n")

mod_v2 <- did_multiplegt_dyn(
  df         = d_ct_demeaned,
  outcome    = "lnemp_dm",
  group      = "county_id",
  time       = "period",
  treatment  = "lnMW_dm",
  effects    = N_EFFECTS,
  placebo    = N_PLACEBO,
  controls   = c("lnpop_dm", "lempTOT_dm"),
  cluster    = "state_fips",
  breps      = 200,
  graph_off  = TRUE
)

cat("\n--- Version 2 results ---\n")
print(mod_v2)

# =============================================================
# SUMMARY COMPARISON
# =============================================================

cat("\n")
cat("=============================================================\n")
cat("SUMMARY COMPARISON ACROSS ALL ESTIMATORS\n")
cat("=============================================================\n")

# Extract overall effect from dCDH output
# mod$results contains the effects table
extract_dcdh <- function(mod, label) {
  tryCatch({
    eff <- mod$results$Effects
    if (!is.null(eff) && nrow(eff) >= 1) {
      # Average effect across post-treatment periods
      avg <- mean(eff[, "Estimate"], na.rm = TRUE)
      cat(sprintf("  %-45s  avg post = %7.4f\n", label, avg))
    }
  }, error = function(e) {
    cat(sprintf("  %-45s  (could not extract)\n", label))
  })
}

cat("\n")
extract_dcdh(mod_v1, "dCDH V1: all-county, continuous lnMW")
extract_dcdh(mod_v2, "dCDH V2: contig pre-demeaned, continuous lnMW")
cat(sprintf("  %-45s  %7.4f\n", "CS (2021) overall ATT (binary, annual):", -0.0463))
cat(sprintf("  %-45s  %7.4f\n", "GB clean comparisons only:", -0.0543))
cat(sprintf("  %-45s  %7.4f\n", "GB TWFE (all comparisons):", -0.0252))
cat(sprintf("  %-45s  %7.4f\n", "Dube Spec 1 (TWFE, quarterly, w/ controls):", -0.1766))
cat(sprintf("  %-45s  %7.4f\n", "Dube Spec 6 (contig pair x period, preferred):", -0.0169))

cat("\nInterpretation guide:\n")
cat("  dCDH V1 vs Dube Spec 1: does using the switcher design\n")
cat("  (vs TWFE) change the employment estimate when we keep\n")
cat("  the all-county sample and continuous treatment?\n")
cat("\n  dCDH V2 vs Dube Spec 6: does adding heterogeneity-robust\n")
cat("  continuous treatment estimation on top of Dube's spatial\n")
cat("  fix change the near-zero employment estimate?\n")
cat("\n  dCDH V2 vs dCDH V1: how much of the remaining estimate\n")
cat("  in V1 is driven by spatial heterogeneity?\n")

# =============================================================
# PLOTS
# =============================================================

# Helper: tidy dCDH results for ggplot
tidy_dcdh <- function(mod, type = c("effects", "placebo")) {
  type <- match.arg(type)
  if (type == "effects") {
    tbl <- mod$results$Effects
    sign <- 1
  } else {
    tbl <- mod$results$Placebos
    sign <- -1   # placebo periods are indexed negatively
  }
  if (is.null(tbl)) return(NULL)
  data.frame(
    event_time = sign * seq_len(nrow(tbl)),
    estimate   = tbl[, "Estimate"],
    se         = tbl[, "SE"],
    lb         = tbl[, "LB CI"],
    ub         = tbl[, "UB CI"]
  )
}

plot_dcdh <- function(mod, title, subtitle) {
  eff <- tidy_dcdh(mod, "effects")
  plac <- tidy_dcdh(mod, "placebo")
  df_plot <- bind_rows(
    if (!is.null(plac)) mutate(plac, type = "Placebo"),
    if (!is.null(eff))  mutate(eff,  type = "Effect")
  )
  if (is.null(df_plot) || nrow(df_plot) == 0) {
    message("No data to plot for: ", title)
    return(NULL)
  }
  ggplot(df_plot, aes(x = event_time, y = estimate,
                      colour = type, fill = type)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -0.0169, linetype = "dotted",
               colour = "darkred", linewidth = 0.7) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2, colour = NA) +
    geom_point(size = 2) +
    geom_line() +
    annotate("text", x = max(df_plot$event_time) * 0.6,
             y = -0.0169 - 0.03,
             label = "Dube Spec 6 (-0.017)",
             size = 3, colour = "darkred") +
    scale_colour_manual(values = c("Effect" = "#2171b5",
                                   "Placebo" = "#74c476")) +
    scale_fill_manual(values   = c("Effect" = "#2171b5",
                                   "Placebo" = "#74c476")) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Periods relative to MW change",
      y        = "Effect on log restaurant employment",
      colour   = NULL, fill = NULL,
      caption  = paste0(
        "Blue = post-treatment effects. Green = placebo (pre-trend test).\n",
        "Treatment: continuous log(MW). Bootstrap SEs (200 reps).\n",
        "Dotted red line = Dube et al. (2010) Spec 6 preferred estimate (-0.017)."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

p_v1 <- plot_dcdh(
  mod_v1,
  "dCDH (2024): All-County Sample",
  "Continuous log(MW) | County FEs + time FEs automatic | Controls: lnpop, lnemp_TOT"
)
if (!is.null(p_v1)) {
  ggsave("dube2010_dcdh_allcounty.png", p_v1,
         width = 8, height = 5.5, dpi = 150)
  cat("\nPlot saved: dube2010_dcdh_allcounty.png\n")
}

p_v2 <- plot_dcdh(
  mod_v2,
  "dCDH (2024): Spatially Pre-Demeaned Contiguous Sample",
  "Within-pair MW gap | Approximates Dube pair×period FEs | Continuous treatment"
)
if (!is.null(p_v2)) {
  ggsave("dube2010_dcdh_spatial.png", p_v2,
         width = 8, height = 5.5, dpi = 150)
  cat("Plot saved: dube2010_dcdh_spatial.png\n")
}

cat("\n=============================================================\n")
cat("Done.\n")
cat("=============================================================\n")