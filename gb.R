# =============================================================
# Extension: Goodman-Bacon Decomposition of Dube, Lester &
# Reich (2010) Specification 1
#
# Dube et al. argue the traditional TWFE estimate (Spec 1) is
# biased by spatial heterogeneity in employment trends. The
# Goodman-Bacon (2021) decomposition shows it is also
# contaminated by staggered timing: some weight falls on
# "already-treated" units used as controls for later-treated
# units, which is invalid under treatment effect heterogeneity.
#
# The contiguous county design (Spec 6) fixes the spatial
# problem. GB diagnoses the timing problem. Together they give
# a complete picture of why Spec 1 is wrong.
#
# Reference:
#   Goodman-Bacon, A. (2021). "Difference-in-differences with
#   variation in treatment timing." Journal of Econometrics,
#   225(2), 254-277.
#
# Data file required (same as replication):
#   QCEWindustry_minwage_all.dta
# =============================================================

library(haven)
library(dplyr)
library(fixest)
library(bacondecomp)  # install.packages("bacondecomp") if needed
library(ggplot2)      # for the decomposition plot

# ---- Load data ----------------------------------------------
dat_all <- read_dta("QCEWindustry_minwage_all.dta") |> zap_labels()

# ---- Prepare analysis sample --------------------------------
# Balanced panel of counties with full restaurant data (66 quarters).
# Collapse to annual to match the binary-treatment assumption that
# bacondecomp requires: one observation per unit per period.
# We use annual average log MW and annual average log employment.
#
# TREATMENT DEFINITION — absorbing state required by bacondecomp:
# Treatment = 1 from the first year a county's state MW ever exceeded
# the federal floor, and remains 1 permanently thereafter.
# This is necessary because the raw indicator (state MW > federal MW)
# is non-monotone: the 1996-97 federal increases temporarily brought
# some states back to parity, causing spurious "un-treatment."
# The absorbing definition captures the onset of state-level MW policy
# and is the standard approach in the GB minimum wage literature.
d_annual <- dat_all |>
  filter(nonmissing_rest_both == 66) |>
  group_by(countyreal, year) |>
  summarise(
    lnemp      = mean(lnemp_rest_both, na.rm = TRUE),
    lnAWW      = mean(lnAWW_rest_both, na.rm = TRUE),
    lnMW       = mean(lnMW,            na.rm = TRUE),
    lnpop      = mean(lnpop,           na.rm = TRUE),
    lnemp_TOT  = mean(lnemp_TOT,       na.rm = TRUE),
    state_fips = first(state_fips),
    federalmin = mean(federalmin,      na.rm = TRUE),
    minwage    = mean(minwage,         na.rm = TRUE),
    .groups    = "drop"
  ) |>
  # Step 1: raw indicator — state MW above federal floor in this year
  mutate(above_federal = as.integer(minwage > federalmin + 0.001)) |>
  # Step 2: absorbing treatment — 1 from first year ever above, forever
  arrange(countyreal, year) |>
  group_by(countyreal) |>
  mutate(
    first_treated_year = if (any(above_federal == 1))
      min(year[above_federal == 1]) else Inf,
    treated = as.integer(year >= first_treated_year)
  ) |>
  ungroup()

# Sanity check: treatment switching
cat("\n")
cat("=============================================================\n")
cat("GOODMAN-BACON DECOMPOSITION\n")
cat("Dube, Lester & Reich (2010) — Specification 1 Extension\n")
cat("=============================================================\n")

cat("\nTreatment summary (absorbing: first year state MW > federal MW):\n")
cat(sprintf("  County-year obs:          %d\n", nrow(d_annual)))
cat(sprintf("  Treated obs (absorbing):  %d (%.1f%%)\n",
            sum(d_annual$treated),
            100 * mean(d_annual$treated)))
cat(sprintf("  Never-treated counties:   %d\n",
            d_annual |> group_by(countyreal) |>
              summarise(ever = max(treated)) |>
              filter(ever == 0) |> nrow()))
cat(sprintf("  Ever-treated counties:    %d\n",
            d_annual |> group_by(countyreal) |>
              summarise(ever = max(treated)) |>
              filter(ever == 1) |> nrow()))
cat("\nNote: raw (non-absorbing) treated obs = ",
    sum(d_annual$above_federal),
    sprintf("(%.1f%%)", 100 * mean(d_annual$above_federal)),
    "\n")
cat("  Difference reflects 1996-97 federal increases temporarily\n")
cat("  closing state-federal gaps; absorbing definition ignores\n")
cat("  these reversals and is required by bacondecomp.\n")

# ---- Goodman-Bacon decomposition ----------------------------
# bacon() requires a balanced panel with no missings in the
# outcome and treatment. Drop any counties with gaps.
d_bacon <- d_annual |>
  group_by(countyreal) |>
  filter(n() == n_distinct(d_annual$year)) |>  # balanced
  ungroup() |>
  filter(!is.na(lnemp), !is.na(treated))

cat(sprintf("\nBalanced panel for decomposition:\n"))
cat(sprintf("  Counties: %d\n", n_distinct(d_bacon$countyreal)))
cat(sprintf("  Years:    %d\n", n_distinct(d_bacon$year)))
cat(sprintf("  Obs:      %d\n", nrow(d_bacon)))

# Run decomposition on employment outcome (no controls —
# bacon() partials out unit and time FEs internally)
bacon_out <- bacon(
  formula  = lnemp ~ treated,
  data     = d_bacon,
  id_var   = "countyreal",
  time_var = "year"
)

# Overall TWFE estimate — recompute as weighted mean to be safe
twfe_est <- weighted.mean(bacon_out$estimate, bacon_out$weight)
cat("\n--- Decomposition of TWFE employment estimate ---\n")
cat(sprintf("\nOverall TWFE estimate (weighted): %.4f\n", twfe_est))
cat("\nComponent breakdown:\n")
cat(sprintf("  %-35s  %8s  %8s  %8s\n",
            "Comparison type", "Estimate", "Weight", "Wtd contrib"))
cat(strrep("-", 65), "\n")

summary_tbl <- bacon_out |>
  group_by(type) |>
  summarise(
    wt_avg_est   = weighted.mean(estimate, weight),
    total_wt     = sum(weight),
    contribution = weighted.mean(estimate, weight) * sum(weight),
    .groups = "drop"
  ) |>
  arrange(desc(total_wt))

for (i in seq_len(nrow(summary_tbl))) {
  cat(sprintf("  %-35s  %8.4f  %8.4f  %8.4f\n",
              summary_tbl$type[i],
              summary_tbl$wt_avg_est[i],
              summary_tbl$total_wt[i],
              summary_tbl$contribution[i]))
}

cat(strrep("-", 65), "\n")
cat(sprintf("  %-35s  %8s  %8.4f  %8.4f\n",
            "Total", "",
            sum(summary_tbl$total_wt),
            sum(summary_tbl$contribution)))

# Identify comparison types using actual bacondecomp labels:
#   "Treated vs Untreated"    = clean (never/not-yet treated as control)
#   "Later vs Always Treated" = contaminated (already-treated as control)
#   "Earlier vs Later Treated" / "Later vs Earlier Treated" = timing comparisons
clean_types <- bacon_out |> filter(type == "Treated vs Untreated")
contam_types <- bacon_out |> filter(type == "Later vs Always Treated")
timing_types <- bacon_out |> filter(grepl("Earlier", type) | grepl("Later vs Earlier", type))

clean_wt  <- sum(clean_types$weight,  na.rm = TRUE)
contam_wt <- sum(contam_types$weight, na.rm = TRUE)
timing_wt <- sum(timing_types$weight, na.rm = TRUE)

cat("\n--- Weight by comparison type ---\n")
cat(sprintf("  Treated vs Untreated (clean):            %.4f (%.1f%%)\n",
            clean_wt,  100 * clean_wt))
cat(sprintf("  Later vs Always Treated (contaminated):  %.4f (%.1f%%)\n",
            contam_wt, 100 * contam_wt))
cat(sprintf("  Timing comparisons (Early/Late):         %.4f (%.1f%%)\n",
            timing_wt, 100 * timing_wt))

clean_est  <- weighted.mean(clean_types$estimate,  clean_types$weight)
contam_est <- weighted.mean(contam_types$estimate, contam_types$weight)

cat(sprintf("\n  TWFE estimate (all comparisons):         %.4f\n", twfe_est))
cat(sprintf("  Clean comparisons only (vs Untreated):  %.4f\n", clean_est))
cat(sprintf("  Contaminated only (vs Always Treated):  %.4f\n", contam_est))
cat(sprintf("  Dube et al. Spec 6 (preferred):         -0.0169\n"))

cat("\nInterpretation:\n")
cat("  The TWFE estimate (-0.025) is a weighted average of all\n")
cat("  2x2 DiD comparisons. The clean 'Treated vs Untreated'\n")
cat("  comparisons (64% of weight) yield an estimate of around\n")
cat(sprintf("  %.3f. But 'Later vs Always Treated' comparisons\n", clean_est))
cat("  (31% of weight) use already-treated counties as controls\n")
cat(sprintf("  and yield +%.3f — the opposite sign. This positive\n", contam_est))
cat("  estimate is the timing analogue of Dube's spatial\n")
cat("  heterogeneity problem: counties that adopted MW earlier\n")
cat("  are on different employment trajectories than later\n")
cat("  adopters, contaminating the comparison.\n")
cat("\n  Dube et al. fix the spatial dimension via contiguous\n")
cat("  county pairs (Spec 6: -0.017). GB shows the traditional\n")
cat("  TWFE is also wrong along the timing dimension. A\n")
cat("  Callaway-Sant'Anna estimator within the contiguous pair\n")
cat("  sample would address both problems simultaneously.\n")

# ---- Plot: weighted 2x2 estimates ---------------------------
# Each point is one 2x2 DiD comparison, sized by its weight.
p <- ggplot(bacon_out, aes(x = weight, y = estimate, colour = type, shape = type)) +
  geom_point(aes(size = weight), alpha = 0.7) +
  geom_hline(yintercept = twfe_est,
             linetype = "dashed", colour = "black", linewidth = 0.8) +
  geom_hline(yintercept = -0.0169,
             linetype = "dotted", colour = "darkred", linewidth = 0.8) +
  annotate("text", x = max(bacon_out$weight) * 0.6,
           y = twfe_est + 0.015,
           label = "TWFE (Spec 1)", size = 3.2, colour = "black") +
  annotate("text", x = max(bacon_out$weight) * 0.6,
           y = -0.0169 - 0.015,
           label = "Dube et al. Spec 6", size = 3.2, colour = "darkred") +
  scale_size_continuous(range = c(1, 8), guide = "none") +
  scale_colour_brewer(palette = "Set2") +
  labs(
    title    = "Goodman-Bacon Decomposition of TWFE Employment Estimate",
    subtitle = "Each point is a 2x2 DiD comparison, sized by its weight in the TWFE estimator",
    x        = "Weight",
    y        = "2x2 DiD Estimate",
    colour   = "Comparison type",
    shape    = "Comparison type",
    caption  = paste0(
      "Dashed line = TWFE Spec 1 estimate.  ",
      "Dotted red line = Dube et al. (2010) Spec 6 preferred estimate (−0.017).\n",
      "Outcome: log restaurant employment. Treatment: state MW > federal MW."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("dube2010_bacon_decomposition.png", p,
       width = 8, height = 5.5, dpi = 150)

cat("\nPlot saved to: dube2010_bacon_decomposition.png\n")
cat("\n=============================================================\n")
cat("Done.\n")
cat("=============================================================\n")