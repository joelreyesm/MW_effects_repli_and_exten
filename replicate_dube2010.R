# =============================================================
# Replication: Dube, Lester, Reich (2010)
# "Minimum Wages and Low-Wage Jobs: Evidence from
#  Cross-State Comparisons Across Time"
# Replication of Tables 1, 2, and 3
#
# Data files (must be in ./minwage_dubeRESTAT/):
#   QCEWindustry_minwage_all.dta
#   QCEWindustry_minwage_contig.dta
# =============================================================

library(haven)
library(dplyr)
library(fixest)

BASE <- "minwage_dubeRESTAT"

# ---- Load and clean data ------------------------------------
dat_all    <- read_dta(file.path(BASE, "QCEWindustry_minwage_all.dta"))    |> zap_labels()
dat_contig <- read_dta(file.path(BASE, "QCEWindustry_minwage_contig.dta")) |> zap_labels()


# =============================================================
# TABLE 1: Descriptive Statistics
# ("descriptives 11-6-09.do")
#
# Two panels: all-county sample and county-pair sample.
# Filter: nonmissing_rest_both == 66 (balanced restaurant panel).
# For industry variables (ACFSRETAIL, RETAIL, MFG, TOT), the
# filter is the industry-specific nonmissing indicator.
# =============================================================

cat("\n")
cat("=============================================================\n")
cat("TABLE 1: DESCRIPTIVE STATISTICS\n")
cat("=============================================================\n")

print_desc <- function(df, panel_label, has_pairs = FALSE) {
  cat(sprintf("\n--- %s ---\n", panel_label))

  d <- df |> filter(nonmissing_rest_both == 66)

  # -- Main county-level variables --
  main_vars <- c("countypop2000", "cntypopdens", "cntyarea",
                 "emp_rest_both", "AWW_rest_both", "minwage")
  cat(sprintf("\n%-22s  %13s  %13s\n", "Variable", "Mean", "Std. Dev."))
  cat(strrep("-", 52), "\n")
  for (v in main_vars) {
    cat(sprintf("%-22s  %13.3f  %13.3f\n", v,
                mean(d[[v]], na.rm = TRUE),
                sd(d[[v]],   na.rm = TRUE)))
  }

  # -- Industry employment and wages --
  # Note: empACFSRETAIL/AWWACFSRETAIL are only in the county-pair dataset;
  # they are absent from the all-county dataset. Skip silently if missing.
  cat(sprintf("\n%-14s  %12s  %12s  %12s  %12s\n",
              "Industry", "Mean Emp", "SD Emp", "Mean AWW", "SD AWW"))
  cat(strrep("-", 66), "\n")
  for (k in c("ACFSRETAIL", "RETAIL", "MFG", "TOT")) {
    nm_col  <- paste0("nonmissing_", k)
    emp_col <- paste0("emp", k)
    aww_col <- paste0("AWW", k)
    # Skip if the required columns do not exist in this dataset
    if (!emp_col %in% names(df) || !aww_col %in% names(df)) {
      cat(sprintf("%-14s  %12s  %12s  %12s  %12s  (not available)\n",
                  k, "", "", "", ""))
      next
    }
    dk <- if (nm_col %in% names(df)) df |> filter(.data[[nm_col]] == 66) else d
    cat(sprintf("%-14s  %12.3f  %12.3f  %12.3f  %12.3f\n", k,
                mean(dk[[emp_col]], na.rm = TRUE),
                sd(dk[[emp_col]],   na.rm = TRUE),
                mean(dk[[aww_col]], na.rm = TRUE),
                sd(dk[[aww_col]],   na.rm = TRUE)))
  }

  # -- Sample sizes --
  n_cty <- d |> distinct(county) |> nrow()
  cat(sprintf("\nNumber of counties:     %d\n", n_cty))

  if (has_pairs) {
    # Each county appears twice in the pairs dataset -> divide by 2
    n_pairs <- df |> filter(nonmissing_both_pair == 66) |>
      distinct(pair_id) |> nrow()
    cat(sprintf("Number of county pairs: %d\n", n_pairs))
  }

  # -- Minimum-wage events per county --
  # The Stata code uses pre-computed cumulative event counts per county
  # (steventspercounty / fedeventspercounty), then divides by 2 for the
  # pairs dataset because each county appears twice.
  divisor <- if (has_pairs) 2 else 1
  if ("event_type" %in% names(d)) {
    fed_ev <- d |> filter(event_type == 1) |>
      count(county) |>
      summarise(m = mean(n / divisor)) |> pull(m)
    sta_ev <- d |> filter(event_type == 2) |>
      count(county) |>
      summarise(m = mean(n / divisor)) |> pull(m)
    cat(sprintf("Mean federal events / county: %.3f\n", fed_ev))
    cat(sprintf("Mean state   events / county: %.3f\n", sta_ev))
  }
}

print_desc(dat_all,    "All-County Sample",  has_pairs = FALSE)
print_desc(dat_contig, "County-Pair Sample", has_pairs = TRUE)


# =============================================================
# TABLE 2: Minimum Wage Effects on Earnings and Employment
# ("main analysis 11-6-09.do")
#
# Six specifications (rows), two control sets (c2 / c3):
#
#   Spec 1: All counties,            county + period FE
#   Spec 2: All counties,            county + period×CBMSA FE
#   Spec 3: All counties,            county + period×CensusDivision FE
#   Spec 4: All counties,            county + period×CensusDivision FE
#                                    + state-specific linear trends
#   Spec 5: Contiguous county pairs, county + period FE
#   Spec 6: Contiguous county pairs, county + pair×period FE
#
# Control sets:
#   c2  AWW: no extra controls  (just lnMW)
#   c2  EMP: lnpop
#   c3  AWW: lnAWW_TOT
#   c3  EMP: lnpop + lnemp_TOT
#
# SE clustering:
#   Specs 1-4: clustered by state_fips
#   Specs 5-6: two-way (Cameron-Gelbach-Miller) by state_fips × bordersegment
# =============================================================

cat("\n")
cat("=============================================================\n")
cat("TABLE 2: MINIMUM WAGE EFFECTS ON EARNINGS AND EMPLOYMENT\n")
cat("=============================================================\n")
cat("Coefficient on log(minimum wage) | SEs in parentheses\n")
cat("c2 = population control; c3 = pop + total emp/wage control\n")

# ---- Prepare all-county analysis dataset --------------------
d_all <- dat_all |>
  filter(nonmissing_rest_both == 66) |>
  mutate(
    # Interaction FE groupings (matching Stata's egen group())
    absorb_2     = as.integer(interaction(period, cbmsa,      drop = TRUE)),
    absorb_3     = as.integer(interaction(period, censusdiv,  drop = TRUE)),
    # Factor for state-specific linear time trends (spec 4)
    state_fips_f = factor(state_fips)
  )

# ---- Helper: pull MW coefficient and SE from a feols model --
mw_coef <- function(mod) {
  b  <- coef(mod)["lnMW"]
  se <- se(mod)["lnMW"]
  c(est = unname(b), se = unname(se))
}

# ---- Print one specification row ----------------------------
print_spec_row <- function(label, aww_c2, aww_c3, emp_c2, emp_c3, n) {
  cat(sprintf("\n  [Spec %s]  N = %d\n", label, as.integer(n)))
  cat(sprintf("    AWW  c2: %8.4f  (%6.4f)   c3: %8.4f  (%6.4f)\n",
              aww_c2[1], aww_c2[2], aww_c3[1], aww_c3[2]))
  cat(sprintf("    EMP  c2: %8.4f  (%6.4f)   c3: %8.4f  (%6.4f)\n",
              emp_c2[1], emp_c2[2], emp_c3[1], emp_c3[2]))
}

# ---- Specs 1–3: all counties, state-clustered SEs -----------
spec_info <- list(
  list(num = "1", label = "All counties, period FE",
       data = d_all, absorb = "period"),
  list(num = "2", label = "All counties, period x CBMSA FE",
       data = d_all, absorb = "absorb_2"),
  list(num = "3", label = "All counties, period x CensDiv FE",
       data = d_all, absorb = "absorb_3")
)

for (sp in spec_info) {
  dd  <- sp$data
  ab  <- sp$absorb

  # Formulas — NOTE: AWW c2 has NO population control (matches Stata)
  f_aww_c2 <- as.formula(paste0("lnAWW_rest_both ~ lnMW             | county + ", ab))
  f_aww_c3 <- as.formula(paste0("lnAWW_rest_both ~ lnMW + lnAWW_TOT | county + ", ab))
  f_emp_c2 <- as.formula(paste0("lnemp_rest_both ~ lnMW + lnpop     | county + ", ab))
  f_emp_c3 <- as.formula(paste0("lnemp_rest_both ~ lnMW + lnpop + lnemp_TOT | county + ", ab))

  m_aww_c2 <- feols(f_aww_c2, data = dd, cluster = ~state_fips)
  m_aww_c3 <- feols(f_aww_c3, data = dd, cluster = ~state_fips)
  m_emp_c2 <- feols(f_emp_c2, data = dd, cluster = ~state_fips)
  m_emp_c3 <- feols(f_emp_c3, data = dd, cluster = ~state_fips)

  print_spec_row(sp$num,
                 mw_coef(m_aww_c2), mw_coef(m_aww_c3),
                 mw_coef(m_emp_c2), mw_coef(m_emp_c3),
                 nobs(m_aww_c2))
}

# ---- Spec 4: Census division FE + state linear time trends --
# Stata: xi: xtreg ... i.state_fips*period ...
# This adds state-specific linear trends: one slope on `period` per state.
# County FEs already absorb state-level intercepts, so only the
# state x period interaction (the slope) matters.

f_aww_c2_s4 <- lnAWW_rest_both ~ lnMW             + state_fips_f:period | county + absorb_3
f_aww_c3_s4 <- lnAWW_rest_both ~ lnMW + lnAWW_TOT + state_fips_f:period | county + absorb_3
f_emp_c2_s4 <- lnemp_rest_both ~ lnMW + lnpop     + state_fips_f:period | county + absorb_3
f_emp_c3_s4 <- lnemp_rest_both ~ lnMW + lnpop + lnemp_TOT + state_fips_f:period | county + absorb_3

m_aww_c2_s4 <- feols(f_aww_c2_s4, data = d_all, cluster = ~state_fips)
m_aww_c3_s4 <- feols(f_aww_c3_s4, data = d_all, cluster = ~state_fips)
m_emp_c2_s4 <- feols(f_emp_c2_s4, data = d_all, cluster = ~state_fips)
m_emp_c3_s4 <- feols(f_emp_c3_s4, data = d_all, cluster = ~state_fips)

print_spec_row("4: All counties, CensDivFE + State Trends",
               mw_coef(m_aww_c2_s4), mw_coef(m_aww_c3_s4),
               mw_coef(m_emp_c2_s4), mw_coef(m_emp_c3_s4),
               nobs(m_aww_c2_s4))

# ---- Prepare contiguous-county analysis dataset -------------
# pair_id format: "SSCCC-SSCCC" where SS = state FIPS (2 chars), CCC = county (3 chars)
# R substr uses (string, start, stop); Stata uses (string, start, length).
# Stata: substr(pair_id, 7, 2) = chars 7-8.  R: substr(pair_id, 7, 8).
d_ct <- dat_contig |>
  filter(nonmissing_rest_both == 66) |>
  mutate(
    state_a       = as.integer(substr(pair_id, 1, 2)),
    state_b       = as.integer(substr(pair_id, 7, 8)),
    st_min        = pmin(state_a, state_b),
    st_max        = pmax(state_a, state_b),
    # Border segment = unique state-pair (for two-way clustering)
    bordersegment = as.integer(interaction(st_min, st_max, drop = TRUE)),
    # pair x period FE (spec 6)
    pair_period   = as.integer(interaction(pair_id, period,  drop = TRUE))
  )

# ---- Specs 5–6: contiguous pairs, two-way clustered SEs -----
# Two-way clustering: state_fips + bordersegment
# (matches Stata's cluster2areg fcluster(state_fips) tcluster(bordersegment))

for (s in 5:6) {
  ab  <- if (s == 5) "period" else "pair_period"
  lbl <- if (s == 5) "5: Contig. pairs, period FE" else
                     "6: Contig. pairs, pair x period FE"

  f_aww_c2 <- as.formula(paste0("lnAWW_rest_both ~ lnMW             | county + ", ab))
  f_aww_c3 <- as.formula(paste0("lnAWW_rest_both ~ lnMW + lnAWW_TOT | county + ", ab))
  f_emp_c2 <- as.formula(paste0("lnemp_rest_both ~ lnMW + lnpop     | county + ", ab))
  f_emp_c3 <- as.formula(paste0("lnemp_rest_both ~ lnMW + lnpop + lnemp_TOT | county + ", ab))

  clust <- ~state_fips + bordersegment   # Cameron-Gelbach-Miller two-way

  m_aww_c2 <- feols(f_aww_c2, data = d_ct, cluster = clust)
  m_aww_c3 <- feols(f_aww_c3, data = d_ct, cluster = clust)
  m_emp_c2 <- feols(f_emp_c2, data = d_ct, cluster = clust)
  m_emp_c3 <- feols(f_emp_c3, data = d_ct, cluster = clust)

  print_spec_row(lbl,
                 mw_coef(m_aww_c2), mw_coef(m_aww_c3),
                 mw_coef(m_emp_c2), mw_coef(m_emp_c3),
                 nobs(m_aww_c2))
}


# =============================================================
# TABLE 3: Pre-Existing Trends in Employment/Earnings and
#          Validity of Controls
# ("PreTrends_nov_6_2009.do")
#
# Strategy: include LEADS of log(MW) as regressors.
#   lnMW_f4  = lnMW 4 quarters (1 year) ahead
#   lnMW_f12 = lnMW 12 quarters (3 years) ahead
#
# The regression is:
#   depvar ~ lnMW + lnMW_f4 + lnMW_f12 + lnpop | county + absorb_X
#
# Reported linear combinations (matching Stata's lincom):
#   (a) lnMW_f12 alone            → 3-year-ahead coefficient
#   (b) lnMW_f12 + lnMW_f4        → cumulative 3-year pre-trend
#   (c) lnMW_f4 alone             → 1-year-ahead coefficient
#
# Sample period for regression: 1990 Q1 – 2006 Q2.
# NOTE: leads are generated on the FULL nonmissing panel first,
#       then the sample-period restriction is applied for the
#       regression. This matches the Stata sequence exactly.
#
# Three specifications:
#   Spec 1: All counties,        county + period FE,       cluster state
#   Spec 2: Cross-state metro,   county + period×CBMSA FE, cluster state
#   Spec 3: Contiguous pairs,    county + pair×period FE,  two-way cluster
# =============================================================

cat("\n")
cat("=============================================================\n")
cat("TABLE 3: PRE-EXISTING TRENDS (LEADS OF LOG MINIMUM WAGE)\n")
cat("=============================================================\n")
cat("Dep vars: lnemp_rest_both, lnemp_TOT, lnAWW_rest_both, lnAWW_TOT\n")
cat("Reported: (a) coeff on F12 alone; (b) lincom F12+F4; (c) coeff on F4 alone\n")

# ---- Generate leads on the full nonmissing panel ------------
# IMPORTANT: leads must be computed before applying the sampleperiod
# filter so that observations near the end of sampleperiod can borrow
# lnMW values from periods outside sampleperiod (e.g. 2007-2009).

add_leads_all <- function(df) {
  df |>
    filter(nonmissing_rest_both == 66) |>
    arrange(countyreal, period) |>
    group_by(countyreal) |>
    mutate(
      lnMW_f4  = dplyr::lead(lnMW,  4),
      lnMW_f12 = dplyr::lead(lnMW, 12)
    ) |>
    ungroup()
}

add_leads_contig <- function(df) {
  df |>
    filter(nonmissing_rest_both == 66, !is.na(pair_id_county)) |>
    arrange(pair_id_county, period) |>
    group_by(pair_id_county) |>
    mutate(
      lnMW_f4  = dplyr::lead(lnMW,  4),
      lnMW_f12 = dplyr::lead(lnMW, 12)
    ) |>
    ungroup()
}

# Generate leads, then restrict to sample period for regressions
samp_all_full <- add_leads_all(dat_all)

samp_all <- samp_all_full |>
  filter(year >= 1990,
         !(year > 2006 | (year == 2006 & quarter >= 3))) |>
  mutate(
    absorb_2 = as.integer(interaction(period, cbmsa, drop = TRUE))
  )

samp_contig_full <- add_leads_contig(dat_contig)

samp_contig <- samp_contig_full |>
  filter(year >= 1990,
         !(year > 2006 | (year == 2006 & quarter >= 3))) |>
  mutate(
    state_a       = as.integer(substr(pair_id, 1, 2)),
    state_b       = as.integer(substr(pair_id, 7, 8)),
    st_min        = pmin(state_a, state_b),
    st_max        = pmax(state_a, state_b),
    bordersegment = as.integer(interaction(st_min, st_max, drop = TRUE)),
    pair_period   = as.integer(interaction(pair_id, period, drop = TRUE))
  )

# ---- Helper: run one pre-trends regression ------------------
# Returns list: (a) F12 alone, (b) F12+F4, (c) F4 alone; each with est + se
run_pretrend_reg <- function(df, depvar, absorb_var, clust_fml) {
  fml <- as.formula(
    paste0(depvar, " ~ lnMW + lnMW_f4 + lnMW_f12 + lnpop | county + ", absorb_var)
  )

  mod <- tryCatch(
    feols(fml, data = df, vcov = clust_fml),
    error = function(e) NULL
  )
  if (is.null(mod)) return(list(a = c(NA, NA), b = c(NA, NA), c = c(NA, NA), N = NA))

  co <- coef(mod)
  vc <- vcov(mod)

  get_lincom <- function(vars) {
    b  <- sum(co[vars])
    se <- sqrt(sum(vc[vars, vars]))   # delta-method for linear combination
    c(est = b, se = se)
  }

  list(
    a = get_lincom("lnMW_f12"),            # F12 alone
    b = get_lincom(c("lnMW_f12", "lnMW_f4")),  # F12 + F4 combined
    c = get_lincom("lnMW_f4"),             # F4 alone
    N = nobs(mod)
  )
}

# ---- Helper: print one Table 3 block ------------------------
depvars_t3 <- c("lnemp_rest_both", "lnemp_TOT", "lnAWW_rest_both", "lnAWW_TOT")
col_lbls   <- c("Emp-Rest", "Emp-TOT", "AWW-Rest", "AWW-TOT")

print_t3_block <- function(df, spec_label, absorb_var, clust_fml) {
  cat(sprintf("\n  Specification: %s\n", spec_label))

  # Run all 4 dep-var regressions
  res <- setNames(
    lapply(depvars_t3, run_pretrend_reg,
           df = df, absorb_var = absorb_var, clust_fml = clust_fml),
    depvars_t3
  )

  N_val <- res[[1]]$N

  # Header
  cat(sprintf("  %-10s", ""))
  for (lbl in col_lbls) cat(sprintf("  %16s", lbl))
  cat(sprintf("  %6s\n", "N"))

  # Three rows: (a), (b), (c)
  row_labels <- c("(a) F12 alone", "(b) F12+F4  ", "(c) F4  alone")
  row_keys   <- c("a", "b", "c")

  for (rk in seq_along(row_keys)) {
    cat(sprintf("  %-13s", row_labels[rk]))
    for (dv in depvars_t3) {
      r <- res[[dv]][[row_keys[rk]]]
      if (is.na(r[1])) {
        cat(sprintf("  %7s (%6s)", "NA", "NA"))
      } else {
        cat(sprintf("  %7.4f (%6.4f)", r[1], r[2]))
      }
    }
    if (rk == 1) cat(sprintf("  %6d", as.integer(N_val))) else cat(sprintf("  %6s", ""))
    cat("\n")
  }
}

# ---- Spec 1: All counties, period FE ------------------------
print_t3_block(samp_all,
               "1: All counties, period FE",
               "period",
               ~state_fips)

# ---- Spec 2: Cross-state metro, period × CBMSA FE -----------
# (sample_flag == 2 selects cross-state MSA counties)
print_t3_block(samp_all |> filter(sample_flag == 2),
               "2: Cross-state metro, period x CBMSA FE",
               "absorb_2",
               ~state_fips)

# ---- Spec 3: Contiguous county pairs, pair × period FE ------
print_t3_block(samp_contig,
               "3: Contiguous pairs, pair x period FE",
               "pair_period",
               ~state_fips + bordersegment)

cat("\n")
cat("=============================================================\n")
cat("Replication complete.\n")
cat("=============================================================\n")
